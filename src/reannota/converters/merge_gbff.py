"""Merge annotation data into GenBank files."""

import logging

import pandas as pd
from Bio import SeqIO
from goatools.obo_parser import GODag

from reannota.parsers import parse_pseudogff_to_dict

logger = logging.getLogger("ReAnnota")

DESCRIPTION_LENGTH_MAX = 50

# Pseudogene policy for Bakta-hypothetical loci that Pseudofinder flags AND
# that eggNOG/InterPro could name:
#   True  -> conservative: flag the locus as a pseudogene, keep /product as
#            "hypothetical protein" and put the candidate name in /note.
#   False -> keep the eggNOG/InterPro name as /product AND flag the locus as
#            a pseudogene.
# Either way the Pseudofinder check now runs on Bakta's ORIGINAL product, so a
# renamed locus can no longer escape the pseudogene flag.
PSEUDO_CANDIDATES_KEEP_HYPOTHETICAL = True

# GO root terms (biological_process, molecular_function, cellular_component)
# carry no information and are never added.
ROOT_GO_TERMS = {"GO:0008150", "GO:0003674", "GO:0005575"}


def _parse_gos(raw_str):
    if str(raw_str).lower() == "nan" or not raw_str:
        return []
    return [g.strip() for g in str(raw_str).replace('"', '').split(',') if 'GO:' in g]


def _canonical_go(go_id, godag):
    """Map alternative GO IDs to their primary ID (unknown IDs are returned unchanged)."""
    term = godag.get(go_id) if godag else None
    return term.id if term is not None else go_id


def _most_specific(go_ids, godag):
    """Return the terms of a set that are not an ancestor of another term in the set."""
    ancestors = set()
    for go_id in go_ids:
        term = godag.get(go_id)
        if term is not None:
            ancestors |= term.get_all_parents()
    return {g for g in go_ids if g not in ancestors}


def resolve_go_terms(egg_go_str, ipr_go_str, ipr_desc, godag, bakta_gos=()):
    """
    Implements the Maximal Specificity Algorithm.

    EggNOG-mapper reports GO terms fully propagated (every ancestor of each
    annotated term), so the candidate set is first collapsed to its most
    specific terms. Then:
      * InterPro terms are tested against EggNOG terms (exact match / parent:
        discard; child: keep; unrelated: keep + bridge note);
      * the union of EggNOG + InterPro terms is collapsed together with the
        GO terms Bakta already provided, so no added term is an ancestor of
        another added term or of a Bakta term, and no term Bakta already
        has is added again;
      * GO root terms are never added.
    Bakta's own GO terms are never removed.

    Returns: (list_of_final_go_terms, list_of_bridge_notes)
    """
    egg_list = _parse_gos(egg_go_str)
    ipr_list = _parse_gos(ipr_go_str)
    bridge_notes = []

    # If no DAG is loaded, fallback to safe addition without notes
    # (redundancy cannot be resolved without the ontology)
    if not godag:
        return list((set(egg_list) | set(ipr_list)) - ROOT_GO_TERMS - set(bakta_gos)), bridge_notes

    egg_list = sorted({_canonical_go(g, godag) for g in egg_list})
    bakta_set = {_canonical_go(g, godag) for g in bakta_gos}

    # Collapse the propagated EggNOG list to its most specific terms.
    # The InterPro comparison below still uses the full EggNOG list (roots
    # included), so the duplicate / parent / child / unrelated classification
    # (and therefore the bridge notes) behaves as before.
    final_gos = _most_specific(set(egg_list) - ROOT_GO_TERMS, godag)
    unrelated_ipr = set()

    # Evaluate every InterPro term against the EggNOG baseline
    for ipr_go in ipr_list:
        ipr_term = godag.get(ipr_go)
        if not ipr_term:
            continue
        ipr_go = ipr_term.id  # primary ID
        if ipr_go in ROOT_GO_TERMS:
            continue

        is_duplicate = False
        is_parent = False
        is_child = False

        for egg_go in egg_list:
            egg_term = godag.get(egg_go)
            if not egg_term:
                continue

            # 1. Exact Match
            if ipr_go == egg_go:
                is_duplicate = True
                break
            # 2. InterPro is Parent (EggNOG is more specific) -> We don't need it
            elif ipr_go in egg_term.get_all_parents():
                is_parent = True
                break
            # 3. InterPro is Child (InterPro is more specific) -> Keep it!
            elif egg_go in ipr_term.get_all_parents():
                is_child = True
                break

        # Apply the logic rules
        if is_duplicate or is_parent:
            continue  # Discard InterPro term
        
        elif is_child:
            final_gos.add(ipr_go)  # Add safely, no note needed (it's just higher res)
            
        else:
            # 4. Completely Unrelated (Scenario 4) -> Swiss Army Knife domain
            final_gos.add(ipr_go)
            unrelated_ipr.add(ipr_go)

    # Collapse EggNOG + InterPro + Bakta together; keep only what is new
    final_gos = _most_specific(final_gos | bakta_set, godag) - bakta_set

    # Bridge note only if an unrelated InterPro term actually survived
    if unrelated_ipr & final_gos and ipr_desc and ipr_desc.lower() != "nan":
        # Explain WHY this unrelated term is here
        bridge_notes.append(f"InterPro detected additional domain/function: {ipr_desc}")

    return sorted(final_gos), bridge_notes

# NOTE: Added 'go_obo_file' as a parameter here with a default fallback
def merge_csv_to_gbff(egg_file, interpro_file, gbff_in, gbff_out, pseudofile=None, go_obo_file="go-basic.obo"):
    """
    Merge functional annotations from EggNOG and InterPro CSV files into a GenBank (GBFF) file.
    """
    logger.debug(f"Merging annotations from EggNOG and InterPro into GBFF: {gbff_in}")

    # --- LOAD GO DATABASE DYNAMICALLY HERE ---
    godag = None
    if go_obo_file:
        logger.info(f"Loading GO database ({go_obo_file})...")
        try:
            godag = GODag(go_obo_file)
        except FileNotFoundError:
            logger.error(f"Could not find {go_obo_file}! GO terms will merge via fallback. Download it with: wget http://purl.obolibrary.org/obo/go/go-basic.obo")

    # Read EggNOG and InterPro CSV files as pandas DataFrames
    if egg_file is not None:
        logger.debug(f"Reading EggNOG annotations from: {egg_file}")
        egg_df = pd.read_csv(egg_file, sep="\t")
        logger.debug(f"Loaded {len(egg_df)} EggNOG annotation entries")
    else:
        egg_df = pd.DataFrame()  # empty DataFrame if missing

    if interpro_file is not None:
        logger.debug(f"Reading InterPro annotations from: {interpro_file}")
        ipr_df = pd.read_csv(interpro_file, sep="\t")
        logger.debug(f"Loaded {len(ipr_df)} InterPro annotation entries")
    else:
        ipr_df = pd.DataFrame()  # empty DataFrame if missing

    if pseudofile is not None :
        logger.debug(f"Reading Pseudogenes from: {pseudofile}")
        pseudo_dict = parse_pseudogff_to_dict(pseudofile)
        logger.debug(f"Loaded {len(pseudo_dict)} Pseudogenes")

    # Convert InterPro DataFrame to a dictionary keyed by Query_ID
    ipr_dictionary = {row["Query_ID"]: row.to_dict() for _, row in ipr_df.iterrows()}

    # Convert EggNOG DataFrame to a dictionary keyed by Query_ID
    egg_dictionary = {row["Query_ID"]: row.to_dict() for _, row in egg_df.iterrows()}

    gbff_read = list(SeqIO.parse(gbff_in, "genbank"))
    logger.debug(f"Loaded {len(gbff_read)} sequences from GBFF file")

    # Counters for statistics
    cds_count = 0
    ipr_annotated = 0
    egg_annotated = 0
    pseudogenes = 0
    stats_genes_affected = {
        "eggnog_product": 0,
        "interpro_product": 0,
        "eggnog_gene": 0,
        "interpro_gene": 0,
        "eggnog_go_terms": 0,
        "interpro_go_terms": 0,
        "db_xrefs": 0,
        "notes": 0,
        "pseudogene": 0
    }

    # Loop through each sequence in the GenBank file
    for sequence in gbff_read:
        for feature in sequence.features:
            if feature.type == "CDS":
                cds_count += 1
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]

                # Fetch data safely (empty dict if not found)
                egg_data = egg_dictionary.get(locus_tag, {})
                ipr_data = ipr_dictionary.get(locus_tag, {})

                # Update counters
                if egg_data:
                    egg_annotated += 1
                if ipr_data:
                    ipr_annotated += 1

                # Ensure fields exist
                feature.qualifiers.setdefault("note", [])
                feature.qualifiers.setdefault("db_xref", [])
                feature.qualifiers.setdefault("gene", [])

                #-----TRACKING STATS FOR AFFECTED GENES-----
                eggnog_product_changed = False
                interpro_product_changed = False
                eggnog_gene_added = False
                interpro_gene_added = False
                pseudogene_added = False
                
                initial_note_count = len(feature.qualifiers["note"])
                initial_go_count = sum(1 for x in feature.qualifiers["db_xref"] if x.startswith("GO:"))
                initial_other_xref = len(feature.qualifiers["db_xref"]) - initial_go_count

                # --- STEP 0: Pseudogene candidacy, decided on Bakta's ORIGINAL product ---
                # (Previously this was checked after renaming, so any hypothetical
                # locus that eggNOG/InterPro could name escaped the pseudogene flag.)
                bakta_product = feature.qualifiers.get("product", [""])[0]
                pseudo_candidate = (
                    pseudofile is not None
                    and locus_tag in pseudo_dict
                    and bakta_product == "hypothetical protein"
                    and not feature.qualifiers.get("pseudogene", [])
                )

                def assign_product(name):
                    """Set /product, or divert the name to /note for kept-hypothetical pseudogene candidates."""
                    if pseudo_candidate and PSEUDO_CANDIDATES_KEEP_HYPOTHETICAL:
                        feature.qualifiers["note"].append(f"Pseudogene candidate; similar to: {name}")
                        return False
                    feature.qualifiers["product"] = [name]
                    return True

                # --- STEP 1: Product Name & Gene Name ---
                egg_desc = str(egg_data.get("Description", "")).strip()
                ipr_desc = str(ipr_data.get("Description", "")).strip()

                junk_names = ["nan", "hypothetical protein", "", "uncharacterized protein"]

                # Assign Product Name
                if egg_desc and egg_desc.lower() not in junk_names:
                    if len(egg_desc) > DESCRIPTION_LENGTH_MAX:
                        feature.qualifiers["note"].append(egg_desc)
                    elif bakta_product == "hypothetical protein":
                        eggnog_product_changed = assign_product(egg_desc)
                elif ipr_desc and ipr_desc.lower() not in junk_names:
                    if bakta_product == "hypothetical protein":
                        interpro_product_changed = assign_product(ipr_desc)

                # Assign Gene Name (EggNOG priority)
                egg_gene = str(egg_data.get("Gene_name", "")).replace('"', '').strip()
                ipr_gene = str(ipr_data.get("Gene_name", "")).replace('"', '').strip()
                
                if egg_gene and egg_gene.lower() != "nan" and not feature.qualifiers["gene"]:
                    feature.qualifiers["gene"].append(egg_gene)
                    eggnog_gene_added = True
                elif ipr_gene and ipr_gene.lower() != "nan" and not feature.qualifiers["gene"]:
                    feature.qualifiers["gene"].append(ipr_gene)
                    interpro_gene_added = True

                # --- STEP 2: The Smart GO Term Merge ---
                egg_gos = str(egg_data.get("GOs", ""))
                ipr_gos = str(ipr_data.get("GO_terms", ""))
                
                # Parse strings into sets to check origin later (primary GO IDs)
                egg_go_set = {_canonical_go(g, godag) for g in _parse_gos(egg_gos)}
                ipr_go_set = {_canonical_go(g, godag) for g in _parse_gos(ipr_gos)}

                bakta_gos = [x for x in feature.qualifiers["db_xref"] if x.startswith("GO:")]
                final_gos, bridge_notes = resolve_go_terms(egg_gos, ipr_gos, ipr_desc, godag, bakta_gos)
                
                eggnog_gos_added = 0
                interpro_gos_added = 0

                for go in final_gos:
                    if go not in feature.qualifiers["db_xref"]:
                        feature.qualifiers["db_xref"].append(go)
                        # Identify source of the newly added GO term
                        if go in egg_go_set:
                            eggnog_gos_added += 1
                        elif go in ipr_go_set:
                            interpro_gos_added += 1
                        
                for note in bridge_notes:
                    feature.qualifiers["note"].append(note)

                # --- STEP 3: Add the remaining specific database cross-references ---
                db_xrefs = feature.qualifiers["db_xref"]

                if egg_data:
                    cog_ref = str(egg_data.get("COG_ref", ""))
                    cog_cat = str(egg_data.get("COG_category", ""))
                    keggs = str(egg_data.get("KEGG_ko", ""))
                    ecs = str(egg_data.get("EC", ""))
                    pfams = str(egg_data.get("PFAM", ""))

                    if cog_ref and cog_ref != "nan" and not any(x.startswith("COG:") for x in db_xrefs):
                        feature.qualifiers["db_xref"].append(cog_ref)
                    if cog_cat and cog_cat != "nan" and not any(x.startswith("COG_cat:") for x in db_xrefs):
                        feature.qualifiers["db_xref"].append(f"COG_cat:{cog_cat}")
                    # KOs and ECs: one db_xref per identifier (a protein can carry
                    # several). Added only when Bakta has none of that type, as before.
                    # removeprefix() also tolerates hit tables written by the old
                    # parser ("EC:x.x.x.x", "ko:Kxxxxx").
                    if keggs and keggs != "nan" and not any(x.startswith("KEGG:") for x in db_xrefs):
                        for ko in dict.fromkeys(k.strip().removeprefix("ko:") for k in keggs.split(",")):
                            if ko:
                                feature.qualifiers["db_xref"].append(f"KEGG:{ko}")
                    if ecs and ecs != "nan" and not any(x.startswith("EC:") for x in db_xrefs):
                        for ec in dict.fromkeys(e.strip().removeprefix("EC:") for e in ecs.split(",")):
                            if ec:
                                feature.qualifiers["db_xref"].append(f"EC:{ec}")
                    if pfams and pfams != "nan":
                        feature.qualifiers["note"].append(f"Eggnog PFAM comment:{pfams}")

                if ipr_data:
                    interpro_terms = str(ipr_data.get("Interpro_terms", "")).replace('"', '')
                    if interpro_terms and interpro_terms != "nan":
                        # Safely split by comma and add each unique IPR ID
                        for ipr in interpro_terms.split(","):
                            ipr = ipr.strip()
                            if ipr and not any(x.startswith(ipr) for x in db_xrefs):
                                feature.qualifiers["db_xref"].append(ipr)

                #--------------Add Pseudogenes-------------------
                if pseudofile is not None and locus_tag in pseudo_dict:
                    pseudogenes += 1

                if pseudo_candidate:
                    values = pseudo_dict[locus_tag]
                    raw_attr = values.get("attributes", "").replace("note=", "")
                    description_pseudo = raw_attr.split(';')[0].strip()
                    description_pseudo = description_pseudo.replace("%25", "%").strip()
                    description_pseudo = " ".join(description_pseudo.split())

                    feature.qualifiers.setdefault("pseudogene", [])
                    feature.qualifiers["pseudogene"].append(description_pseudo)
                    pseudogene_added = True

                # ----------------- TRACKING EVALUATION & LOGGING -----------------
                final_note_count = len(feature.qualifiers["note"])
                final_go_count = sum(1 for x in feature.qualifiers["db_xref"] if x.startswith("GO:"))
                final_other_xref = len(feature.qualifiers["db_xref"]) - final_go_count

                notes_added = final_note_count - initial_note_count
                xrefs_added = final_other_xref - initial_other_xref

                if any([eggnog_product_changed, interpro_product_changed, eggnog_gene_added, interpro_gene_added, pseudogene_added, notes_added > 0, eggnog_gos_added > 0, interpro_gos_added > 0, xrefs_added > 0]):
                    changes_msg = []
                    if eggnog_product_changed: 
                        changes_msg.append("EggNOG Product updated")
                        stats_genes_affected["eggnog_product"] += 1
                    if interpro_product_changed: 
                        changes_msg.append("InterPro Product updated")
                        stats_genes_affected["interpro_product"] += 1
                    if eggnog_gene_added: 
                        changes_msg.append("EggNOG Gene name added")
                        stats_genes_affected["eggnog_gene"] += 1
                    if interpro_gene_added: 
                        changes_msg.append("InterPro Gene name added")
                        stats_genes_affected["interpro_gene"] += 1
                    if eggnog_gos_added > 0: 
                        changes_msg.append(f"{eggnog_gos_added} EggNOG GO terms added")
                        stats_genes_affected["eggnog_go_terms"] += 1
                    if interpro_gos_added > 0: 
                        changes_msg.append(f"{interpro_gos_added} InterPro GO terms added")
                        stats_genes_affected["interpro_go_terms"] += 1
                    if xrefs_added > 0: 
                        changes_msg.append(f"{xrefs_added} db_xrefs added")
                        stats_genes_affected["db_xrefs"] += 1
                    if notes_added > 0: 
                        changes_msg.append(f"{notes_added} notes added")
                        stats_genes_affected["notes"] += 1
                    if pseudogene_added: 
                        changes_msg.append("Pseudogene added")
                        stats_genes_affected["pseudogene"] += 1
                    
                    # Emits a debug log specific to this one gene locus
                    logger.debug(f"Locus [{locus_tag}] modified: {', '.join(changes_msg)}")
                # -----------------------------------------------------------------

    # Clean up empty gene qualifiers that GenBank hates
    for sequence in gbff_read:
        for feature in sequence.features:
            if "gene" in feature.qualifiers and not feature.qualifiers["gene"]:
                del feature.qualifiers["gene"]

    # Write the updated sequences to the output GenBank file
    logger.debug(f"Writing annotated GBFF to: {gbff_out}")
    with open(gbff_out, "w") as out_handle:
        SeqIO.write(gbff_read, out_handle, "genbank")

    logger.info(
        f"Merged annotations: {cds_count} CDS features processed, "
        f"{ipr_annotated} annotated from InterPro, {egg_annotated} annotated from EggNOG"
    )

    logger.info(
        f"Genes modified per information type -> "
        f"EggNOG Product Name: {stats_genes_affected['eggnog_product']} | "
        f"InterPro Product Name: {stats_genes_affected['interpro_product']} | "
        f"EggNOG Gene Name: {stats_genes_affected['eggnog_gene']} | "
        f"InterPro Gene Name: {stats_genes_affected['interpro_gene']} | "
        f"EggNOG GO Terms: {stats_genes_affected['eggnog_go_terms']} | "
        f"InterPro GO Terms: {stats_genes_affected['interpro_go_terms']} | "
        f"db_xrefs (COG/KEGG/EC/IPR): {stats_genes_affected['db_xrefs']} | "
        f"Notes: {stats_genes_affected['notes']} | "
        f"Pseudogene: {stats_genes_affected['pseudogene']}"
    )
    return gbff_out