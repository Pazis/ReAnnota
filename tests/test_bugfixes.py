"""Regression tests for the October 2026 merge fixes (bugs 1, 2, 3, 5, 6, 7).

Runs the parsers and the merge once on the bundled example genome
(tests/fixtures) and checks invariants that the old code violated.
The GO checks need examples/output/go-basic.obo (download it with
examples/example_run.sh); they are skipped if the file is missing.

Run:  pytest tests/test_bugfixes.py -p no:cacheprovider --no-cov
"""

from pathlib import Path

import pytest
from Bio import SeqIO

from reannota.converters.merge_gbff import PSEUDO_CANDIDATES_KEEP_HYPOTHETICAL, ROOT_GO_TERMS
from reannota.converters import merge_csv_to_gbff
from reannota.parsers import (
    build_egg_dictionary_clean,
    egg_dict_to_tsv,
    ipr_dictotsv,
    ipr_termfinder,
    load_antismash_gbk,
    parse_pseudogff_to_dict,
)

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "fixtures"
OBO = ROOT / "examples" / "output" / "go-basic.obo"


def _cds(path):
    return {
        f.qualifiers["locus_tag"][0]: f.qualifiers
        for rec in SeqIO.parse(str(path), "genbank")
        for f in rec.features
        if f.type == "CDS" and "locus_tag" in f.qualifiers
    }


def _raw_eggnog():
    rows = {}
    with open(FIX / "eggnogmapper_test.annotations") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            rows[p[0]] = p
    return rows


@pytest.fixture(scope="module")
def merged(tmp_path_factory):
    out = tmp_path_factory.mktemp("merge")
    egg_tsv, ipr_tsv, gbff = out / "egg.tsv", out / "ipr.tsv", out / "enhanced.gbff"
    egg_dict_to_tsv(build_egg_dictionary_clean(str(FIX / "eggnogmapper_test.annotations")), str(egg_tsv))
    ipr_dictotsv(ipr_termfinder(str(FIX / "interproscan_test.gff3")), str(ipr_tsv))
    merge_csv_to_gbff(
        str(egg_tsv), str(ipr_tsv), str(FIX / "bakta_test.gbff"), str(gbff),
        pseudofile=str(FIX / "pseudofinder_test.gff"),
        go_obo_file=str(OBO) if OBO.exists() else None,
    )
    return {"egg_tsv": egg_tsv, "before": _cds(FIX / "bakta_test.gbff"), "after": _cds(gbff)}


# --- Bug 1: EC numbers written as "EC:EC:x.x.x.x" -------------------------
def test_no_double_ec_prefix(merged):
    ecs = [x for q in merged["after"].values() for x in q.get("db_xref", []) if x.startswith("EC:")]
    assert ecs, "expected EC cross-references"
    assert not any(x.startswith("EC:EC:") for x in ecs)


# --- Bug 3: only the first KO / EC of a protein was kept ------------------
def test_all_kos_and_ecs_written(merged):
    raw, before, after = _raw_eggnog(), merged["before"], merged["after"]
    checked = 0
    for tag, p in raw.items():
        if tag not in before or float(p[2]) >= 1:
            continue
        bakta = before[tag].get("db_xref", [])
        new = set(after[tag].get("db_xref", []))
        if p[11] != "-" and "," in p[11] and not any(x.startswith("KEGG:") for x in bakta):
            kos = {"KEGG:" + k.strip().replace("ko:", "") for k in p[11].split(",")}
            assert kos <= new, (tag, kos - new)
            checked += 1
        if p[10] != "-" and "," in p[10] and not any(x.startswith("EC:") for x in bakta):
            ecs = {"EC:" + e.strip() for e in p[10].split(",")}
            assert ecs <= new, (tag, ecs - new)
            checked += 1
    assert checked > 0, "fixture should contain multi-KO/EC proteins"


# --- Bug 2: propagated eggNOG GO lists were kept whole --------------------
@pytest.mark.skipif(not OBO.exists(), reason="go-basic.obo not downloaded")
def test_added_go_terms_are_most_specific(merged):
    from goatools.obo_parser import GODag

    dag = GODag(str(OBO))
    before, after = merged["before"], merged["after"]
    for tag, q in after.items():
        bakta = {x for x in before[tag].get("db_xref", []) if x.startswith("GO:")}
        added = {x for x in q.get("db_xref", []) if x.startswith("GO:")} - bakta
        assert not (added & ROOT_GO_TERMS), tag
        # no added term may be an ancestor of another added term or of a Bakta term
        everything = added | bakta
        ancestors = set()
        for g in everything:
            if g in dag:
                ancestors |= dag[g].get_all_parents()
        assert not (added & ancestors), (tag, added & ancestors)


# --- Bugs 5 + 6: pseudogene candidates escaped the flag -------------------
def test_pseudofinder_fragments_are_split():
    d = parse_pseudogff_to_dict(str(FIX / "pseudofinder_test.gff"))
    assert not any("," in k for k in d), "comma-joined locus tags must be split"
    frag = [v for v in d.values() if len(v["locus_tags"]) > 1]
    assert frag, "fixture should contain fragmented-gene calls"


def test_every_hypothetical_candidate_is_flagged(merged):
    before, after = merged["before"], merged["after"]
    pseudo = parse_pseudogff_to_dict(str(FIX / "pseudofinder_test.gff"))
    candidates = [
        t for t in pseudo
        if t in before
        and before[t].get("product", [""])[0] == "hypothetical protein"
        and not before[t].get("pseudogene")
    ]
    assert candidates
    for t in candidates:
        assert after[t].get("pseudogene"), t
        if PSEUDO_CANDIDATES_KEEP_HYPOTHETICAL:
            assert after[t]["product"][0] == "hypothetical protein", t


# --- Bug 7: antiSMASH protoclusters were counted instead of regions -------
def test_antismash_regions_not_protoclusters():
    n_regions = sum(
        1 for rec in SeqIO.parse(str(FIX / "antismash_test.gbk"), "genbank")
        for f in rec.features if f.type == "region"
    )
    bgcs = load_antismash_gbk(str(FIX / "antismash_test.gbk"))
    assert len(bgcs) == n_regions
    for b in bgcs:
        tags = [o.locus_tag for o in b.orfs]
        assert len(tags) == len(set(tags)), b.bgc_name
