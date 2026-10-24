def parse_pseudogff_to_dict(gff_file):
    """
    Parse a Pseudofinder GFF and return a dictionary keyed by the original
    (Bakta) locus tag, with the basic fields of the pseudogene call as value.

    Pseudofinder reports fragmented genes ("Predicted fragmentation of a
    single gene") with several comma-separated tags in old_locus_tag, e.g.
    old_locus_tag=JDBBEK_05692,JDBBEK_05693. Each tag is registered
    separately so that every fragment can be matched to its CDS.
    """
    gff_dict = {}

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid, source, type_, start, end, score, strand, phase, attributes = parts

            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            if 'old_locus_tag' in attr_dict:
                tags = [t.strip() for t in attr_dict['old_locus_tag'].split(',') if t.strip()]
                record = {
                    'seqid': seqid,
                    'source': source,
                    'type': type_,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'attributes': attributes,
                    'locus_tags': tags,
                }
                for tag in tags:
                    gff_dict[tag] = record

    return gff_dict
