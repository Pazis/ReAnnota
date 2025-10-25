

def parse_pseudogff_to_dict(gff_file):
    """
    Επιστρέφει dictionary με key = ID και value = dict με βασικά πεδία
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
                gff_dict[attr_dict['old_locus_tag']] = {
                    'seqid': seqid,
                    'source': source,
                    'type': type_,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'attributes': attributes
                }

    return gff_dict


