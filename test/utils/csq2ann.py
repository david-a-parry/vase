#!/usr/bin/env python3
import sys
import pysam
from vase.vcf_reader import VcfReader

csq_fields = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene',
              'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON']
ann_fields = ['Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name',
              'Gene_ID', 'Feature_Type', 'Feature_ID', 'Transcript_BioType',
              'Rank']
csq2ann = dict(zip(csq_fields, ann_fields))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("Usage: {} input.vcf output.vcf".format(sys.argv[0]))
    vcf = VcfReader(sys.argv[1])
    meta_data = []
    vcf.header.info.add('ANN',
                        '.',
                        'String',
                        "Functional annotations: '" + ' | '.join(ann_fields)
                        + "'")
    out = pysam.VariantFile(sys.argv[2], 'w', header=vcf.header.header)
    out.header.info.remove_header('CSQ')
    for record in vcf:
        annots = []
        for csq in record.CSQ:
            ann = dict((csq2ann[k], v) for k, v in csq.items() if k in csq2ann)
            if not ann['Rank']:
                ann['Rank'] = csq['INTRON']
            ann['Allele'] = record.alleles[record._vep_to_alt(csq)]
            annots.append("|".join(ann[k] for k in ann_fields))
        _ = record.info.pop('CSQ')
        record.info['ANN'] = annots
        out.write(record.record)
    out.close()
