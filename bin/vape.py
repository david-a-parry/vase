#!/usr/bin/env python3

import sys
import argparse
sys.path.insert(0, '')
from vape.vape_runner import VapeRunner

def parse_args():
    parser = argparse.ArgumentParser(
        description='Filter, annotate and prioritize variants.',
        formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required_args = parser.add_argument_group('Required Arguments')
    optional_args = parser.add_argument_group('Optional Arguments')
    file_args = parser.add_argument_group('Annotation File Arguments')
    filter_args = parser.add_argument_group('Variant Filtering  Arguments')

    #required arguments
    required_args.add_argument(
'-i', '--input', required=True, metavar='VCF', help=
'''Input VCF filename

''')
    #misc optional arguments
    optional_args.add_argument(
'-o', '--output',  help=
#'-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help=
'''Filename for VCF output. If this ends in .gz or 
.bgz the output will be BGZIP compressed. 
Default = STDOUT

''')
    #args for filtering/retaining variants based on features
    filter_args.add_argument('--gq', type=int, 
                             help=
'''Minimum genotype quality score threshold. Genotype 
calls with a score lower than this threshold will 
be treated as no-calls

''')

    filter_args.add_argument(
'-v', '--variant_quality', type=float, metavar='QUAL', help=
'''Minimum variant quality score ('QUAL' field).
Variants with a QUAL score below this value will be 
filtered/ignored.

''') 
    filter_args.add_argument(
'-p', '--pass_filters', action='store_true', help=
'''Only keep variants that have passed filters 
(i.e. FILTER field must be "PASS")

''' )
    filter_args.add_argument(
'-c', '--csq', nargs='+', help=
'''One or more VEP consequence classes to  retain. 
Variants which do not result in  one of these VEP 
consequence classes will be filtered.  You may pass
the value  "default" in order to include the  
following classes: 
                  TFBS_ablation
                  TFBS_amplification
                  inframe_deletion
                  inframe_insertion
                  frameshift_variant
                  initiator_codon_variant
                  missense_variant
                  protein_altering_variant
                  regulatory_region_ablation
                  regulatory_region_amplification
                  splice_acceptor_variant
                  splice_donor_variant
                  stop_gained
                  stop_lost
                  transcript_ablation
                  transcript_amplification

''' ) 
    filter_args.add_argument(
'--canonical', action='store_true', help=
'''When used in conjunction with --csq argument, 
ignore consequences for non-canonical transcripts.

''')
    filter_args.add_argument(
'--keep_nmd_transcripts', action='store_true', default=False, help=
'''When used in conjunction with --csq  argument, 
also include variants labelled as an
'NMD_transcript_variant' if they also have a valid
consequence (as determined by the --csq argument).

''')


    #args for specifying files for annotations/filtering
    file_args.add_argument(
'-d', '--dbsnp', metavar='VCF', nargs='+', default=[], help=
'''dbSNP file for variant annotating/filteirng

''')
    file_args.add_argument(
'-g', '--gnomad', '--exac',  metavar='VCF', nargs='+', default=[], help=
'''gnomAD/ExAC file for variant annotating/filtering
using population allele frequencies

''')
    file_args.add_argument('-f', '--freq', type=float, help=
'''Allele frequency cutoff (between 0 and 1). Used 
for extenal allele frequency sources such as 
--dbsnp or --gnomad files. Alleles/variants with
an allele frequency equal to or greater than 
this value in these sources will be filtered 
from your input.

''')
    file_args.add_argument('--min_freq', type=float, help=
'''Minimum allele frequency cutoff (between 0 and 1).
Used for extenal allele frequency sources such as 
--dbsnp or --gnomad files. Alleles/variants with 
a frequency lower than this value will be filtered.

''')
    file_args.add_argument('--filter_novel', action='store_true', 
                           default=False, help=
'''Filter any allele/variant not present in 
any of the files supplied to --gnomad or --dbsnp
arguments.

''')
    return parser.parse_args()

if __name__ == '__main__':
    vape_args = parse_args()
    vape_runner = VapeRunner(vape_args)
    vape_runner.run()
