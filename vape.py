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
    filter_args = parser.add_argument_group('Variant Filtering  Arguments', 
        'Arguments for filtering based on variant features')
    sample_args = parser.add_argument_group('Sample Based Filtering Arguments',
        'Arguments for filtering variants based on presence/absence in ' + 
        'samples and/or inheritance patterns')

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
    optional_args.add_argument(
'--quiet', action='store_true', help=
'''Do not output progress information.

''')
    optional_args.add_argument(
'--debug', action='store_true', help=
'''Output debugging information.

''')
    #args for filtering/retaining variants based on features
    filter_args.add_argument(
'-v', '--variant_quality', type=float, metavar='QUAL', help=
'''Minimum variant quality score ('QUAL' field).
Variants with a QUAL score below this value will be 
filtered/ignored.

''') 
    filter_args.add_argument(
'-p', '--pass_filters', action='store_true', default=False, help=
'''Only keep variants that have passed filters 
(i.e. FILTER field must be "PASS")

''' )
    filter_args.add_argument(
'-c', '--csq', nargs='*', help=
'''One or more VEP consequence classes to  retain. 
Variants which do not result in  one of these VEP 
consequence classes will be filtered. If no 
values are passed then the following default 
classes will be used:

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

You may also pass the value "default" in order to 
include these default classes in addition to other 
specified classes.

''' ) 
    filter_args.add_argument(
'--canonical', action='store_true', help=
'''When used in conjunction with --csq argument, 
ignore consequences for non-canonical transcripts.

''')
    filter_args.add_argument(
'--biotypes', nargs='+', default=[], metavar='BIOTYPE', help=
'''When used in conjunction with --csq  argument, 
ignore consequences in biotypes other than those 
specified here. By default only consequences in 
features with the following biotypes are 
considered:
    
            3prime_overlapping_ncrna
            antisense
            CTCF_binding_site
            enhancer
            IG_C_gene
            IG_D_gene
            IG_J_gene
            IG_V_gene
            lincRNA
            miRNA
            misc_RNA
            Mt_rRNA
            Mt_tRNA
            open_chromatin_region
            polymorphic_pseudogene
            processed_transcript
            promoter
            promoter_flanking_region
            protein_coding
            rRNA
            sense_intronic
            sense_overlapping
            snoRNA
            snRNA
            TF_binding_site
            translated_processed_pseudogene
            TR_C_gene
            TR_D_gene
            TR_J_gene
            TR_V_gene

Use this argument to specify one or more biotypes
to consider instead of those listed above. You may
also include the value 'default' in your list to
include the default values listed above in
addition to others provided to this argument.
Alternatively you may use the value 'all' to
disable filtering on biotypes.

''')
    filter_args.add_argument( 
'-m', '--missense_filters', default=[], nargs='+', help=
'''A list of in silico prediction programs to use
use for filtering missense variants (must be used
in conjunction with --csq argument). The programs 
provided her must have been annotated on the 
input VCF file either directly by VEP or via the 
dbNSFP VEP plugin. Recognised program names and 
default 'damaging' values are provided in the
"data/vep_insilico_pred.tsv" file. 

You may optionally specify score criteria for 
filtering as in the the following examples:
                    
    FATHMM_pred=D
    MutationTaster_pred=A
    MetaSVM_rankscore=0.8

Or you may just provide the program names 
and the default 'damaging' prediction values 
will be used, as listed in the file
"../data/vep_insilico_pred.tsv".

By default, a missense consequence is filtered 
unless each of the programs listed here have an 
appropriate or missing prediction/score. This 
behaviour can be changed using the 
--filter_unpredicted or --keep_if_any_damaging 
flags.

''')
    filter_args.add_argument(
'--filter_unpredicted', action='store_true', default=False, help=
'''For use in conjunction with --missense_filters. 
The default behaviour when using 
--missense_filters is to ignore a program if 
there is no prediction given (i.e. the score/pred 
is empty). That is, if there are no predictions 
for any of the programs annotating a missense 
consequence, it will not be filtered, while if 
predictions are missing for only some, filtering 
will proceed as normal with the other programs. If
this option is given, missense variants will be 
filtered if any program does not have a
prediction/score.

''')
    filter_args.add_argument(
'--keep_if_any_damaging', action='store_true', default=False, help=
'''For use in conjunction with --missense_filters. 
If this option is provided, a missense consequence
is only filtered if ALL of the programs provided 
to --missense_filters do not have an appropriate 
prediction/score - that is, the missense 
consequence will be retained if ANY of the given 
programs has an appropriate value for the 
prediction/score. This behaviour is overridden by
'--filter_unpredicted' when a prediction/score is 
missing for any program.

''')

    #args for specifying files for annotations/filtering
    file_args.add_argument(
'-d', '--dbsnp', '--clinvar', metavar='VCF', nargs='+', default=[], help=
'''dbSNP or ClinVar VCF file for variant 
annotating/filtering

''')
    file_args.add_argument(
'-g', '--gnomad', '--exac',  metavar='VCF', nargs='+', default=[], help=
'''gnomAD/ExAC file for variant annotating/filtering
using population allele frequencies

''')
    file_args.add_argument(
'-f', '--freq', type=float, help=
'''Allele frequency cutoff (between 0 and 1). Used 
for extenal allele frequency sources such as 
--dbsnp or --gnomad files. Alleles/variants with
an allele frequency equal to or greater than 
this value in these sources will be filtered 
from your input.

''')
    file_args.add_argument(
'--min_freq', type=float, help=
'''Minimum allele frequency cutoff (between 0 and 1).
Used for extenal allele frequency sources such as 
--dbsnp or --gnomad files. Alleles/variants with 
a frequency lower than this value will be filtered.

''')
 
    file_args.add_argument(
'-b', '--build', type=int, metavar='dbSNP_build', help=
'''dbSNP build version cutoff. For use with --dbsnp 
files. Alleles/variants present in this dbSNP 
build or earlier will be filtered from input.
from your input.

''')
    file_args.add_argument(
'--max_build', type=int, metavar='dbSNP_build', help=
'''Maximum dbSNP build version cutoff. For use with 
--dbsnp files. Alleles/variants present in dbSNP
builds later than this version will be filtered.

''')
    file_args.add_argument(
'--filter_novel', action='store_true', default=False, help=
'''Filter any allele/variant not present in 
any of the files supplied to --gnomad or --dbsnp
arguments.

''')
    file_args.add_argument(
'--clinvar_path', '-path', action='store_true', default=False, help=
'''Retain variants with ClinVar 'likely pathogenic' 
or 'pathogenic' flags regardless of frequency or 
other settings provided to other Annotation File 
Arguments. This requires one of the files 
provided to --dbsnp to have CLNSIG annotations 
from ClinVar.

''')

    #args for sample based filtering
    sample_args.add_argument(
'-cases', '--cases', default=[], nargs='+', metavar='SAMPLE_ID', help=
'''One or more sample IDs to treat as cases. Default 
behaviour is to retain variants/alleles present in 
all of these samples as long as they are not 
present in any sample specified using the 
'--controls' option. This behaviour can be 
adjusted using other options detailed below.

''')
    sample_args.add_argument(
'-controls', '--controls', default=[], nargs='+', metavar='SAMPLE_ID', help=
'''One or more sample IDs to treat as controls. 
Default behaviour is to filter variants/alleles 
present in any of these samples. This behaviour 
can be adjusted using other options detailed 
below.

''')
    sample_args.add_argument(
'-ped', '--ped', help=
'''A ped file containing information about samples in
your VCF for use for filtering on affectation 
status and inheritance patterns.

A PED file is a white-space (space or tab) 
delimited file with the first six mandatory 
columns:

     Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype

Affection status should be coded:

    -9 missing 
     0 missing
     1 unaffected
     2 affected

''')
    sample_args.add_argument(
'-gq', '--gq', type=int, default=20, help=
'''Minimum genotype quality score threshold. Sample 
genotype calls with a score lower than this 
threshold will be treated as no-calls. 
Default = 20.

''')

    sample_args.add_argument(
'-n_cases', '--n_cases', type=int, help=
'''Instead of requiring a variant to be present in
ALL samples specified by --cases, require at least
this many cases.

''')
    sample_args.add_argument(
'-n_controls', '--n_controls', type=int, help=
'''Instead of filtering an allele/variant if present
in ANY sample specified by --controls, require at 
least this many controls to carry a variant before
it is filtered.

''')
    sample_args.add_argument(
'--biallelic', action='store_true', help=
'''Identify variants matching a recessive inheritance
pattern in cases present in the PED file specified
by the --ped argument. Input must be VEP 
annotated. If the --csq argument is given, only
variants/alleles resulting in the given functional
consequences will be used to identify qualifying 
variants/alleles, otherwise the default set of 
VEP consequences (see --csq argument for details)
will be used.

''')
    sample_args.add_argument(
'--de_novo', action='store_true', help=
'''Idenfify apparent de novo variants in cases
present in the PED file specified by the --ped
argument. This requires that at least one 
parent-child trio exists in the given PED file.

''')
    #end of args
    return parser.parse_args()

if __name__ == '__main__':
    vape_args = parse_args()
    vape_runner = VapeRunner(vape_args)
    vape_runner.run()
