import sys
from .vcf_filter import *


class GnomadFilter(VcfFilter):
    '''
        An object that filters VCF records based on variant frequency
        annotations in a gnomAD or ExAC VCF file.
    '''

    def __init__(self, vcf, prefix, freq=None, min_freq=None, pops=None):
        '''
            Initialize object with a VCF file and optional filtering
            arguments.

            Args:
                vcf:      VCF containing variants to use to filter or
                          annotate records.

                prefix:   Prefix to prepend to added INFO field
                          annotations. Required.

                freq:     Filter alleles if allele frequency is greater
                          than this value. Optional.

                min_freq: Filter alleles if allele frequency is less
                          than this value. Optional.

                pops:     gnomAD population annotations to use. Default
                          are AFR, AMR, EAS, FIN, NFE and SAS.

        '''
        if pops is None:
            pops = ["AFR", "AMR", "EAS", "FIN", "NFE", "SAS"]
        freq_info = ["AF_" + p for p in pops]
        ac_info = ["AC_" + p for p in pops]
        an_info = ["AN_" + p for p in pops]
        super().__init__(vcf, prefix, freq, min_freq, freq_info, ac_info,
                         an_info)


