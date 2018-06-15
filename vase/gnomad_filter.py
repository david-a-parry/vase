import sys
from .vcf_filter import *


class GnomadFilter(VcfFilter):
    '''
        An object that filters VCF records based on variant frequency
        annotations in a gnomAD or ExAC VCF file.
    '''

    def __init__(self, vcf, prefix, freq=None, min_freq=None):
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

        '''

        freq_info = ("AF_POPMAX", "AF_AFR", "AF_AMR", "AF_EAS", "AF_FIN",
                       "AF_NFE", "AF_SAS")
        ac_info = ("AC_POPMAX", "AC_AFR", "AC_AMR", "AC_EAS", "AC_FIN",
                      "AC_NFE", "AC_SAS",)
        an_info = ("AN_POPMAX", "AN_AFR", "AN_AMR", "AN_EAS", "AN_FIN",
                    "AN_NFE", "AN_SAS",)
        super().__init__(vcf, prefix, freq, min_freq, freq_info, ac_info,
                         an_info)


