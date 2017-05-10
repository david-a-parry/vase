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

        super().__init__(vcf, prefix, freq, min_freq)


    def get_annot_fields(self):
        '''
            Creates dicts of INFO field names to dicts of 'Type', 
            'Number'and 'Description' as found in the VCF metadata for 
            INFO field names for frequency (freq_fields = ['AF']), 
            or allele annotations fields (annot_fields = ['AN', 'AC'])
        '''

        freq_fields = ["AF_POPMAX", "AF_AFR", "AF_AMR", "AF_EAS", "AF_FIN", 
                       "AF_NFE", "AF_SAS"]
        annot_fields = ["AC_POPMAX", "AC_AFR", "AC_AMR", "AC_EAS", "AC_FIN", 
                        "AC_NFE", "AC_SAS", "AN_POPMAX", "AN_AFR", "AN_AMR", 
                        "AN_EAS", "AN_FIN", "AN_NFE", "AN_SAS"]
        for f in freq_fields:
            if f in self.vcf.metadata['INFO']:
                self.freq_fields[f] = self.vcf.metadata['INFO'][f][-1]
        for f in annot_fields:
            if f in self.vcf.metadata['INFO']:
                self.annot_fields[f] = self.vcf.metadata['INFO'][f][-1]

        if not self.freq_fields and (self.freq is not None or 
                                   self.min_freq is not None):
            raise Exception("ERROR: no frequency fields identified in VCF " + 
                            "header for file '{}'.".format(self.vcf.filename) +
                            " Unable to use freq/min_freq arguments for " + 
                            "variant filtering.")


