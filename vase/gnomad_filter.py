from .vcf_filter import VcfFilter


class GnomadFilter(VcfFilter):
    '''
        An object that filters VCF records based on variant frequency
        annotations in a gnomAD or ExAC VCF file.
    '''

    def __init__(self, vcf, prefix, freq=None, min_freq=None, pops=None,
                 max_homozygotes=None):
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

                max_homozygotes:
                          Filter alleles if the total number of homozygotes
                          or hemizygotes is equal to or greater than this
                          value.

        '''
        if pops is None:
            pops = ["AFR", "AMR", "EAS", "FIN", "NFE", "SAS"]
        #v2.1 put pops in lowercase, previously were in uppercase - allow both
        g_pops = [x.lower() for x in pops] + [x.upper() for x in pops]
        freq_info = ["AF_" + p for p in g_pops]
        ac_info = ["AC_" + p for p in g_pops]
        an_info = ["AN_" + p for p in g_pops]
        hom_info = ["Hom_" + p for p in g_pops] + ["nhomalt_" + p for p in
                                                   g_pops]
        hemi_info = ["Hemi_" + p for p in g_pops] #v2.1 uses nhomalt for hemi
        super().__init__(vcf=vcf, prefix=prefix, freq=freq, min_freq=min_freq,
                         freq_fields=freq_info, ac_fields=ac_info,
                         an_fields=an_info, annotations=hom_info+hemi_info,
                         allow_missing_annotations=True)
        self.hom_annots = [self.prefix + "_" + f for f in self.extra if f in
                           self.annot_fields]
        self.max_homozygotes = max_homozygotes
        if self.max_homozygotes is not None:
            self.annotate_and_filter_record = self.filter_including_homozygotes
        else:
            self.annotate_and_filter_record = super().annotate_and_filter_record


    def filter_including_homozygotes(self, record):
        filt, keep, matched = super().annotate_and_filter_record(record)
        if all(filt):
            #no need to check homozygotes - filtering anyway
            return filt, keep, matched
        hom_info = record.parsed_info_fields(self.hom_annots)
        if hom_info:
            for i in range(len(record.DECOMPOSED_ALLELES)):
                if filt[i]:
                    continue
                hom_counts = sum(hom_info[x][i] for x in self.hom_annots if
                                 x in hom_info and hom_info[x][i] is not None)
                if hom_counts > self.max_homozygotes:
                    filt[i] = True
        return filt, keep, matched
