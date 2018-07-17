class SvGtFilter(object):
    '''
        Given a dict of GT information from the 'parsed_gts' function of
        VcfRecord from parse_vcf.py, this provides a function 'filter'
        which returns True if the call meets the criteria (e.g. GQ, SR,
        etc.) established on initialisation.
    '''

    __slots__ = ['gq', 'dp', 'max_dp', 'het_ab', 'hom_ab', 'gt_is_ok',
                 'ab_filter', 'ref_ab_filter', 'ad_over_threshold', 'fields',
                 'enough_support']

    def __init__(self, vcf, gq=0, dp=0, max_dp=0, het_ab=0., hom_ab=0.,
                 ref_ab_filter=None):
        '''
            Args:
                vcf:    Input VCF, which will be checked to ensure any
                        specified filters can be used by ensuring
                        appropriate FORMAT fields are defined in the
                        VCF header.

                gq:     Minimum genotype quality (GQ) for a genotype.
                        Default=0 (i.e. not checked).

                dp:     Semi-equivalent to minimum depth for short
                        variant analysis. Actually equates to the number
                        of supporting reads (SR + PR) for an ALT allele.
                        Default=0.

                max_dp: Semi-equivalent to maximum depth for short
                        variant analysis. Actually equates to the number
                        of supporting reads (SR + PR) for an ALT allele.
                        Default=0.

                het_ab: Minimum allele balance for a heterozygous
                        genotype. This is calculated using SR + PR FORMAT
                        fields, such as provided by Manta. Default=0.

                hom_ab: Minimum allele balance for a homozygous
                        genotype. This is calculated using SR + PR FORMAT
                        fields, such as provided by Manta. Default=0.

                ref_ab_filter:
                        Maximum ALT allele balance of homozygous
                        reference calls. If a 0/0 genotype call has
                        this fraction of supporting reads (as calculated
                        using SR + PR FORMAT fields) it will be filtered
                        (i.e. self.gt_is_ok will return False).
                        Default=None (not used).

        '''
        self.gq = gq
        self.dp = dp
        self.max_dp = max_dp
        self.het_ab = het_ab
        self.hom_ab = hom_ab
        self.ref_ab_filter = ref_ab_filter
        self.fields = ['GT']
        self.ab_filter = None
        self.ad_over_threshold = None
        self.enough_support = None
        ab_fields = None
        if not gq and not dp and not het_ab and not hom_ab:
            #if no parameters are set then every genotype passes
            self.gt_is_ok = lambda gt, smp, al: True
        else:
            ab_fields = self._check_header_fields(vcf)
            if ab_fields == ('PR', 'SR'):
                #only option now, but may support other annotations in future
                self.ab_filter = self._ab_filter_prsr
                self.enough_support = self._enough_support_prsr
            self.gt_is_ok = self._gt_is_ok
        if ref_ab_filter:
            if ab_fields is None:
                ab_fields = self._check_header_fields(vcf)
            if ab_fields == ('PR', 'SR'):
                #only option now, but may support other annotations in future
                self.ad_over_threshold = self._alt_prsr_over_threshold

    def _alt_prsr_over_threshold(self, gts, sample, allele):
        support = self._get_pr_sr(gts, sample)
        al_dp = support[allele]
        dp = sum(support)
        if dp > 0 and al_dp is not None:
            ab = float(al_dp)/dp
            if ab > self.ref_ab_filter:
                #ALT/REF supporting read counts > threshold
                return True #filter
        return False

    def _ab_filter_prsr(self, gts, sample, allele):
        support = self._get_pr_sr(gts, sample)
        al_dp = support[allele]
        dp = sum(support)
        is_hom_alt = False
        is_het_alt = False
        if len(set(gts['GT'][sample])) == 1:
            if allele in gts['GT'][sample]:
                is_hom_alt = True
        elif allele in gts['GT'][sample]:
            is_het_alt = True
        if al_dp is not None and dp > 0 and (is_het_alt or is_hom_alt):
            ab = float(al_dp)/dp
            if is_het_alt and ab < self.het_ab:
                return False #filter
            if is_hom_alt and ab < self.hom_ab:
                return False #filter
        return True #do not filter

    def _enough_support_prsr(self, gts, sample, allele):
        '''
            Returns False if 'SR' + 'PR' fields for given sample
            genotype is < self.dp.
        '''
        support = self._get_pr_sr(gts, sample)
        if self.dp and sum(support) < self.dp:
            return False
        if self.max_dp and sum(support) > self.max_dp:
            return False
        return True

    def _get_pr_sr(self, gts, sample):
        ''' Returns a tuple of SR + PR counts for a sample.'''
        pr = gts['PR'][sample]
        sr = gts['SR'][sample]
        if pr == (None,): #no PR values - just check SR
            pr = (0, 0)
        if sr == (None,): #no SR values - just check PR
            sr = (0, 0)
        return tuple(sum(t) for t in zip(sr, pr))

    def _gt_is_ok(self, gts, sample, allele):
        '''
            Returns True if genotype (from parse_vcf.py parsed_gts
            function) passes all parameters set on initialisation.
        '''
        if self.dp or self.max_dp:
            if not self.enough_support(gts, sample, allele):
                return False
        if self.gq:
            #if GQ is None presumably is a no call
            if gts['GQ'][sample] is None or gts['GQ'][sample] < self.gq:
                return False
        if self.ab_filter is not None:
            if not self.ab_filter(gts, sample, allele):
                return False
        return True #passes all filters

    def _check_header_fields(self, vcf):
        ''' Ensure the required annotations are present in VCF header. '''
        #DP and GQ are common fields that may not be defined in header
        if self.gq:
            self.fields.append('GQ')
        if self.dp or self.het_ab or self.hom_ab or self.ref_ab_filter:
            if ('PR' in vcf.header.metadata['FORMAT'] and
                'SR' in vcf.header.metadata['FORMAT']):
                self.fields.append('PR')
                self.fields.append('SR')
                return ('PR', 'SR')
            else:
                raise RuntimeError("Genotype filtering on SV allele balance " +
                                   "or supporting read depth is set but 'PR'" +
                                   " and/or 'SR' FORMAT fields are not " +
                                   "defined in your VCF header.")
        return None

