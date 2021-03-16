class SvGtFilter(object):
    '''
        Given sample calls from a VariantRecord, this class provides a
        function 'filter' which returns True if the call meets the
        criteria (e.g. GQ, SR etc.) established on initialisation.
    '''

    __slots__ = ['gq', 'dp', 'max_dp', 'ad', 'max_ad', 'het_ab', 'hom_ab',
                 'gt_is_ok', 'ab_filter', 'ref_ab_filter', 'ad_over_threshold',
                 'fields', 'enough_support', 'del_dhffc', 'dup_dhbfc',
                 'duphold_filter']

    def __init__(self, vcf, gq=0, dp=0, max_dp=0, ad=0, max_ad=0, het_ab=0.,
                 hom_ab=0., ref_ab_filter=None, del_dhffc=None,
                 dup_dhbfc=None):
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
                        of supporting reads (SR + PR) for a genotype.
                        Default=0.

                max_dp: Semi-equivalent to maximum depth for short
                        variant analysis. Actually equates to the number
                        of supporting reads (SR + PR) for a genotype.
                        Default=0.

                ad:     Semi-equivalent to minimum allele depth for short
                        variant analysis. Actually equates to the number
                        of supporting reads (SR + PR) for an ALT allele.
                        Default=0.

                max_ad: Semi-equivalent to maximum depth for short
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

                del_dhffc:
                        Maximum fold-change for deletion calls relative
                        to flanking regions as annotated by duphold
                        (https://github.com/brentp/duphold). Deletion
                        calls will be filtered if the DHFFC annotation
                        from duphold is greater than this value.

                dup_dhbfc:
                        Minimum fold-change for duplicaton calls relative
                        to similar GC-content bins as annotated by
                        duphold (https://github.com/brentp/duphold).
                        Duplication calls will be filtered if the DHBFC
                        annotation from duphold is less than this value.

        '''
        self.gq = gq
        self.dp = dp
        self.max_dp = max_dp
        self.ad = ad
        self.max_ad = max_ad
        self.het_ab = het_ab
        self.hom_ab = hom_ab
        self.ref_ab_filter = ref_ab_filter
        self.del_dhffc = del_dhffc
        self.dup_dhbfc = dup_dhbfc
        self.fields = ['GT']
        self.ab_filter = None
        self.ad_over_threshold = None
        self.enough_support = None
        self.duphold_filter = None
        ab_fields = None
        if not any((gq, dp, max_dp, ad, max_ad, het_ab, hom_ab, del_dhffc,
                    dup_dhbfc)):
            # if no parameters are set then every genotype passes
            self.gt_is_ok = lambda gt, smp, al, svtype: True
        else:
            ab_fields = self._check_header_fields(vcf)
            if ab_fields == ('PR', 'SR'):
                # only option now, but may support other annotations in future
                self.ab_filter = self._ab_filter_prsr
                self.enough_support = self._enough_support_prsr
            if dup_dhbfc or del_dhffc:
                self.duphold_filter = self._duphold_filter
            self.gt_is_ok = self._gt_is_ok
        if ref_ab_filter:
            if ab_fields is None:
                ab_fields = self._check_header_fields(vcf)
            if ab_fields == ('PR', 'SR'):
                # only option now, but may support other annotations in future
                self.ad_over_threshold = self._alt_prsr_over_threshold

    def _duphold_filter(self, gts, sample, allele, svtype):
        is_hom_alt = False
        is_het_alt = False
        if len(set(gts[sample]['GT'])) == 1:
            if allele in gts[sample]['GT']:
                is_hom_alt = True
        elif allele in gts[sample]['GT']:
            is_het_alt = True
        if not is_hom_alt and not is_het_alt:
            return True  # do not filter
        if svtype == 'DUP' and self.dup_dhbfc:
            fc = gts[sample]['DHBFC']
            if fc is not None:
                return fc > self.dup_dhbfc
        if svtype == 'DEL' and self.del_dhffc:
            fc = gts[sample]['DHFFC']
            if fc is not None:
                return fc < self.del_dhffc
        return True  # do not filter

    def _alt_prsr_over_threshold(self, gts, sample, allele):
        support = self._get_pr_sr(gts, sample)
        al_dp = support[allele]
        dp = sum(support)
        if dp > 0 and al_dp is not None:
            ab = float(al_dp)/dp
            if ab > self.ref_ab_filter:
                # ALT/REF supporting read counts > threshold
                return True  # filter
        return False

    def _ab_filter_prsr(self, gts, sample, allele):
        support = self._get_pr_sr(gts, sample)
        al_dp = support[allele]
        dp = sum(support)
        is_hom_alt = False
        is_het_alt = False
        if len(set(gts[sample]['GT'])) == 1:
            if allele in gts[sample]['GT']:
                is_hom_alt = True
        elif allele in gts[sample]['GT']:
            is_het_alt = True
        if al_dp is not None and dp > 0 and (is_het_alt or is_hom_alt):
            ab = float(al_dp)/dp
            if is_het_alt and ab < self.het_ab:
                return False  # filter
            if is_hom_alt and ab < self.hom_ab:
                return False  # filter
        return True  # do not filter

    def _enough_support_prsr(self, gts, sample):
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
        pr = gts[sample].get('PR', (None,))
        sr = gts[sample].get('SR', (None,))
        if pr is None or pr == (None,):  # no PR values - just check SR
            pr = (0, 0)
        if sr is None or sr == (None,):  # no SR values - just check PR
            sr = (0, 0)
        return tuple(sum(t) for t in zip(sr, pr))

    def _gt_is_ok(self, gts, sample, allele, svtype):
        '''
            Returns True if genotype (from VariantRecord.samples) passes
            all parameters set on initialisation.
        '''
        if self.dp or self.max_dp:
            if not self.enough_support(gts, sample):
                return False
        if (self.ad or self.max_ad) and allele in gts[sample]['GT']:
            support = self._get_pr_sr(gts, sample)
            if self.ad and support[allele] < self.ad:
                return False
            if self.max_ad and support[allele] > self.max_ad:
                return False
        if self.gq:
            # if GQ is None presumably is a no call
            gq = gts[sample].get('GQ', None)
            if gq is None or gq < self.gq:
                return False
        if self.ab_filter is not None:
            if not self.ab_filter(gts, sample, allele):
                return False
        if self.duphold_filter:
            if not self.duphold_filter(gts, sample, allele, svtype):
                return False
        return True  # passes all filters

    def _check_header_fields(self, vcf):
        ''' Ensure the required annotations are present in VCF header. '''
        # DP and GQ are common fields that may not be defined in header
        if self.gq:
            self.fields.append('GQ')
        if self.dup_dhbfc:
            if 'DHBFC' not in vcf.header.header.formats:
                raise RuntimeError("Genotype filtering on SV duphold DHBFC " +
                                   "annotation is set but 'DHBFC' FORMAT " +
                                   "field is not defined in your VCF header.")
            self.fields.append("DHBFC")
        if self.del_dhffc:
            if 'DHFFC' not in vcf.header.header.formats:
                raise RuntimeError("Genotype filtering on SV duphold DHFFC " +
                                   "annotation is set but 'DHFFC' FORMAT " +
                                   "field is not defined in your VCF header.")
            self.fields.append("DHFFC")
        if (self.dp or self.max_dp or self.het_ab or self.hom_ab or
                self.ref_ab_filter):
            if ('PR' in vcf.header.header.formats and
                    'SR' in vcf.header.header.formats):
                self.fields.append('PR')
                self.fields.append('SR')
                return ('PR', 'SR')
            else:
                raise RuntimeError("Genotype filtering on SV allele balance " +
                                   "or supporting read depth is set but 'PR'" +
                                   " and/or 'SR' FORMAT fields are not " +
                                   "defined in your VCF header.")
        return None

