import warnings


class GtFilter(object):
    '''
        Given a dict of GT information from the 'parsed_gts' function of
        VcfRecord from vcf_record.py, this provides a function 'filter'
        which returns True if the call meets the criteria (e.g. GQ, AD
        etc.) established on initialisation.
    '''

    __slots__ = ['gq', 'dp', 'max_dp', 'ad', 'max_ad', 'het_ab', 'hom_ab',
                 'gt_is_ok', 'get_ad', 'ab_filter', 'ref_ab_filter',
                 'ref_ad_filter', 'ab_over_threshold', 'fields', '_get_gq',
                 '_format_ref_dp', '_format_alt_dp']

    def __init__(self, vcf, gq=0, dp=0, max_dp=0, ad=0, max_ad=0, het_ab=0.,
                 hom_ab=0., ref_ab_filter=None, ref_ad_filter=None):
        '''
            Args:
                vcf:    Input VCF, which will be checked to ensure any
                        specified filters can be used by ensuring
                        appropriate FORMAT fields are defined in the
                        VCF header.

                gq:     Minimum genotype quality (GQ) for a genotype.
                        Default=0 (i.e. not checked).

                dp:     Minimum depth (DP) for a genotype. Default=0.

                max_dp: Maximum depth (DP) for a genotype. Default=0 (i.e. not
                        used).

                ad:     Minimum allele depth for a genotype. Default=0.

                max_ad: Maximum allele depth for a genotype. Default=0 (i.e.
                        not used).

                het_ab: Minimum allele balance for a heterozygous
                        genotype. This is calculated using AD FORMAT
                        field if present. If AD is not present, AO and
                        RO will be used instead if present, otherwise
                        throws a RuntimeError. Default=0.

                hom_ab: Minimum allele balance for a homozygous
                        genotype. This is calculated using AD FORMAT
                        field if present. If AD is not present, AO and
                        RO will be used instead if present, otherwise
                        throws a RuntimeError. Default=0.

                ref_ab_filter:
                        ALT allele balance threshold for homozygous
                        reference calls. If a 0/0 genotype call has
                        more than this fraction of ALT alleles it will
                        be filtered (i.e. self.gt_is_ok will return
                        False). Default=None (not used).

                ref_ad_filter:
                        ALT allele count threshold for homozygous
                        reference calls. If a 0/0 genotype call has
                        this many of ALT alleles it will be
                        filtered (i.e. self.gt_is_ok will return
                        False). Default=None (not used).

        '''
        self.gq = gq
        self.dp = dp
        self.max_dp = max_dp
        self.ad = ad
        self.max_ad = max_ad
        self.het_ab = het_ab
        self.hom_ab = hom_ab
        self.ref_ab_filter = ref_ab_filter
        self.ref_ad_filter = ref_ad_filter
        self.fields = ['GT']
        self._get_gq = self._get_gq_standard
        self.get_ad = None
        self.ab_filter = None
        self.ab_over_threshold = None
        self._format_ref_dp = None
        self._format_alt_dp = None
        ab_field = None
        if not any((gq, dp, max_dp, ad, max_ad, het_ab, hom_ab)):
            # if no parameters are set then every genotype passes
            self.gt_is_ok = lambda gt, smp, al: True
        else:
            ab_field = self._check_header_fields(vcf)
            if self.ad or self.max_ad:
                if ab_field == 'AD':
                    self.get_ad = self._get_ad
                else:
                    self.get_ad = self._get_nonstandard_ad
            if het_ab or hom_ab:
                if ab_field == 'AD':
                    self.ab_filter = self._ab_filter_ad
                elif ab_field == 'RO' or ab_field == 'NR':
                    self.ab_filter = self._ab_filter_ro
            self.gt_is_ok = self._gt_is_ok
        if ref_ab_filter or ref_ad_filter:
            if ab_field is None:
                ab_field = self._check_header_fields(vcf)
            if ab_field == 'AD':
                self.ab_over_threshold = self._alt_ab_ad_over_threshold
                self.get_ad = self._get_ad
            elif ab_field == 'RO' or ab_field == 'NR':
                self.ab_over_threshold = self._alt_ab_ao_over_threshold
                self.get_ad = self._get_nonstandard_ad

    def _get_gq_standard(self, gts, sample):
        return gts[sample].get('GQ', None)

    def _get_gq_from_tuple(self, gts, sample):
        gq = gts[sample].get('GQ', (None, ))
        return gq[0]

    def _get_ad(self, gts, sample, allele):
        ad = gts[sample].get('AD', (None,))
        if ad != (None,):
            return ad[allele]
        return None

    def _get_nonstandard_ad(self, gts, sample, allele):
        ad = gts[sample].get(self._format_alt_dp)
        if ad is not None:
            return ad[allele - 1]
        return None

    def ad_over_threshold(self, gts, sample, allele):
        ad = self.get_ad(gts, sample, allele)
        if ad is None:  # no AD values - filter?
            return True
        return ad >= self.ref_ad_filter

    def _alt_ab_ad_over_threshold(self, gts, sample, allele):
        ad = gts[sample].get('AD', (None,))
        if ad == (None,):  # no AD values - filter?
            return True
        al_dp = ad[allele]
        dp = sum(filter(None, ad))
        if dp > 0 and al_dp is not None:
            ab = float(al_dp)/dp
            if ab >= self.ref_ab_filter:  # ALT/REF read counts > threshold
                return True  # filter
        return False

    def _alt_ab_ao_over_threshold(self, gts, sample, allele):
        aos = gts[sample].get(self._format_alt_dp, (None,))
        ro = gts[sample].get(self._format_ref_dp, None)
        if isinstance(ro, tuple):  # platypus
            ro = ro[0]
        if aos != (None,) and ro is not None:
            dp = sum(filter(None, aos)) + ro
            if dp > 0:
                ao = aos[allele-1]
                ab = float(ao)/ro
                if ab >= self.ref_ab_filter:
                    return True
        return False

    def _ab_filter_ad(self, gts, sample, allele):
        ad = gts[sample].get('AD', (None,))
        if ad == (None,):  # no AD values - assume OK?
            return True
        al_dp = ad[allele]
        dp = sum(filter(None, ad))
        is_hom_alt = False
        is_het_alt = False
        if len(set(gts[sample]['GT'])) == 1:
            if allele in gts[sample]['GT']:
                is_hom_alt = True
        elif allele in gts[sample]['GT']:
            is_het_alt = True
        if al_dp is not None and dp > 0 and (is_het_alt or is_hom_alt):
            ab = float(al_dp)/dp
            if is_het_alt and self.het_ab and ab < self.het_ab:
                return False  # filter
            if is_hom_alt and self.hom_ab and ab < self.hom_ab:
                return False  # filter
        return True  # do not filter

    def _ab_filter_ro(self, gts, sample, allele):
        aos = gts[sample].get(self._format_alt_dp, (None,))
        ro = gts[sample].get(self._format_ref_dp, None)
        if isinstance(ro, tuple):  # platypus
            ro = ro[0]
        if ro is None or aos == (None,):  # no AD values - assume OK?
            return True
        is_hom_alt = False
        is_het_alt = False
        if len(set(gts[sample]['GT'])) == 1:
            if allele in gts[sample]['GT']:
                is_hom_alt = True
        elif allele in gts[sample]['GT']:
            is_het_alt = True
        if aos is not None and ro is not None and (is_hom_alt or is_het_alt):
            dp = sum(filter(None, aos)) + ro
            if allele > 0:
                ao = aos[allele-1]
            else:
                ao = ro
            if dp > 0:
                ab = float(ao)/dp
                if is_het_alt and ab < self.het_ab:
                    return False  # filter
                if is_hom_alt and ab < self.hom_ab:
                    return False  # filter
        return True  # do not filter

    def _gt_is_ok(self, gts, sample, allele):
        '''
            Returns True if genotype (from VariantRecord.samples) passes
            all parameters set on initialisation.
        '''
        if self.dp or self.max_dp:
            dp = gts[sample].get('DP', None)
            if self.dp:
                if dp is not None and dp < self.dp:
                    return False
            if self.max_dp:
                if dp is not None and dp > self.max_dp:
                    return False
        if (self.ad or self.max_ad) and allele in gts[sample]['GT']:
            ad = self.get_ad(gts, sample, allele)
            if self.ad:
                if ad is not None and ad < self.ad:
                    return False
            if self.max_ad:
                if ad is not None and ad > self.max_ad:
                    return False
        if self.gq:  # if GQ is None do not filter(?)
            gq = self._get_gq(gts, sample)
            if gq is not None and gq < self.gq:
                return False
        if self.ab_filter is not None:
            if not self.ab_filter(gts, sample, allele):
                return False
        return True  # passes all filters

    def _check_header_fields(self, vcf):
        ''' Ensure the required annotations are present in VCF header. '''
        # DP and GQ are common fields that may not be defined in header
        if self.dp or self.max_dp:
            self.fields.append('DP')
        if self.gq:
            self.fields.append('GQ')
            if 'GQ' in vcf.header.formats:
                if vcf.header.formats['GQ'].number == '.':  # Platypus formats
                    self._get_gq = self._get_gq_from_tuple
        if any((self.ad, self.max_ad, self.het_ab, self.hom_ab,
                self.ref_ab_filter, self.ref_ad_filter)):
            if 'AD' in vcf.header.formats:
                self.fields.append('AD')
                return 'AD'
            elif ('AO' in vcf.header.formats and 'RO' in vcf.header.formats):
                self.fields.append('AO')
                self.fields.append('RO')
                self._format_ref_dp = 'RO'
                self._format_alt_dp = 'AO'
                return 'RO'
            elif ('NV' in vcf.header.formats and 'NR' in vcf.header.formats):
                self.fields.append('NV')
                self.fields.append('NR')
                self._format_ref_dp = 'NR'
                self._format_alt_dp = 'NV'
                return 'NR'
            else:
                warnings.warn("Genotype filtering on allele balance is " +
                              "set but neither 'AD' nor standard supported " +
                              "FORMAT fields from Freebayes or Platypus are " +
                              "defined in your VCF header.")
        return None
