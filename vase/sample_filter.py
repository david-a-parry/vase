from .sv_gt_filter import SvGtFilter
import warnings


class SampleFilter(object):
    ''' A class for filtering VCF records on sample calls'''

    def __init__(self, vcf, cases=[], controls=[], n_cases=0, n_controls=0,
                 confirm_missing=False, gq=0, dp=0, max_dp=0, het_ab=0.,
                 hom_ab=0., min_control_gq=None, min_control_dp=None,
                 max_control_dp=None, control_het_ab=None, control_hom_ab=None,
                 con_ref_ab=None, sv_gq=0, sv_dp=0, sv_max_dp=0, sv_het_ab=0.,
                 sv_hom_ab=0., sv_min_control_gq=None, sv_min_control_dp=None,
                 sv_max_control_dp=None, sv_control_het_ab=None,
                 sv_control_hom_ab=None, sv_con_ref_ab=None, del_dhffc=None,
                 dup_dhbfc=None, control_del_dhffc=None,
                 control_dup_dhbfc=None):
        '''
            Initialize filtering options.

            Args:
                vcf:    VcfReader object for input. Sample IDs will be
                        checked against header information.

                cases:  IDs of samples that you want to check for the
                        presence of an ALT allele.

                controls:
                        IDs of samples that you want to confirm absence
                        of an ALT allele.

                n_cases:
                        Minimum number of cases required to carry an ALT
                        allele. Alleles will be filtered unless carried
                        by this number of cases or more. Default=0.

                n_controls:
                        Minimum number of controls required to carry an
                        ALT allele for it to be filtered. Alleles will
                        only be filtered if carried by this number of
                        controls or more. Default=0.

                gq:     Minimum genotype quality score (GQ). Genotype
                        calls with a GQ lower than this value will be
                        treated as no-calls. Default=0.

                dp:     Minimum genotype depth (DP). Genotype calls with
                        a DP lower than this value will be treated as
                        no-calls. Default=0.

                max_dp: Maximum genotype depth (DP). Genotype calls with
                        a DP higher than this value will be treated as
                        no-calls. Default=0 (i.e. not used).

                het_ab: Minimum genotype allele balance for heterozygous
                        calls. Genotype calls with an allele balance
                        lower than this value will be treated as
                        no-calls. The allele balance is calculated using
                        'AD' fields if present, otherwise 'AO' and 'RO'
                        fields (e.g. from freebayes). If none of these
                        fields are present in the VCF header and ab is
                        not 0.0 a RuntimeError will be thrown.
                        Default=0.0.

                hom_ab: As above but for homozygous genotype calls.

                min_control_gq:
                        Same as 'gq' but specific to control samples.
                        Defaults to the same as 'gq'.

                min_control_dp:
                        Same as 'dp' but specific to control samples.
                        Defaults to the same as 'dp'.

                max_control_dp:
                        Same as 'max_dp' but specific to control samples.
                        Defaults to the same as 'max_dp'.

                control_het_ab:
                        Same as 'het_ab' but specific to control samples.
                        Defaults to the same as 'het_ab'.

                control_hom_ab:
                        Same as 'hom_ab' but specific to control samples.
                        Defaults to the same as 'hom_ab'.

                con_ref_ab:
                        If a control sample has an ALT allele balance
                        equal to or greater than this value, consider
                        this control sample as carrying this allele
                        despite being called as 0/0.

                sv_gq:  Minimum genotype quality score (GQ) for
                        structural variants only. Defaults to the same
                        value as 'gq'.

                sv_dp:  Minimum number of supporting reads (SR + PR) for
                        structural variant calls. Genotype calls with
                        a fewer supporting reads than this value will be
                        treated as no-calls. Defaults to same as 'dp'.

                sv_max_dp:
                        Maximum number of supporting reads (SR + PR) for
                        structural variant calls. Genotype calls with
                        a more supporting reads than this value will be
                        treated as no-calls. Defaults to same as 'dp'.

                sv_het_ab:
                        Minimum allele balance for heterozygous
                        genotypefor structural variants. This is
                        calculated using SR + PR FORMAT fields, such as
                        provided by Manta. Defaults to the same as
                        het_ab.

                sv_hom_ab:
                        Minimum allele balance for homozygous
                        genotypes for structural variants. This is
                        calculated using SR + PR FORMAT fields, such as
                        provided by Manta. Defaults to the same as
                        hom_ab.

                sv_min_control_gq:
                        Same as 'sv_gq' but specific to control samples.
                        Defaults to the same as 'sv_gq'.

                sv_min_control_dp:
                        Same as 'sv_dp' but specific to control samples.
                        Defaults to the same as 'sv_dp'.

                sv_max_control_dp:
                        Same as 'sv_max dp' but specific to control
                        samples. Defaults to the same as 'sv_max_dp'.

                sv_control_het_ab:
                        Same as 'sv_het_ab' but specific to control
                        samples. Defaults to the same as 'sv_het_ab'.

                sv_control_hom_ab:
                        Same as 'sv_hom_ab' but specific to control
                        samples. Defaults to the same as 'sv_hom_ab'.

                sv_con_ref_ab:
                        If a control sample has an ALT allele balance
                        equal to or greater than this value (as
                        calculated using SR + PR FORMAT fields),
                        consider this control sample as carrying this
                        allele despite being called as 0/0.

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

                con_del_dhffc:
                        Same as del_dhffc but specific to control
                        samples.

                con_dup_dhbfc:
                        Same as dup_dhbfc but specific to control
                        samples.

                confirm_missing:
                        If True, only keep a variant if all controls are
                        not no-calls and are above the GQ threshold set
                        by the 'gq' option. Default=False.

        '''

        self.vcf = vcf
        self.confirm_missing = confirm_missing
        self._parse_sample_args(cases=cases, controls=controls,
                                n_cases=n_cases, n_controls=n_controls, gq=gq,
                                het_ab=het_ab, hom_ab=hom_ab, dp=dp,
                                max_dp=max_dp, con_gq=min_control_gq,
                                con_dp=min_control_dp,
                                con_max_dp=max_control_dp,
                                con_het_ab=control_het_ab,
                                con_hom_ab=control_hom_ab,
                                con_ref_ab=con_ref_ab, sv_gq=sv_gq,
                                sv_het_ab=sv_het_ab, sv_hom_ab=sv_hom_ab,
                                sv_dp=sv_dp, sv_max_dp=sv_max_dp,
                                sv_con_gq=sv_min_control_gq,
                                sv_con_dp=sv_min_control_dp,
                                sv_con_max_dp=sv_max_control_dp,
                                sv_con_het_ab=sv_control_het_ab,
                                sv_con_hom_ab=sv_control_hom_ab,
                                sv_con_ref_ab=sv_con_ref_ab,
                                del_dhffc=del_dhffc, dup_dhbfc=dup_dhbfc,
                                con_del_dhffc=control_del_dhffc,
                                con_dup_dhbfc=control_dup_dhbfc)

    def filter(self, v_record, allele):
        '''
            For a given VaseRecord and ALT allele, return True if that
            allele should be filtered based on presence/absence in cases
            and controls.
        '''
        record = v_record.record
        gts = record.samples
        case_matches = 0
        control_matches = 0
        svtype = None
        if v_record.IS_SV:
            gt_filter = self.sv_gt_filter
            control_filter = self.sv_con_gt_filter
            svtype = record.info['SVTYPE']
        else:
            gt_filter = self.gt_filter
            control_filter = self.con_gt_filter
        # check controls first
        for s in self.controls:
            gt_ok_args = [gts, s, allele]
            if svtype:
                gt_ok_args.append(svtype)
            if not control_filter.gt_is_ok(*gt_ok_args):
                if self.confirm_missing:
                    if self.n_controls:
                        control_matches += 1
                    else:
                        return True
                continue
            sgt = gts[s]['GT']
            if self.confirm_missing and sgt == (None,) * len(sgt):
                # no-call and we require confirmed gts for controls
                if self.n_controls:
                    control_matches += 1
                    continue
                else:
                    return True
            if allele in sgt:  # checks for presence, not whether het/hom
                if self.n_controls:
                    control_matches += 1
                else:
                    return True
            elif control_filter.ad_over_threshold is not None:
                # check hom ref for ALT allele counts
                if control_filter.ad_over_threshold(record.samples, s, allele):
                    if 'AD' not in record.format or gts[s]['AD'] != (None,):
                        if self.n_controls:
                            control_matches += 1
                        else:
                            return True
        if self.n_controls and control_matches >= self.n_controls:
            return True
        # check for presence in cases
        for s in self.cases:
            gt_ok_args = [gts, s, allele]
            if svtype:
                gt_ok_args.append(svtype)
            if not gt_filter.gt_is_ok(*gt_ok_args):
                sgt = None
            else:
                sgt = gts[s]['GT']
            if sgt is None:
                if self.n_cases:
                    continue
                else:
                    return True
            if allele in sgt:
                case_matches += 1
            elif not self.n_cases:
                return True
        if self.n_cases:
            if case_matches < self.n_cases:
                return True
        return False

    def _parse_sample_args(self, cases, controls, n_cases=0, n_controls=0,
                           gq=0, dp=0, max_dp=0, het_ab=0., hom_ab=0.,
                           con_gq=None, con_dp=None, con_max_dp=None,
                           con_het_ab=None, con_hom_ab=None, con_ref_ab=None,
                           sv_gq=0, sv_dp=0, sv_max_dp=None, sv_het_ab=0.,
                           sv_hom_ab=0., sv_con_gq=None, sv_con_dp=None,
                           sv_con_max_dp=None, sv_con_het_ab=None,
                           sv_con_hom_ab=None, sv_con_ref_ab=None,
                           del_dhffc=None, dup_dhbfc=None, con_del_dhffc=None,
                           con_dup_dhbfc=None):
        not_found = set()
        case_set = set()
        control_set = set()
        for c in cases:
            if c in self.vcf.header.samples:
                case_set.add(c)
            elif c.lower() == 'all':
                if len(cases) > 1:
                    raise RuntimeError("'all' can not be used in --cases " +
                                       "argument in conjunction with sample " +
                                       "names - please specify either 'all' " +
                                       "or a list of sample IDs, not both.")
                case_set.update(self.vcf.header.samples)
            else:
                not_found.add(c)
        if not_found:
            s = list(not_found)
            s.sort()
            raise RuntimeError("The following samples specified by --cases " +
                               "were not found in the input VCF: "
                               + str.join(", ", s))
        for c in controls:
            if c in self.vcf.header.samples:
                control_set.add(c)
            elif c.lower() == 'all':
                if len(controls) > 1:
                    raise RuntimeError("'all' can not be used in --controls " +
                                       "argument in conjunction with sample " +
                                       "names - please specify either 'all' " +
                                       "or a list of sample IDs, not both.")
                for samp in self.vcf.header.samples:
                    if samp not in case_set:
                        control_set.add(samp)
            else:
                not_found.add(c)
        if not_found:
            s = list(not_found)
            s.sort()
            raise RuntimeError("The following samples specified by " +
                               "--controls were not found in the input VCF: " +
                               str.join(", ", s))
        self.cases = list(case_set)
        self.controls = list(control_set)
        self.n_cases = None
        self.n_controls = None
        if n_cases and len(self.cases) < n_cases:
            raise RuntimeError("Number of cases specified by --n_cases is " +
                               "greater than the number of cases specified " +
                               "by --cases")
        if n_controls and len(self.controls) < n_controls:
            raise RuntimeError("Number of controls specified by --n_controls" +
                               "is greater than the number of controls " +
                               "specified by --controls")
        self.samples = self.cases + self.controls
        self.gt_filter = GtFilter(self.vcf, gq=gq, dp=dp, max_dp=max_dp,
                                  het_ab=het_ab, hom_ab=hom_ab)
        self.gt_fields = set(self.gt_filter.fields)
        if con_gq is None:
            con_gq = gq
        if con_dp is None:
            con_dp = dp
        if con_max_dp is None:
            con_max_dp = max_dp
        if con_het_ab is None:
            con_het_ab = het_ab
        if con_hom_ab is None:
            con_hom_ab = hom_ab
        self.con_gt_filter = GtFilter(self.vcf, gq=con_gq, dp=con_dp,
                                      max_dp=con_max_dp, het_ab=con_het_ab,
                                      hom_ab=hom_ab, ref_ab_filter=con_ref_ab)
        self.gt_fields.update(self.con_gt_filter.fields)
        if sv_gq is None:
            sv_gq = gq
        if sv_dp is None:
            sv_dp = dp
        if sv_max_dp is None:
            sv_max_dp = max_dp
        if sv_het_ab is None:
            sv_het_ab = het_ab
        if sv_hom_ab is None:
            sv_hom_ab = hom_ab
        if sv_con_gq is None:
            sv_con_gq = sv_gq
        if sv_con_dp is None:
            sv_con_dp = sv_dp
        if sv_con_max_dp is None:
            sv_con_max_dp = sv_max_dp
        if sv_con_het_ab is None:
            sv_con_het_ab = sv_het_ab
        if sv_con_hom_ab is None:
            sv_con_hom_ab = sv_hom_ab
        if con_del_dhffc is None:
            con_del_dhffc = del_dhffc
        if con_dup_dhbfc is None:
            con_dup_dhbfc = dup_dhbfc
        self.sv_gt_filter = SvGtFilter(self.vcf, gq=sv_gq, dp=sv_dp,
                                       max_dp=sv_max_dp, het_ab=sv_het_ab,
                                       hom_ab=sv_hom_ab, del_dhffc=del_dhffc,
                                       dup_dhbfc=dup_dhbfc)
        self.sv_gt_fields = set(self.sv_gt_filter.fields)
        self.sv_con_gt_filter = SvGtFilter(self.vcf, gq=sv_con_gq,
                                           dp=sv_con_dp, max_dp=sv_con_max_dp,
                                           het_ab=sv_con_het_ab,
                                           hom_ab=sv_hom_ab,
                                           del_dhffc=con_del_dhffc,
                                           dup_dhbfc=con_dup_dhbfc,
                                           ref_ab_filter=sv_con_ref_ab)
        self.sv_gt_fields.update(self.sv_con_gt_filter.fields)
        if n_cases:
            self.n_cases = n_cases
        if n_controls:
            self.n_controls = n_controls


class GtFilter(object):
    '''
        Given a dict of GT information from the 'parsed_gts' function of
        VcfRecord from vcf_record.py, this provides a function 'filter'
        which returns True if the call meets the criteria (e.g. GQ, AD
        etc.) established on initialisation.
    '''

    __slots__ = ['gq', 'dp', 'max_dp', 'het_ab', 'hom_ab', 'gt_is_ok',
                 'ab_filter', 'ref_ab_filter', 'ad_over_threshold', 'fields',
                 '_get_gq', '_format_ref_dp', '_format_alt_dp']

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

                dp:     Minimum depth (DP) for a genotype. Default=0.

                max_dp: Maximum depth (DP) for a genotype. Default=0 (i.e. not
                        used).

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
                        Maximum ALT allele balance of homozygous
                        reference calls. If a 0/0 genotype call has
                        this fraction of ALT alleles it will be
                        filtered (i.e. self.gt_is_ok will return
                        False). Default=None (not used).

        '''
        self.gq = gq
        self.dp = dp
        self.max_dp = max_dp
        self.het_ab = het_ab
        self.hom_ab = hom_ab
        self.ref_ab_filter = ref_ab_filter
        self.fields = ['GT']
        self._get_gq = self._get_gq_standard
        self.ab_filter = None
        self.ad_over_threshold = None
        self._format_ref_dp = None
        self._format_alt_dp = None
        ab_field = None
        if not gq and not dp and not max_dp and not het_ab and not hom_ab:
            # if no parameters are set then every genotype passes
            self.gt_is_ok = lambda gt, smp, al: True
        else:
            ab_field = self._check_header_fields(vcf)
            if het_ab or hom_ab:
                if ab_field == 'AD':
                    self.ab_filter = self._ab_filter_ad
                elif ab_field == 'RO' or ab_field == 'NR':
                    self.ab_filter = self._ab_filter_ro
            self.gt_is_ok = self._gt_is_ok
        if ref_ab_filter:
            if ab_field is None:
                ab_field = self._check_header_fields(vcf)
            if ab_field == 'AD':
                self.ad_over_threshold = self._alt_ad_over_threshold
            elif ab_field == 'RO' or ab_field == 'NR':
                self.ad_over_threshold = self._alt_ao_over_threshold

    def _get_gq_standard(self, gts, sample):
        return gts[sample].get('GQ', None)

    def _get_gq_from_tuple(self, gts, sample):
        gq = gts[sample].get('GQ', (None, ))
        return gq[0]

    def _alt_ad_over_threshold(self, gts, sample, allele):
        ad = gts[sample].get('AD', (None,))
        if ad == (None,):  # no AD values - assume OK?
            return True
        al_dp = ad[allele]
        dp = sum(filter(None, ad))
        if dp > 0 and al_dp is not None:
            ab = float(al_dp)/dp
            if ab > self.ref_ab_filter:  # ALT/REF read counts > threshold
                return True  # filter
        return False

    def _alt_ao_over_threshold(self, gts, sample, allele):
        aos = gts[sample].get(self._format_alt_dp, (None,))
        ro = gts[sample].get(self._format_ref_dp, None)
        if isinstance(ro, tuple):  # platypus
            ro = ro[0]
        if aos != (None,) and ro is not None:
            dp = sum(filter(None, aos)) + ro
            if dp > 0:
                ao = aos[allele-1]
                ab = float(ao)/ro
                if ab > self.ref_ab_filter:
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
        if self.het_ab or self.hom_ab or self.ref_ab_filter:
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
