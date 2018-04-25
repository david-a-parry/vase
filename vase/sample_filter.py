

class SampleFilter(object):
    ''' A class for filtering VCF records on sample calls'''

    def __init__(self, vcf, cases=[], controls=[], n_cases=0, n_controls=0, 
                 confirm_missing=False, gq=0, dp=0, het_ab=0., hom_ab=0., 
                 min_control_gq=None, min_control_dp=None, control_het_ab=None, 
                 control_hom_ab=None, con_ref_ab=None):
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
                        a GQ lower than this value will be treated as 
                        no-calls. Default=0.

                ab:     Minimum genotype allele balance. Genotype calls
                        an allele balance lower than this value will be 
                        treated as no-calls. The allele balance is 
                        calculated using 'AD' fields if present, 
                        otherwise 'AO' and 'RO' fields (e.g. from 
                        freebayes). If none of these fields are present
                        in the VCF header and ab is not 0.0 a 
                        RuntimeError will be thrown. Default=0.0.

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
                                con_gq=min_control_gq, 
                                con_dp=min_control_dp,
                                con_het_ab=control_het_ab, 
                                con_hom_ab=control_hom_ab, 
                                con_ref_ab=con_ref_ab)

    def filter(self, record, allele):
        '''
            For a given VcfRecord and ALT allele, return True if that 
            allele should be filtered based on presence/absence in cases
            and controls.
        '''
        case_matches = 0
        control_matches = 0
        gts = record.parsed_gts(fields=self.gt_fields, samples=self.samples)
        #check controls first
        for s in self.controls:
            if not self.con_gt_filter.gt_is_ok(gts, s, allele):
                if self.confirm_missing:
                    return True 
                continue
            sgt =  gts['GT'][s]
            if allele in sgt: #checks for presence, not whether het/hom
                if self.n_controls:
                    control_matches += 1
                else:
                    return True
            elif (sgt == (0, 0) and 
                  self.con_gt_filter.ad_over_threshold is not None): 
                #check hom ref for ALT allele counts
                if self.con_gt_filter.ad_over_threshold(gts, s, allele):
                    if self.n_controls:
                        control_matches += 1
                    else:
                        return True
        if self.n_controls and control_matches >= self.n_controls:
            return True
        #check for presence in cases
        for s in self.cases:
            if not self.gt_filter.gt_is_ok(gts, s, allele):
                sgt = None
            else:
                sgt =  gts['GT'][s]
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
                           gq=0, dp=0, het_ab=0., hom_ab=0., con_gq=None, 
                           con_dp=None, con_het_ab=None, con_hom_ab=None,
                           con_ref_ab=None):
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
                               + str.join(", ",s))
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
            raise RuntimeError("The following samples specified by --controls"+
                               " were not found in the input VCF: " + 
                               str.join(", ", s))
        self.cases = list(case_set)
        self.controls = list(control_set)
        self.n_cases = None
        self.n_controls = None
        if n_cases and len(self.cases) < n_cases:
            raise RuntimeError("Number of cases specified by --n_cases is " + 
                               "greater than the number of cases specified by " +
                               "--cases")
        if n_controls and len(self.controls) < n_controls:
            raise RuntimeError("Number of controls specified by --n_controls" + 
                               "is greater than the number of controls " +
                               "specified by --controls")
        self.samples = self.cases + self.controls
        self.gt_filter = GtFilter(self.vcf, gq=gq, dp=dp, het_ab=het_ab, 
                                  hom_ab=hom_ab)
        self.gt_fields = set(self.gt_filter.fields)
        if con_gq is None:
            con_gq = gq
        if con_dp is None:
            con_dp = dp
        if con_het_ab is None:
            con_het_ab = het_ab
        if con_hom_ab is None:
            con_hom_ab = hom_ab
        self.con_gt_filter = GtFilter(self.vcf, gq=con_gq, dp=con_dp, 
                                      het_ab=con_het_ab, hom_ab=hom_ab,
                                      ref_ab_filter=con_ref_ab)
        self.gt_fields.update(self.con_gt_filter.fields)
        if n_cases:
            self.n_cases = n_cases    
        if n_controls:
            self.n_controls = n_controls    
        

class GtFilter(object):
    '''
        Given a dict of GT information from the 'parsed_gts' function of
        VcfRecord from parse_vcf.py, this provides a function 'filter' 
        which returns True if the call meets the criteria (e.g. GQ, AD 
        etc.) established on initialisation.
    '''

    __slots__ = ['gq', 'dp', 'het_ab', 'hom_ab', 'gt_is_ok', 'ab_filter', 
                 'ref_ab_filter', 'ad_over_threshold', 'fields']
    
    def __init__(self, vcf, gq=0, dp=0, het_ab=0., hom_ab=0., 
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
        self.het_ab = het_ab
        self.hom_ab = hom_ab
        self.ref_ab_filter = ref_ab_filter
        self.fields = ['GT']
        self.ab_filter = None
        self.ad_over_threshold = None
        ab_field = None
        if not gq and not dp and not het_ab and not hom_ab: 
            #if no parameters are set then every genotype passes
            self.gt_is_ok = lambda gt, smp, al: True
        else:
            ab_field = self._check_header_fields(vcf)
            if het_ab or hom_ab:
                if ab_field == 'AD':
                    self.ab_filter = self._ab_filter_ad
                elif ab_field == 'RO':
                    self.ab_filter = self._ab_filter_ro
            self.gt_is_ok = self._gt_is_ok
        if ref_ab_filter:
            if ab_field is None:
                ab_field = self._check_header_fields(vcf)
            if ab_field == 'AD':
                self.ad_over_threshold = self._alt_ad_over_threshold
            elif ab_field == 'RO':
                self.ad_over_threshold = self._alt_ao_over_threshold

    def _alt_ad_over_threshold(self, gts, sample, allele):
        ad = gts['AD'][sample]
        if ad == (None,): #no AD values - assume OK?
            return True
        al_dp = ad[allele]
        dp = sum(ad)
        if dp > 0 and al_dp is not None:
            ab = float(al_dp)/dp
            if ab > self.ref_ab_filter:
                #ALT/REF read counts > threshold 
                return True#filter
        return False

    def _alt_ao_over_threshold(self, gts, sample, allele):
        aos = gts['AO'][sample]
        ro = gts['RO'][sample]
        if aos is not None and ro is not None:
            dp = sum(aos) + ro
            if dp > 0:
                ao = aos[allele-1]
                ab = float(ao)/ro
                if ab > self.ref_ab_filter:
                    return True
        return False
        
    def _ab_filter_ad(self, gts, sample, allele):
        ad = gts['AD'][sample]
        if ad == (None,): #no AD values - assume OK?
            return True
        al_dp = ad[allele]
        dp = sum(ad)
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

    def _ab_filter_ro(self, gts, sample, allele):
        aos = gts['AO'][sample]
        ro = gts['RO'][sample]
        is_hom_alt = False
        is_het_alt = False
        if len(set(gts['GT'][sample])) == 1:
            if allele in gts['GT'][sample]:
                is_hom_alt = True
        elif allele in gts['GT'][sample]:
                is_het_alt = True
        if aos is not None and ro is not None and (is_hom_alt or is_het_alt):
            dp = sum(aos) + ro
            if allele > 0:
                ao = aos[allele-1]
            else:
                ao = ro
            if dp > 0:
                ab =float(ao)/dp
                if is_het_alt and ab < self.het_ab:
                    return False #filter
                if is_hom_alt and ab < self.hom_ab:
                    return False #filter
        return True #do not filter

    def _gt_is_ok(self, gts, sample, allele):
        '''
            Returns True if genotype (from parse_vcf.py parsed_gts 
            function) passes all parameters set on initialisation.
        '''
        if self.dp:
            if gts['DP'][sample] is None or gts['DP'][sample] < self.dp:
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
        if self.dp:
            self.fields.append('DP')
        if self.gq:
            self.fields.append('GQ')
        if self.het_ab or self.hom_ab or self.ref_ab_filter:
            if 'AD' in vcf.header.metadata['FORMAT']:
                self.fields.append('AD')
                return 'AD'
            elif ('AO' in vcf.header.metadata['FORMAT'] and 
                  'RO' in vcf.header.metadata['FORMAT']):
                self.fields.append('AO')
                self.fields.append('RO')
                return 'RO'
            else:
                raise RuntimeError("Genotype filtering on allele balance is " +
                                   "set but neither 'AD' nor 'RO' plus 'AO' " +
                                   "FORMAT fields are defined in your VCF " + 
                                   "header.")
        return None

