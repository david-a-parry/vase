

class SampleFilter(object):
    ''' A class for filtering VCF records on sample calls'''

    def __init__(self, vcf, cases=[], controls=[], n_cases=0, n_controls=0, 
                 gq=0, confirm_missing=False):
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

                confirm_missing:
                        If True, only keep a variant if all controls are
                        not no-calls and are above the GQ threshold set
                        by the 'gq' option. Default=False.

        '''
                        
        self.vcf = vcf
        self.confirm_missing = confirm_missing
        self._parse_sample_args(cases, controls, n_cases, n_controls, gq)

    def filter(self, record, allele):
        '''
            For a given VcfRecord and ALT allele, return True if that 
            allele should be filtered based on presence/absence in cases
            and controls.
        '''
        case_matches = 0
        control_matches = 0
        gt_fields = ['GT']
        if self.min_gq:
            gt_fields.append('GQ')
        gts = record.parsed_gts(fields=gt_fields, samples=self.samples)
        #check controls first
        for s in self.controls:
            if self.min_gq and (gts['GQ'][s] is None 
                or gts['GQ'][s] < self.min_gq):
                if self.confirm_missing:
                    return True 
                continue
            sgt =  gts['GT'][s]
            if allele in sgt: #checks for presence, not whether het/hom
                if self.n_controls:
                    control_matches += 1
                else:
                    return True
        
        if self.n_controls and control_matches >= self.n_controls:
            return True
        #check for presence in cases
        for s in self.cases:
            if self.min_gq and (gts['GQ'][s] is None 
                or gts['GQ'][s] < self.min_gq):
                #if GQ is None presumably is a no call
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
                           gq=0):
        not_found = set()
        case_set = set()
        control_set = set()
        for c in cases:
            if c in self.vcf.header.samples:
                case_set.add(c)
            elif c.lower() == 'all':
                if len(cases) > 1:
                    raise Exception("'all' can not be used in --cases " + 
                                    "argument in conjunction with sample " + 
                                    "names - please specify either 'all' or " +
                                    "a list of sample IDs, not both.")
                case_set.update(self.vcf.header.samples)
            else:
                not_found.add(c)
        if not_found:
            s = list(not_found)
            s.sort()
            raise Exception("The following samples specified by --cases were" +
                            " not found in the input VCF: " + str.join(", ",s))
        for c in controls:
            if c in self.vcf.header.samples:
                control_set.add(c)
            elif c.lower() == 'all':
                if len(controls) > 1:
                    raise Exception("'all' can not be used in --controls " + 
                                    "argument in conjunction with sample " + 
                                    "names - please specify either 'all' or " +
                                    "a list of sample IDs, not both.")
                for samp in self.vcf.header.samples:
                    if samp not in case_set:
                        control_set.add(samp)
            else:
                not_found.add(c)
        if not_found:
            s = list(not_found)
            s.sort()
            raise Exception("The following samples specified by --controls " +
                            "were not found in the input VCF: " + 
                            str.join(", ", s))
        self.cases = list(case_set)
        self.controls = list(control_set)
        self.n_cases = None
        self.n_controls = None
        if n_cases and len(self.cases) < n_cases:
            raise Exception("Number of cases specified by --n_cases is " + 
                            "greater than the number of cases specified by " +
                            "--cases")
        if n_controls and len(self.controls) < n_controls:
            raise Exception("Number of controls specified by --n_controls " + 
                            "is greater than the number of controls " +
                            "specified by --controls")
        self.samples = self.cases + self.controls
        if gq:
            self.min_gq = gq
        if n_cases:
            self.n_cases = n_cases    
        if n_controls:
            self.n_controls = n_controls    
        

