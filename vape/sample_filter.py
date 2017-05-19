

class SampleFilter(object):
    ''' A class for filtering VCF records on sample calls'''

    def __init__(self, args, vcf):
        '''
            Initialize with parsed args from vape.py and input VcfReader
            object
        '''
        self.cases = []
        self.controls = []
        self.samples = []
        self.n_cases = 0
        self.n_controls = 0
        self.min_gq = 0
        self.vcf = vcf
        self._parse_sample_args(args)

    def filter(self, record):
        '''
            For a given VcfRecord, return a list of booleans indicating
            for each allele whether it should be filtered or not.
        '''
        filter_alleles = [False] * (len(record.ALLELES) -1)
        case_allele_matches = [0] * (len(record.ALLELES) -1) 
        control_allele_matches = [0] * (len(record.ALLELES) -1) 
        gt_fields = ['GT']
        if self.min_gq:
            gt_fields.append('GQ')
        gts = record.parsed_gts(fields=gt_fields, samples=self.samples)
        #check controls first
        for s in self.controls:
            if self.min_gq and (gts['GQ'][s] is None 
                or gts['GQ'][s] < self.min_gq):
                continue
            sgt =  gts['GT'][s]
            for i in range(1, len(record.ALLELES)):
                if i in sgt: #checks for presence, not whether het/hom
                    if not self.n_controls:
                        filter_alleles[i-1] = True
                    control_allele_matches[i-1] += 1
            #bail out if all ALTs are present in controls and !self.n_controls
            if sum(filter_alleles) == len(filter_alleles):
                return filter_alleles
        
        if self.n_controls:
            for i in range(1, len(record.ALLELES)):
                if control_allele_matches[i-1] > self.n_controls:
                    filter_alleles[i-1] = True
            if sum(filter_alleles) == len(filter_alleles):
                return filter_alleles

        #check for presence in cases
        for s in self.cases:
            if self.min_gq and (gts['GQ'][s] is None 
                or gts['GQ'][s] < self.min_gq):
                #if GQ is None presumably is a no call
                sgt = None
            else:
                sgt =  gts['GT'][s]
            if sgt is None:
                if not self.n_cases:
                    return [True] * (len(record.ALLELES) -1)
                else:
                    continue
            for i in range(1, len(record.ALLELES)):
                if i in sgt:
                    case_allele_matches[i-1] += 1
                elif not self.n_cases:
                    filter_alleles[i-1] = True
            #bail if all ALTs are absent in at least one case and !self.n_cases
            if sum(filter_alleles) == len(filter_alleles):
                return filter_alleles
        if self.n_cases:
            for i in range(1, len(record.ALLELES)):
                if case_allele_matches[i-1] < self.n_cases:
                    filter_alleles[i-1] = True
            
        return filter_alleles
            
                    
            

    def _parse_sample_args(self, args):
        not_found = set()
        cases = set()
        controls = set()
        for c in args.cases:
            if c in self.vcf.header.samples:
                cases.add(c)
            elif c.lower() == 'all':
                if len(args.cases) > 1:
                    raise Exception("'all' can not be used in --cases " + 
                                    "argument in conjunction with sample " + 
                                    "names - please specify either 'all' or " +
                                    "a list of sample IDs, not both.")
                cases.update(self.vcf.header.samples)
            else:
                not_found.add(c)
        if not_found:
            s = list(not_found)
            s.sort()
            raise Exception("The following samples specified by --cases were" +
                            " not found in the input VCF: " + str.join(", ",s))
        for c in args.controls:
            if c in self.vcf.header.samples:
                controls.add(c)
            elif c.lower() == 'all':
                if len(args.controls) > 1:
                    raise Exception("'all' can not be used in --controls " + 
                                    "argument in conjunction with sample " + 
                                    "names - please specify either 'all' or " +
                                    "a list of sample IDs, not both.")
                for samp in self.vcf.header.samples:
                    if samp not in cases:
                        controls.add(samp)
            else:
                not_found.add(c)
        if not_found:
            s = list(not_found)
            s.sort()
            raise Exception("The following samples specified by --controls " +
                            "were not found in the input VCF: " + 
                            str.join(", ", s))
        self.cases = list(cases)
        self.controls = list(controls)
        if args.n_cases and len(self.cases) < args.n_cases:
            raise Exception("Number of cases specified by --n_cases is " + 
                            "greater than the number of cases specified by " +
                            "--cases")
        if args.n_controls and len(self.controls) < args.n_controls:
            raise Exception("Number of controls specified by --n_controls " + 
                            "is greater than the number of controls " +
                            "specified by --controls")
        self.samples = self.cases + self.controls
        if args.gq:
            self.min_gq = args.gq
        if args.n_cases:
            self.n_cases = args.n_cases    
        if args.n_controls:
            self.n_controls = args.n_controls    
        


