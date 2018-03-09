from .sample_filter import GtFilter
from collections import defaultdict

class BurdenCounter(object):
    ''' For a set of variants count the number of qualifying alleles
        per transcript.
    '''
    
    def __init__(self, vcf, output, gq=0, dp=0, het_ab=0., hom_ab=0., 
                use_ac=False, gnomad_pop=None, cases=[], controls=[]):
        self.vcf = vcf
        if 'SYMBOL' in vcf.header.csq_fields:
            self.gene_field = 'SYMBOL'
        elif 'Gene' in vcf.header.csq_fields:
            self.gene_field = 'Gene'
        else:
            self.gene_field = 'Feature'
        self.min_gq = gq
        self.min_dp = dp
        self.use_ac = use_ac
        self.samples = None
        self.total_case_alleles = 0
        self.total_control_alleles = 0
        if cases or controls:
            if not vcf.header.samples:
                raise RuntimeError("No samples defined in VCF header - " +
                                   "cannot use --cases or --controls " + 
                                    "arguments.")
            #warn that use_ac is disabled?
            self.use_ac = False
            self.samples = cases + controls
            self.total_case_alleles = len(cases) * 2
            self.total_control_alleles = len(controls) * 2
        self.controls = controls
        self.cases = cases
        self.gt_filter = GtFilter(vcf, gq=gq, dp=dp, het_ab=het_ab, 
                                  hom_ab=hom_ab)
        self.gt_fields = self.gt_filter.fields
        .or x in cases + controls:
            if x not in vcf.header.samples:
               raise RuntimeError("Burden counter sample '{}' not found in "
                                  .format(x) + "VCF input.")
        self.gnomad_pop = None
        if gnomad_pop:
            self._check_gnomad_pop(vcf, gnomad_pop)
            self.gnomad_pop = gnomad_pop
        self.feat_to_cases = defaultdict(set)
        self.feat_to_controls = defaultdict(set)
        self.transcript_to_gene = dict()
        self.case_counts = defaultdict(int)
        self.control_counts = defaultdict(int)
        self.current_features = set()
        self.out_fh = open(output, 'wt')
        self.write_header()

    def write_header(self):
        cols = ["Feature", "Gene", ]
        if self.controls:
            cols.extend(["Controls", "N_Controls"])
            if self.cases
                cols.extend(["Cases", "N_Cases"])
        elif self.gnomad_pop:
            cols.extend([self.gnomad_pop, "N_" + self.gnomad_pop])
        else:
            cols.extend(["Cases", "N_Cases"])
        self.out_fh.write(str.join("\t", cols) + "\n")

    def _check_gnomad_pop(vcf, pop):
        ac = 'AC_' + pop
        an = 'AN_' + pop
        expected = {ac: 'A', an: 1}
        for f in [ac, an]:
            if f in self.vcf.metadata['INFO']:
                num = self.vcf.metadata['INFO'][f][-1]['Number']
                typ = self.vcf.metadata['INFO'][f][-1]['Type']
                if typ != 'Integer':
                    raise RuntimeError("Unexpected 'Type' identifier '{}' "
                                       .format(typ) + "for INFO field {}"
                                       .format(f) + "in input VCF header.\n")
                if num != expected[f]:
                    raise RuntimeError("Unexpected 'Number' identifier '{}' "
                                       .format(num) + "for INFO field {}"
                                       .format(f) + "in input VCF header.\n")
    
    def count(record, ignore_alleles=[], ignore_csq=[]):
        if not self.use_ac and not self.gnomad_pop:
            these_feats = set([x['Feature'] for x in record.CSQ])
            if (self.self.current_features and these_feats.isdisjoint(
                self.current_features)):
                # if we've moved on to next set of features clear feat_to_cases 
                # etc. and add sample counts 
                for feat in self.feat_to_cases:
                    self.case_counts[feat] += len(self.feat_to_cases[feat])
                    del self.feat_to_cases[feat]
                    self.control_counts[feat] += len(self.feat_to_controls[feat])
                    del self.feat_to_controls[feat]
                self.current_features.clear()
            self.current_features.update(these_feats)
        for i in range(len(record.ALLELES) - 1):
            if ignore_alleles and ignore_alleles[i]:
                continue
            for j in range(len(record.CSQ)):
                if ignore_csq and ignore_csq[j]:
                    continue
                if record.CSQ[j]['alt_index'] == allele:
                    feat = record.CSQ[j]['Feature']
                    if not feat:
                        continue
                    self.count_samples(record, feat, i)
                    if feat not in self.transcript_to_gene:
                        gene = record.CSQ[j][self.gene_field]
                        if not gene:
                            gene = feat
                        self.transcript_to_gene = gene

    def count_samples(self, record, feat, allele):
        ''' 
            If using cases and controls add IDs to 
            self.feat_to_cases/controls set, otherwise add number of 
            alleles.
        '''
        if self.gnomad_pop:
            g_ac = 'AC_' + self.gnomad_pop
            g_an = 'AN_' + self.gnomad_pop
            info = record.parsed_info_fields(fields=[g_ac, g_an])
            if info[g_ac][allele] is not None:
                self.case_counts[feat] += info[g_ac][allele]
            if info[g_an][allele] is not None:
                if info[g_an][allele] > self.total_case_alleles:
                    self.total_case_alleles = info[g_an][allele]
        elif self.use_ac:
            info = record.parsed_info_fields(fields=['AC', 'AN'])
            self.case_counts[feat] += info['AC'][allele]
            if info['AN'][allele] is not None:
                if info['AN'][allele] > self.total_case_alleles:
                    self.total_case_alleles = info['AN'][allele]
        else:
            gts = record.parsed_gts(fields=self.gt_fields, samples=self.samples)
            # TODO!
            # major limitation here is that we don't count how many alleles per
            # sample - i.e. if homozygous. Maybe create a dict of samples to 
            # allele counts?
            if cases or controls:
                for s in cases:
                    if self.gt_filter.gt_is_ok(gts, s, allele):
                        if allele in gts['GT'][s]:
                            self.feat_to_cases[feat].add(s)
                for s in controls:
                    if self.gt_filter.gt_is_ok(gts, s, allele):
                        if allele in gts['GT'][s]:
                            self.feat_to_controls[feat].add(s)
            else:
                self.feat_to_cases.update(s for s in gts if 
                  self.gt_filter.gt_is_ok(gts, s, allele) and allele in gts[s])
    
    def output_counts(self):
        for feat in self.case_counts:
            row = [feat, self.transcript_to_gene[feat]]
            if self.controls:
                if feat in self.control_counts:
                    row.append(str(self.control_counts[feat]))
                else:
                    row.append("0")
                row.append(self.total_control_alleles)
            row.append(self.case_counts[feat])
            row.append(self.total_case_alleles)
            self.out_fh.write(str.join("\t", row) + "\n")
        for feat in self.control_counts:
            row.append(str(self.control_counts[feat]))
            row.append(self.total_control_alleles)
            if self.cases:
                if feat in self.case_counts:
                    row.append(self.case_counts[feat])
                else:
                    row.append("0")
                row.append(self.total_case_alleles)
            self.out_fh.write(str.join("\t", row) + "\n")
        self.out_fh.close()
        
