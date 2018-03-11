from .sample_filter import GtFilter
from collections import defaultdict

class BurdenCounter(object):
    ''' For a set of variants count the number of qualifying alleles
        per transcript.
    '''
    
    def __init__(self, vcf, output, gq=0, dp=0, het_ab=0., hom_ab=0., 
                 is_gnomad=False, cases=[], controls=[]):
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
        self.use_ac = False
        self.total_alleles = {'Cases': 0, 'Controls': 0}
        if cases or controls:
            if not vcf.header.samples:
                raise RuntimeError("No samples defined in VCF header - " +
                                   "cannot use --cases or --controls " + 
                                    "arguments.")
            self.samples = cases + controls
            self.total_alleles['Cases'] = len(cases) * 2
            self.total_alleles['Controls'] = len(controls) * 2
        elif not vcf.header.samples and not is_gnomad:
            self.use_ac = True
        self.controls = controls
        self.cases = cases
        self.gt_fields = self.gt_filter.fields
        for x in cases + controls:
            if x not in vcf.header.samples:
               raise RuntimeError("Burden counter sample '{}' not found in "
                                  .format(x) + "VCF input.")
        self.gnomad_pops = []
        if is_gnomad:
            self.gnomad_pops = self._check_gnomad_pops(vcf)
            for pop in self.gnomad_pops:
                self.total_alleles.clear()
                self.total_alleles[pop] = 0
        elif not self.use_ac:
            self.gt_filter = GtFilter(vcf, gq=gq, dp=dp, het_ab=het_ab, 
                                      hom_ab=hom_ab)
        self.feat_to_cases = defaultdict(dict)
        self.feat_to_controls = defaultdict(dict)
        self.transcript_to_gene = dict()
        self.counts = defaultdict(dict)
        self.current_features = set()
        self.out_fh = open(output, 'wt')
        self.write_header()

    def write_header(self):
        cols = ["Feature", "Gene", ]
        if self.gnomad_pops:
            for pop in self.gnomad_pops:
                cols.extend([pop, "N_" + pop])
        if self.cases or self.controls:
            if self.cases:
                cols.extend(["Cases", "N_Cases"])
            if self.controls:
                cols.extend(["Controls", "N_Controls"])
        else: 
            #all samples labelled as 'Cases' if neither cases or controls given 
            cols.extend(["Cases", "N_Cases"])
        self.out_fh.write(str.join("\t", cols) + "\n")

    def _check_gnomad_pop(vcf, pop):
        ac = 'AC_' + pop
        an = 'AN_' + pop
        pop_ac_re = re.compile(r'''^AC_([A-Z]+)$''')
        expected = {ac: 'A', an: 1}
        pops = []
        for f in self.vcf.metadata['INFO']:
            match = pop_ac_re.match(f)
            if match:
                p = match.groups(1)
                an = 'AN_' + p
                ac_num = self.vcf.metadata['INFO'][f][-1]['Number']
                ac_typ = self.vcf.metadata['INFO'][f][-1]['Type']
                an_num = self.vcf.metadata['INFO'][an][-1]['Number']
                an_typ = self.vcf.metadata['INFO'][an][-1]['Type']
                if (ac_num == 'A' and ac_type == 'Integer' and 
                    an_num == 1 and ac_type == 'Integer'):
                    pops.append(p)
        if not pops:
            raise RuntimeError("No gnomAD populations found for VCF input!")
        pops.sort()
        return pops
    
    def count(record, ignore_alleles=[], ignore_csq=[]):
        if not self.use_ac and not self.gnomad_pops:
            these_feats = set([x['Feature'] for x in record.CSQ])
            if (self.self.current_features and these_feats.isdisjoint(
                self.current_features)):
                # if we've moved on to next set of features clear feat_to_cases 
                # etc. and add sample counts 
                for feat in self.feat_to_cases:
                    self.counts[feat]['Cases'] = sum(
                        self.feat_to_cases[feat][x] for x in 
                        self.feat_to_cases[feat])
                    del self.feat_to_cases[feat]
                    self.counts[feat]['Controls'] = sum(
                        self.feat_to_controls[feat][x] for x in 
                        self.feat_to_controls[feat])
                    del self.feat_to_controls[feat]
                for feat in (x for x in self.feat_to_controls if x not in
                             self.feat_to_cases):
                    self.counts[feat]['Controls'] = sum(
                        self.feat_to_controls[feat][x] for x in 
                        self.feat_to_controls[feat])
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
        if self.gnomad_pops:
            for pop in self.gnomad_pops:
                g_ac = 'AC_' + self.gnomad_pop
                g_an = 'AN_' + self.gnomad_pop
                info = record.parsed_info_fields(fields=[g_ac, g_an])
                if info[g_ac][allele] is not None:
                    self.counts[feat][pop] += info[g_ac][allele]
                if info[g_an] is not None:
                    if info[g_an] > self.total_alleles[pop]:
                        self.total_alleles[pop] = info[g_an]
        elif self.use_ac:
            info = record.parsed_info_fields(fields=['AC', 'AN'])
            self.counts[feat]['Cases'] += info['AC'][allele]
            if info['AN'] is not None:
                if info['AN'] > self.total_alleles['Cases']:
                    self.total_alleles['Cases'] = info['AN']
        else:
            gts = record.parsed_gts(fields=self.gt_fields, samples=self.samples)
            if cases or controls:
                for s in cases:
                    if self.gt_filter.gt_is_ok(gts, s, allele):
                        ac = gts['GT'][s].count(allele)
                        if s in self.feat_to_cases[feat]:
                            self.feat_to_cases[feat][s] += ac
                            #do not count more than 2 alleles per sample
                            if self.feat_to_cases[feat][s] > 2:
                                self.feat_to_cases[feat][s] = 2
                        else:
                            self.feat_to_cases[feat][s] = ac
                for s in controls:
                    if self.gt_filter.gt_is_ok(gts, s, allele):
                        ac = gts['GT'][s].count(allele)
                        if s in self.feat_to_cases[feat]:
                            self.feat_to_controls[feat][s] += ac
                            if self.feat_to_controls[feat][s] > 2:
                                #do not count more than 2 alleles per sample
                                self.feat_to_controls[feat][s] = 2
                        else:
                            self.feat_to_controls[feat][s] = ac
            else:
                self.feat_to_cases.update(gts['GT'][s].count(allele) for s in 
                                          gts if self.gt_filter.gt_is_ok(
                                          gts, s, allele) and allele in gts[s])

    def output_counts(self):
        groups = []
        if self.gnomad_pops:
            groups = self.gnomad_pops
        elif self.use_ac:
            groups = ['Cases']
        elif self.cases or self.controls:
            if self.cases:
                groups.append('Cases')
            if self.controls:
                groups.append('Controls')
        else:
            groups = ['Cases']
        for feat in self.counts:
            row = [feat, self.transcript_to_gene[feat]]
            for g in groups:
                if g in self.counts[feat]:
                    row.append(str(self.counts[feat][g]))
                else:
                    row.append("0")
                row.append(self.total_alleles[g])
            self.out_fh.write(str.join("\t", row) + "\n")
        self.out_fh.close()
        
