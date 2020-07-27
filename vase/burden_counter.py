import re
from .sample_filter import GtFilter
from collections import defaultdict


class BurdenCounter(object):
    ''' For a set of variants count the number of qualifying alleles
        per transcript.
    '''

    def __init__(self, vcf, output, gq=0, dp=0, max_dp=0, het_ab=0., hom_ab=0.,
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
        self.max_dp = max_dp
        self.samples = None
        self.use_ac = False
        self.total_alleles = {'Cases': 0, 'Controls': 0}
        self.controls = controls
        self.cases = cases
        vcf_samples = list(vcf.header.samples)
        if cases or controls:
            if not vcf_samples:
                raise RuntimeError("No samples defined in VCF header - " +
                                   "cannot use --cases or --controls " +
                                   "arguments.")
            self.samples = cases + controls
            self.total_alleles['Cases'] = len(cases) * 2
            self.total_alleles['Controls'] = len(controls) * 2
        elif vcf_samples:
            self.cases = vcf_samples
            self.total_alleles['Cases'] = len(self.cases) * 2
            self.samples = self.cases
        elif not is_gnomad:
            self.use_ac = True
        for x in cases + controls:
            if x not in vcf.header.samples:
                raise ValueError("Burden counter sample '{}' not found in "
                                 .format(x) + "VCF input.")
        self.gnomad_pops = []
        if is_gnomad:
            self.gnomad_pops = self._check_gnomad_pops(vcf)
            self.total_alleles.clear()
            for pop in self.gnomad_pops:
                self.total_alleles[pop] = 0
        elif not self.use_ac:
            self.gt_filter = GtFilter(vcf, gq=gq, dp=dp, max_dp=max_dp,
                                      het_ab=het_ab, hom_ab=hom_ab)
        self.feat_to_cases = defaultdict(dict)
        self.feat_to_controls = defaultdict(dict)
        self.transcript_to_gene = dict()
        self.counts = defaultdict(dict)
        self.current_features = set()
        self.out_fh = open(output, 'wt')
        self.write_header()

    def write_header(self):
        cols = ["Feature", "Gene"]
        if self.gnomad_pops:
            for pop in self.gnomad_pops:
                cols.extend([pop, "N_" + pop])
        elif self.cases or self.controls:
            if self.cases:
                cols.extend(["Cases", "N_Cases"])
            if self.controls:
                cols.extend(["Controls", "N_Controls"])
        else:
            # all samples labelled as 'Cases' if cases and controls are empty
            cols.extend(["Cases", "N_Cases"])
        self.out_fh.write(str.join("\t", cols) + "\n")

    def _check_gnomad_pops(self, vcf):
        pop_ac_re = re.compile(r'''^AC_([A-Za-z]{3})$''')
        pops = []
        for f in self.vcf.header.info:
            match = pop_ac_re.match(f)
            if match:
                p = match.group(1)
                an = 'AN_' + p
                ac_num = self.vcf.header.info[f].number
                ac_typ = self.vcf.header.info[f].type
                an_num = self.vcf.header.info[an].number
                an_typ = self.vcf.header.info[an].type
                if (ac_num == 'A' and ac_typ == 'Integer' and
                        an_num == '1' and an_typ == 'Integer'):
                    pops.append(p)
        if not pops:
            raise RuntimeError("No gnomAD populations found for VCF input!")
        pops.sort()
        return pops

    def count(self, record, ignore_alleles=[], ignore_csq=[]):
        '''
            If using cases and controls add IDs to
            self.feat_to_cases/controls set, otherwise add number of
            alleles.

            Args:
                record: VcfRecord from VcfReader

                ignore_alleles:
                        List of 0-based ALT allele indices that should
                        NOT be counted.

                ignore_csq:
                        List of 0-based CSQ annotations from VEP that
                        should not be counted (i.e. do not assign
                        counts for the 'features' associated with
                        these consequences).

        '''
        if not self.use_ac and not self.gnomad_pops:
            these_feats = set([x['Feature'] for x in record.CSQ])
            if (self.current_features and these_feats.isdisjoint(
                    self.current_features)):
                # if we've moved on to next set of features clear feat_to_cases
                # etc. and add sample counts
                feats_to_delete = set()
                # can't delete keys while iterating over dict
                for feat in self.feat_to_cases:
                    feats_to_delete.add(feat)
                    self.counts[feat]['Cases'] = sum(
                            self.feat_to_cases[feat][x] for x in
                            self.feat_to_cases[feat])
                    self.counts[feat]['Controls'] = sum(
                            self.feat_to_controls[feat][x] for x in
                            self.feat_to_controls[feat])
                for feat in (x for x in self.feat_to_controls if x not in
                             self.feat_to_cases):
                    self.counts[feat]['Controls'] = sum(
                            self.feat_to_controls[feat][x] for x in
                            self.feat_to_controls[feat])
                    feats_to_delete.add(feat)
                for feat in feats_to_delete:
                    self.feat_to_cases.pop(feat, None)
                    self.feat_to_controls.pop(feat, None)
                self.current_features.clear()
            self.current_features.update(these_feats)
        for i in range(len(record.alts)):
            features = []
            if ignore_alleles and ignore_alleles[i]:
                continue
            for j in range(len(record.CSQ)):
                if ignore_csq and ignore_csq[j]:
                    continue
                if record.CSQ[j]['alt_index'] == i + 1:
                    features.append(record.CSQ[j]['Feature'])
            if not features:
                continue
            self.count_samples(record, features, i)

    def count_samples(self, record, features, allele, max_alleles=2):
        '''
            If using cases and controls add IDs to
            self.feat_to_cases/controls set, otherwise add number of
            alleles.

            Args:
                record: VcfRecord from VcfReader

                features:
                        Names of the features/transcripts for
                        registering allele counts.

                allele: 0-based ALT allele index (i.e. 0 = first ALT
                        allele)

                max_alleles:
                        Maximum number of allele counts that can be
                        assigned to one sample. Default=2.

        '''
        a_counts = defaultdict(int)
        if self.gnomad_pops:
            for pop in self.gnomad_pops:
                g_ac = 'AC_' + pop
                g_an = 'AN_' + pop
                if record.info[g_ac][allele] is not None:
                    a_counts[pop] += record.info[g_ac][allele]
                if record.info[g_an] is not None:
                    if record.info[g_an] > self.total_alleles[pop]:
                        self.total_alleles[pop] = record.info[g_an]
        elif self.use_ac:
            a_counts['Cases'] += record.info['AC'][allele]
            if record.info['AN'] is not None:
                if record.info['AN'] > self.total_alleles['Cases']:
                    self.total_alleles['Cases'] = record.info['AN']
        else:
            a_counts = dict((s, record.samples[s]['GT'].count(allele+1)) for s
                            in self.samples if self.gt_filter.gt_is_ok(
                                record.samples, s, allele) and (allele+1) in
                            record.samples[s]['GT'])
        for feat in features:
            if not feat:  # skip any intergenic variants
                continue
            self._check_gene_name(feat, record)
            if self.gnomad_pops or self.use_ac:
                for group in a_counts:
                    if group in self.counts[feat]:
                        self.counts[feat][group] += a_counts[group]
                    else:
                        self.counts[feat][group] = a_counts[group]
            else:  # if we have samples we can ensure we don't count twice
                if self.cases or self.controls:
                    for s in self.cases:
                        if s in a_counts:
                            if s in self.feat_to_cases[feat]:
                                self.feat_to_cases[feat][s] += a_counts[s]
                                # do not count more than 2 alleles per sample
                                if self.feat_to_cases[feat][s] > max_alleles:
                                    self.feat_to_cases[feat][s] = max_alleles
                            else:
                                self.feat_to_cases[feat][s] = a_counts[s]
                    for s in self.controls:
                        if s in a_counts:
                            if s in self.feat_to_cases[feat]:
                                self.feat_to_controls[feat][s] += a_counts[s]
                                if self.feat_to_controls[feat][s] > max_alleles:
                                    # do not count more than 2 alleles per samp
                                    self.feat_to_controls[feat][s] = max_alleles
                            else:
                                self.feat_to_controls[feat][s] = a_counts[s]
                else:
                    self.feat_to_controls[feat] = a_counts

    def output_counts(self):
        if not self.use_ac and not self.gnomad_pops:
            for feat in self.feat_to_cases:
                self.counts[feat]['Cases'] = sum(
                        self.feat_to_cases[feat][x] for x in
                        self.feat_to_cases[feat])
                self.counts[feat]['Controls'] = sum(
                        self.feat_to_controls[feat][x] for x in
                        self.feat_to_controls[feat])
            for feat in (x for x in self.feat_to_controls if x not in
                         self.feat_to_cases):
                self.counts[feat]['Controls'] = sum(
                        self.feat_to_controls[feat][x] for x in
                        self.feat_to_controls[feat])
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
                    row.append(str(self.total_alleles[g] -
                                   self.counts[feat][g]))
                else:
                    row.append("0")
                    row.append(str(self.total_alleles[g]))
            self.out_fh.write(str.join("\t", row) + "\n")
        self.out_fh.close()

    def _check_gene_name(self, feat, record):
        if feat not in self.transcript_to_gene:
            csq = next(x for x in record.CSQ if x['Feature'] == feat)
            gene = csq[self.gene_field]
            if not gene:
                gene = feat
            self.transcript_to_gene[feat] = gene
