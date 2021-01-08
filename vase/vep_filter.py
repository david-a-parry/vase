import os
import logging
from .insilico_filter import InSilicoFilter
from .csq_filter import CsqFilter

lof_csq = {'frameshift_variant', 'stop_gained', 'splice_acceptor_variant',
           'splice_donor_variant'}


class VepFilter(CsqFilter):
    '''An object that filters VCF records based on annotated VEP data.'''

    def __init__(self, vcf, csq=[], impact=[], canonical=False, biotypes=[],
                 in_silico=[], filter_unpredicted=False,
                 keep_any_damaging=False, splice_in_silico=[],
                 loftee=False, splice_filter_unpredicted=False,
                 splice_keep_any_damaging=False, retain_labels=[],
                 filter_flagged_features=False, freq=None, min_freq=None,
                 afs=[], gene_filter=None, blacklist=None,filter_known=False,
                 filter_novel=False, pathogenic=False, no_conflicted=False,
                 g2p=None, check_g2p_consequence=False,
                 logging_level=logging.WARNING):
        '''
            Args:
                vcf:    input VcfReader object

                csq:    list of consequence types to keep. If 'default'
                        appears anywhere in this list then the default
                        consequence set (as indicated in
                        data/vep_classes.tsv) will be used. Similarly if
                        'all' appears anywhere in this list no filtering
                        on consequence type will occur.

                impact: list of variant impacts to retain.

                canonical:
                        Filter consequences on non-canonical transcirpts.

                biotypes:
                        Filter consequences for features not of the
                        given biotypes. If not provided the default set
                        of biotypes (as indicated in data/biotypes.tsv)
                        will be used for biotype filtering.

                in_silico:
                        List of programs and optionally score criteria
                        for filtering of missense variants using the
                        InSilicoFilter class in vase.insilico_filter.

                filter_unpredicted:
                        If using 'in_silico' option, filter missense
                        variants that have missing values for any of the
                        specified filtering programs.

                keep_any_damaging:
                        If using 'in_silico' option, retain variants if
                        any of the criteria are met for any of the
                        specified filtering programs.

                loftee: Only retain LoF (i.e. high impact variants)
                        variants if the LoF annotation from loftee is
                        'HC' (high confidence).

                splice_in_silico:
                        Similar to 'in_silico' but the prediction
                        programs are checked for splice_donor_variants,
                        splice_acceptor_variants and
                        splice_region_variants rather than missense.
                        Currently only dbscSNV (rf_score and ada_score),
                        MaxEntScan and SpliceDistance
                        (https://github.com/david-a-parry/SpliceDistance)
                        annotations are supported. This option can be
                        used to, for example, retain
                        splice region variants that are have
                        an 'ada_score' > 0.6 by specifying
                        'ada_score=0.6' with this option.

                splice_filter_unpredicted:
                        If using 'splice_in_silico' option, filter
                        splice region variants that have missing values
                        for any of the specified filtering programs.

                splice_keep_any_damaging:
                        If using 'splice_in_silico' option, retain
                        variants if any of the criteria are met for any
                        of the specified filtering programs.

                retain_labels:
                        Do not filter on consequence type if the
                        following values are present for a label. Labels
                        and values must be separated by '=' sign. For
                        example, to retain any consequence which has
                        a VEP annotation named 'FOO' with  value 'BAR'
                        use 'FOO=BAR'.

                filter_flagged_features:
                        Filter consequences on features which are
                        flagged by VEP.

                freq:   Filter consequences if the annotated allele
                        frequency is equal to or greater than this value.
                        By default all allele frequency annotations as
                        listed in "data/vep_maf.tsv" are used, but this
                        can be altered using the 'afs' option.

                min_freq:
                        As for 'freq' argument but filters consequences
                        if the allele frequency annotation is less than
                        this value.

                filter_known:
                        Filter consequences if allele frequency is given
                        for any of the available VEP frequency
                        annotations.

                filter_novel:
                        Filter consequences if no allele frequency is
                        given for any of the available VEP frequency
                        annotations.

                afs:    Only use the listed allele frequency annotations
                        for freq/min_freq/novelty filtering.

                gene_filter:
                        VarByRegion object from vase.var_by_region. If
                        provided, consequences will be filtered if they
                        do not alter the features specified in the
                        VarByRegion object for the current region.

                blacklist:
                        File containing a list of Feature IDs to ignore.

                pathogenic:
                        If True, retain consequences regardless of type
                        if annotated as 'pathogenic' or 'likely
                        pathogenic' in 'CLIN_SIG' or 'clinvar_clnsig'
                        VEP fields. Frequency, biotype and canonical
                        filtering will still be applied.

                no_conflicted:
                        If 'pathogenic' option is True, only retain
                        'likely pathogenic' and 'pathogenic'
                        consequences if there are no conflicting
                        'benign' or 'likely benign' assertions.

                g2p:
                        G2P object from vase.g2p for filtering on
                        presence and/or requirements from a G2P file.

                check_g2p_consequence:
                        If a G2P object is provided above, require that
                        that the observed consequence matches the
                        'mutation consequence' in the G2P file.

                logging_level:
                        Logging level to use. Default=logging.WARNING.

        '''
        self.logger = self._get_logger(logging_level)
        self.canonical = canonical
        self.loftee = loftee
        self.filter_flagged = filter_flagged_features
        self.freq = freq
        self.min_freq = min_freq
        self.afs = afs
        self.filter_known = filter_known
        self.filter_novel = filter_novel
        self._check_freq_fields(vcf)
        self.in_silico = False
        self.splice_in_silico = False
        if in_silico:
            in_silico = set(in_silico)
            self.in_silico = InSilicoFilter(in_silico, filter_unpredicted,
                                            keep_any_damaging)
        if splice_in_silico:
            splice_in_silico = set(splice_in_silico)
            self.splice_in_silico = InSilicoFilter(
                        programs=splice_in_silico,
                        filter_unpredicted=splice_filter_unpredicted,
                        keep_if_any_damaging=splice_keep_any_damaging,
                        pred_file=os.path.join(os.path.dirname(__file__),
                                               "data",
                                               "vep_splice_insilico_pred.tsv"))
        self.pathogenic = pathogenic
        self.no_conflicted = no_conflicted
        if pathogenic:
            self.path_fields = self._get_path_fields(vcf)
        super().__init__(vcf=vcf, csq_attribute='CSQ', csq=csq, impact=impact,
                         biotypes=biotypes, retain_labels=retain_labels,
                         filter_flagged_features=filter_flagged_features,
                         gene_filter=gene_filter, blacklist=blacklist, g2p=g2p,
                         check_g2p_consequence=check_g2p_consequence)

    def filter_csq(self, csq):
        '''
        Returns two boolean values. The first indicates whether the consequence
        annotation should be filtered. The second indicates whether the ALT
        allele should be filtered irrespective of the given or any other
        consequence annotation.
        '''
        if self.canonical:
            try:
                if csq['CANONICAL'] != 'YES':
                    return True, False
            except KeyError:
                pass
        if self.filter_flagged:
            try:
                if csq['FLAGS']:
                    return True, False
            except KeyError:
                pass
        if (self.biotypes is not None and csq['BIOTYPE'].lower() not in
                self.biotypes):
            return True, False
        if self.gene_filter:
            if not self.gene_filter.target_in_csq(csq):
                return True, False
        if self.g2p:
            if csq['SYMBOL'] not in self.g2p.g2p:
                return True, False
        if self.blacklist and csq['Feature'] in self.blacklist:
            return True, False
        if (self.freq or self.min_freq or self.filter_known or
                self.filter_novel):
            known = False
            for af in self.freq_fields:
                if csq[af] == '' or csq[af] == '.':
                    continue
                try:
                    c_af = float(csq[af])
                except ValueError:
                    try:
                        c_af = max(float(x) for x in csq[af].split('&') if x
                                   != '.')
                    except ValueError:
                        continue
                known = True
                if self.filter_known:
                    return True, True
                if self.freq:
                    if c_af >= self.freq:
                        return True, True
                if self.min_freq:
                    if c_af < self.min_freq:
                        return True, True
            if self.filter_novel and not known:
                return True, True
        if (self.csq is None and self.impact is None and
                not self.check_g2p_consequence):
            # if only using biotypes/MAF for filtering
            return False, False
        if self.pathogenic and self._has_pathogenic_annotation(csq):
            return False, False
        if self._retain_label_matched(csq):
            return False, False
        if self.check_g2p_consequence and self.g2p:
            filt_csq = self.g2p.consequences_from_gene(csq['SYMBOL'])
        else:
            filt_csq = self.csq
        for s_csq in [x.lower() for x in csq['Consequence'].split('&')]:
            matches_csq = False
            matches_impact = False
            if filt_csq is not None and s_csq in filt_csq:
                matches_csq = True
            if self.impact is not None and csq['IMPACT'] in self.impact:
                matches_impact = True
            if matches_csq or matches_impact:
                if self.in_silico and s_csq == 'missense_variant':
                    do_filter = self.in_silico.filter(csq)
                    if not do_filter:
                        return False, False
                elif self.splice_in_silico and s_csq.startswith("splice"):
                    do_filter = self.splice_in_silico.filter(csq)
                    if not do_filter:
                        return False, False
                elif self.loftee and (s_csq in lof_csq or matches_impact
                                      and csq['IMPACT'] == 'HIGH'):
                    if csq['LoF'] == 'HC':
                        return False, False
                else:
                    return False, False
        return True, False

    def _has_pathogenic_annotation(self, csq):
        path = []
        benign = []
        for annot in self.path_fields:
            if not csq[annot]:
                continue
            assertions = csq[annot].split('&')
            if annot == 'clinvar_clnsig':
                # benign = 2, likely benign = 3
                # likely pathogenic = 4, pathogenic = 5
                try:
                    benign.extend((4 > int(x) > 1 for x in assertions))
                    path.extend((6 > int(x) > 3 for x in assertions))
                except ValueError:
                    self.logger.warn("Error parsing 'clinvar_clnsig' field " +
                                     "'{}' - expected numeric values.".format(
                                         csq[annot]))
            else:
                benign.extend(('benign' in x for x in assertions))
                path.extend(('pathogenic' in x for x in assertions))
        if self.no_conflicted:
            return any(path) and not any(benign)
        return any(path)

    def _read_maf_file(self):
        data_file = os.path.join(os.path.dirname(__file__),
                                 "data",
                                 "vep_maf.tsv")
        values = []
        with open(data_file, encoding='UTF-8') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                values.append(line.rstrip())
        return values

    def _check_freq_fields(self, vcf):
        self.freq_fields = []
        if (not self.freq and not self.min_freq and not self.filter_novel and
                not self.filter_known):
            return
        if self.afs:
            for fq in self.afs:
                if fq in vcf.header.csq_fields:
                    self.freq_fields.append(fq)
                    self.logger.info("Found '{}' VEP allele ".format(fq) +
                                     "frequency annotation")
                else:
                    raise RuntimeError("Could not find '{}' ".format(fq) +
                                       "VEP AF field in VEP annotations.")
        else:
            for fq in self._read_maf_file():
                if fq in vcf.header.csq_fields:
                    self.freq_fields.append(fq)
                    self.logger.info("Found '{}' VEP allele ".format(fq) +
                                     "frequency annotation")
        if not self.freq_fields:
            self.logger.warn("No compatible (>= v90) allele frequency fields" +
                             " in VEP annotations.")

    def _get_path_fields(self, vcf):
        cln_fields = ['CLIN_SIG', 'clinvar_clnsig']
        path_fields = [f for f in vcf.header.csq_fields if f in cln_fields]
        if not path_fields:
            self.logger.warn("No compatible ClinVar VEP annotations found " +
                             "for use with pathogenic allele identification.")
        return path_fields

    def get_required_header_fields(self):
        '''
        Check which CSQ/ANN annotation fields are required given arguments
        passed to __init__
        '''
        required = ['Consequence', 'BIOTYPE']
        if self.impact:
            required.append('IMPACT')
        if self.canonical:
            required.append('CANONICAL')
        if self.loftee:
            required.append('LoF')
        if self.filter_flagged:
            required.append('FLAGS')
        return required
