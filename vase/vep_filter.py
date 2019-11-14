import os
import logging
from collections import defaultdict
from parse_vcf import HeaderError
from .insilico_filter import InSilicoFilter

lof_csq = {'frameshift_variant', 'stop_gained', 'splice_acceptor_variant',
           'splice_donor_variant'}

class VepFilter(object):
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
        default_csq, valid_csq = self._read_csq_file()
        default_biotypes, valid_biotypes = self._read_biotype_file()
        self.csq = set()
        self.impact = None
        self.biotypes = set()
        if not csq and not impact:
            csq = ['default']
        if csq is None:
            csq = []
        for c in csq:
            lc = c.lower()
            if lc == 'default':
                self.csq.update(default_csq)
            elif lc == 'all':
                self.csq = None
                break
            else:
                if lc in valid_csq:
                    self.csq.add(lc)
                else:
                    raise RuntimeError("ERROR: Unrecognised VEP consequence " +
                                       "class '{}'".format(c))
        if impact:
            self.impact = set((x.upper() for x in impact))
            valid_impacts = set(['HIGH', 'MODERATE', 'LOW', 'MODIFIER'])
            if not self.impact.issubset(valid_impacts):
                raise RuntimeError("ERROR: Unrecognised VEP IMPACT provided."+
                                   "Valid values are 'HIGH', 'MODERATE', " +
                                   "'LOW' or 'MODIFIER'")
        if len(biotypes) == 0:
            biotypes = ['default']
        for b in biotypes:
            lb = b.lower()
            if lb == 'all':
                self.biotypes = None
                break
            elif lb == 'default':
                self.biotypes.update(default_biotypes)
            else:
                if lb in valid_biotypes:
                    self.biotypes.add(lb)
                else:
                    raise RuntimeError("ERROR: Unrecognised VEP biotype " +
                                       "'{}'".format(b))
        self.canonical = canonical
        self.loftee = loftee
        self.filter_flagged = filter_flagged_features
        required = ['Consequence', 'BIOTYPE']
        if self.impact:
            required.append('IMPACT')
        if self.canonical:
            required.append('CANONICAL')
        if self.loftee:
            required.append('LoF')
        if self.filter_flagged:
            required.append('FLAGS')
        for rq in required:
            try:
                if rq not in vcf.header.csq_fields:
                    raise RuntimeError("Could not find required VEP " +
                                       "annotation '{}' in VCF".format(rq))
            except HeaderError:
                raise RuntimeError("Could not identify CSQ or ANN fields in " +
                                "VCF header. Please ensure your input is " +
                                "annotated with Ensembl's VEP")
        self.gene_filter = gene_filter
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
        self.retain_labels = defaultdict(set)
        if retain_labels:
            for lbl in retain_labels:
                split_lbl = lbl.split('=', 1)
                if len(split_lbl) != 2:
                    raise RuntimeError("Error in --retain_label - VEP " +
                                       "annotation and value must be " +
                                       "separated by a '=' character")
                self.retain_labels[split_lbl[0]].add(split_lbl[1])
        self.blacklist = self._read_blacklist(blacklist)
        self.pathogenic = pathogenic
        self.no_conflicted = no_conflicted
        self.g2p = g2p
        self.check_g2p_consequence = check_g2p_consequence
        if pathogenic:
            self.path_fields = self._get_path_fields(vcf)

    def filter(self, record):
        filter_alleles = [True] * (len(record.ALLELES) -1)
        #whether an ALT allele should be filtered or not
        filter_af = [False] * (len(record.ALLELES) -1)
        try:
            filter_csq = [True] * len(record.CSQ)
            #whether each csq should be filtered or not
        except HeaderError:
            raise RuntimeError("Could not identify CSQ or ANN fields in VCF " +
                               "header. Please ensure your input is annotated " +
                               "with Ensembl's VEP")
        i = -1
        for c in record.CSQ:
            i += 1
            alt_i = c['alt_index'] -1
            if filter_af[alt_i]: #already filtered on freq for this allele
                continue
            if self.canonical:
                try:
                    if c['CANONICAL'] != 'YES':
                        continue
                except KeyError:
                    pass
            if self.filter_flagged:
                try:
                    if c['FLAGS']:
                        continue
                except KeyError:
                    pass
            if (self.biotypes is not None and
                c['BIOTYPE'].lower() not in self.biotypes):
                continue
            if self.gene_filter:
                if not self.gene_filter.target_in_csq(c):
                    continue
            if self.g2p:
                if c['SYMBOL'] not in self.g2p.g2p:
                    continue
            if self.blacklist and c['Feature'] in self.blacklist:
                continue
            if (self.freq or self.min_freq or self.filter_known or
                    self.filter_novel):
                known = False
                for af in self.freq_fields:
                    if c[af] == '' or c[af] == '.':
                        continue
                    try:
                        c_af = float(c[af])
                    except ValueError:
                        try:
                            c_af = max(float(x) for x in c[af].split('&') if x
                                       != '.')
                        except ValueError:
                            continue
                    known = True
                    if self.filter_known:
                        filter_af[alt_i] = True
                        break
                    if self.freq:
                        if c_af >= self.freq:
                            filter_af[alt_i] = True
                            break
                    if self.min_freq:
                        if c_af < self.min_freq:
                            filter_af[alt_i] = True
                            break
                if self.filter_novel and not known:
                    filter_af[alt_i] = True
                if filter_af[alt_i]:
                    continue
            if (self.csq is None and self.impact is None and
                    not self.check_g2p_consequence):
                #if only using biotypes/MAF for filtering
                filter_alleles[alt_i] = False
                filter_csq[i] = False
                continue
            if self.pathogenic and self._has_pathogenic_annotation(c, record):
                filter_alleles[alt_i] = False
                filter_csq[i] = False
                continue
            if self._retain_label_matched(c):
                filter_alleles[alt_i] = False
                filter_csq[i] = False
                continue
            if self.check_g2p_consequence and self.g2p:
                filt_csq = self.g2p.consequences_from_gene(c['SYMBOL'])
            else:
                filt_csq = self.csq
            for s_csq in [x.lower() for x in c['Consequence'].split('&')]:
                matches_csq = False
                matches_impact = False
                if filt_csq is not None and s_csq in filt_csq:
                    matches_csq = True
                if self.impact is not None and c['IMPACT'] in self.impact:
                    matches_impact = True
                if matches_csq or matches_impact:
                    if self.in_silico and s_csq == 'missense_variant':
                        do_filter = self.in_silico.filter(c)
                        if not do_filter:
                            filter_alleles[alt_i] = False
                            filter_csq[i] = False
                            break
                        filter_csq[i] = do_filter
                    elif self.splice_in_silico and s_csq.startswith("splice"):
                        do_filter = self.splice_in_silico.filter(c)
                        if not do_filter:
                            filter_alleles[alt_i] = False
                            filter_csq[i] = False
                            break
                        filter_csq[i] = do_filter
                    elif self.loftee and (s_csq in lof_csq or matches_impact
                                          and c['IMPACT'] == 'HIGH'):
                        if c['LoF'] == 'HC':
                            filter_alleles[alt_i] = False
                            filter_csq[i] = False
                            break
                    else:
                        filter_alleles[alt_i] = False
                        filter_csq[i] = False
                        break
        return filter_alleles, filter_csq

    def _retain_label_matched(self, csq):
        for k,v in self.retain_labels.items():
            for lbl in csq[k].split('&'):
                if lbl in v:
                    return True
        return False

    def _has_pathogenic_annotation(self, csq, record):
        path = []
        benign = []
        for annot in self.path_fields:
            if not csq[annot]:
                continue
            assertions = csq[annot].split('&')
            if annot == 'clinvar_clnsig':
                #benign = 2, likely benign = 3
                #likely pathognic = 4, pathogenic = 5
                try:
                    benign.extend((4 > int(x) > 1 for x in assertions))
                    path.extend((6 > int(x) > 3 for x in assertions))
                except ValueError:
                    self.logger.warn("Error parsing 'clinvar_clnsig' fields " +
                                     "at {}:{} - expected".format(record.CHROM,
                                                                  record.POS) +
                                     " numeric values.")
            else:
                benign.extend(('benign' in x for x in assertions))
                path.extend(('pathogenic' in x for x in assertions))
        if self.no_conflicted:
            return any(path) and not any(benign)
        return any(path)

    def _read_csq_file(self):
        data_file = os.path.join(os.path.dirname(__file__),
                                 "data",
                                 "vep_classes.tsv")
        return self._get_valid_and_default(data_file)

    def _read_biotype_file(self):
        data_file = os.path.join(os.path.dirname(__file__),
                                 "data",
                                 "biotypes.tsv")
        return self._get_valid_and_default(data_file)

    def _get_valid_and_default(self, data_file):
        defaults = list()
        valid = list()
        with open(data_file,encoding='UTF-8') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                cols = line.rstrip().split('\t')
                if len(cols) < 2:
                    continue
                valid.append(cols[0].lower())
                if cols[1] == 'default':
                    defaults.append(cols[0].lower())
        return defaults, valid

    def _read_maf_file(self):
        data_file = os.path.join(os.path.dirname(__file__),
                                 "data",
                                 "vep_maf.tsv")
        values = []
        with open(data_file,encoding='UTF-8') as fh:
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

    def _read_blacklist(self, blacklist):
        '''
            Reads a file of feature IDs to ignore. Only reads first
            non-whitespace text from each line.
        '''
        if blacklist is None:
            return None
        with open (blacklist, 'rt') as bfile:
            features = set(line.split()[0] for line in bfile)
        self.logger.info("{:,} unique feature IDs blacklisted from {}.".format(
                         len(features), blacklist))
        return features

    def _get_path_fields(self, vcf):
        cln_fields = ['CLIN_SIG', 'clinvar_clnsig']
        path_fields = [f for f in vcf.header.csq_fields if f in cln_fields]
        if not path_fields:
            self.logger.warn("No compatible ClinVar VEP annotations found " +
                             "for use with pathogenic allele identification.")
        return path_fields

    def _get_logger(self, logging_level):
        logger = logging.getLogger(__name__)
        if not logger.hasHandlers():
            logger.setLevel(logging_level)
            formatter = logging.Formatter(
                        '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
            ch = logging.StreamHandler()
            ch.setLevel(logger.level)
            ch.setFormatter(formatter)
            logger.addHandler(ch)
        return logger
