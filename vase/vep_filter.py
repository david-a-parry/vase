import sys
import logging
from parse_vcf import *
from .insilico_filter import *


class VepFilter(object):
    '''An object that filters VCF records based on annotated VEP data.'''

    def __init__(self, vcf, csq=[], canonical=False, biotypes=[], in_silico=[],
                 filter_unpredicted=False, keep_any_damaging=False,
                 filter_flagged_features=False, freq=None, min_freq=None,
                 afs=[], gene_filter=None, blacklist=None,filter_known=False,
                 filter_novel=False, logging_level=logging.WARNING):
        '''
            Args:
                vcf:    input VcfReader object

                csq:    list of consequence types to keep. If 'default'
                        appears anywhere in this list then the default
                        consequence set (as indicated in
                        data/vep_classes.tsv) will be used. Similarly if
                        'all' appears anywhere in this list no filtering
                        on consequence type will occur.

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

                logging_level:
                        Logging level to use. Default=logging.WARNING.

        '''
        self.logger = self._get_logger(logging_level)
        default_csq, valid_csq = self._read_csq_file()
        default_biotypes, valid_biotypes = self._read_biotype_file()
        self.csq = set()
        self.biotypes = set()
        if len(csq) == 0:
            csq = ['default']
        for c in csq:
            if c.lower() == 'default':
                self.csq.update(default_csq)
            elif c.lower() == 'all':
                self.csq = None
                break
            else:
                if c.lower() in valid_csq:
                    self.csq.add(c.lower())
                else:
                    raise RuntimeError("ERROR: Unrecognised VEP consequence " +
                                       "class '{}'".format(c))
        if len(biotypes) == 0:
            biotypes = ['default']
        for b in biotypes:
            if b.lower() == 'all':
                self.biotypes = None
                break
            elif b.lower() == 'default':
                self.biotypes.update(default_biotypes)
            else:
                if b in valid_biotypes:
                    self.biotypes.add(b)
                else:
                    raise RuntimeError("ERROR: Unrecognised VEP biotype " +
                                       "'{}'".format(b))
        required = ['Consequence', 'BIOTYPE']
        self.canonical = canonical
        if self.canonical:
            required.append('CANONICAL')
        self.filter_flagged = filter_flagged_features
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
        if in_silico:
            in_silico = set(in_silico)
            self.in_silico = InSilicoFilter(in_silico, filter_unpredicted,
                                            keep_any_damaging)
        self.blacklist = self._read_blacklist(blacklist)


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
            if self.biotypes is not None and c['BIOTYPE'] not in self.biotypes:
                continue
            if self.gene_filter:
                if not self.gene_filter.target_in_csq(c):
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
            if self.csq is None: #if only using biotypes/MAF for filtering
                filter_alleles[alt_i] = False
                filter_csq[i] = False
                continue
            consequence = c['Consequence'].split('&')
            for s_csq in consequence:
                if s_csq in self.csq:
                    if self.in_silico and s_csq == 'missense_variant':
                        do_filter = self.in_silico.filter(c)
                        if not do_filter:
                            filter_alleles[alt_i] = False
                        filter_csq[i] = do_filter
                    else:
                        filter_alleles[alt_i] = False
                        filter_csq[i] = False
        return filter_alleles, filter_csq

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
                valid.append(cols[0])
                if cols[1] == 'default':
                    defaults.append(cols[0])
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
        if blacklist is None:
            return None
        features = set()
        with open (blacklist, 'rt') as bfile:
            for line in bfile:
                f = line.rstrip().split()[0]
                if f:
                    features.add(f)
        return features

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
