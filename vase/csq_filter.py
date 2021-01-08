import os
import logging
from collections import defaultdict
from .utils import get_valid_and_default

lof_csq = {'frameshift_variant', 'stop_gained', 'splice_acceptor_variant',
           'splice_donor_variant'}
_csq_class_data = os.path.join(os.path.dirname(__file__),
                               "data",
                               "vep_classes.tsv")
_biotype_class_data = os.path.join(os.path.dirname(__file__),
                                   "data",
                                   "biotypes.tsv")


class CsqFilter(object):
    '''Parent class for functional-annotation filtering classes.'''

    def __init__(self, vcf, csq_attribute='CSQ', csq=[], impact=[],
                 biotypes=[], csq_classes_file=None,
                 biotype_classes_file=None, retain_labels=[],
                 filter_flagged_features=False, gene_filter=None,
                 blacklist=None, g2p=None, check_g2p_consequence=False):
        '''
            Args:
                vcf:    input VcfReader object

                csq_attribute:
                        Attribute of VcfRecord to use for consequence
                        annotations (e.g. CSQ for VEP annotations or ANN for
                        SnpEff annotations). Default='CSQ'.

                csq:    list of consequence types to keep. If 'default'
                        appears anywhere in this list then the default
                        consequence set (as indicated in
                        data/vep_classes.tsv) will be used. Similarly if
                        'all' appears anywhere in this list no filtering
                        on consequence type will occur.

                impact: list of variant impacts to retain.

                biotypes:
                        Filter consequences for features not of the
                        given biotypes. If not provided the default set
                        of biotypes (as indicated in data/biotypes.tsv)
                        will be used for biotype filtering.

                retain_labels:
                        Do not filter on consequence type if the
                        following values are present for a label. Labels
                        and values must be separated by '=' sign. For
                        example, to retain any consequence which has
                        a SnpEff annotation named 'FOO' with  value 'BAR'
                        use 'FOO=BAR'.

                filter_flagged_features:
                        Filter consequences on features which have warnings
                        from SnpEff's 'ERRORS / WARNINGS / INFO' field.

                gene_filter:
                        VarByRegion object from vase.var_by_region. If
                        provided, consequences will be filtered if they
                        do not alter the features specified in the
                        VarByRegion object for the current region.

                blacklist:
                        File containing a list of Feature IDs to ignore.

                g2p:
                        G2P object from vase.g2p for filtering on
                        presence and/or requirements from a G2P file.

                check_g2p_consequence:
                        If a G2P object is provided above, require that
                        that the observed consequence matches the
                        'mutation consequence' in the G2P file.

        '''
        self.csq_attribute = csq_attribute
        self.csq_header_attribute = csq_attribute.lower() + "_fields"
        if not csq_classes_file:
            csq_classes_file = _csq_class_data
        default_csq, valid_csq = get_valid_and_default(csq_classes_file)
        if not biotype_classes_file:
            biotype_classes_file = _biotype_class_data
        default_biotypes, valid_biotypes = get_valid_and_default(
            biotype_classes_file)
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
                    raise RuntimeError("ERROR: Unrecognised consequence " +
                                       "class '{}'".format(c))
        if impact:
            self.impact = set((x.upper() for x in impact))
            valid_impacts = set(['HIGH', 'MODERATE', 'LOW', 'MODIFIER'])
            if not self.impact.issubset(valid_impacts):
                raise RuntimeError("ERROR: Unrecognised IMPACT provided. " +
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
        self.filter_flagged = filter_flagged_features
        self.gene_filter = gene_filter
        self.g2p = g2p
        self.check_g2p_consequence = check_g2p_consequence
        required = self.get_required_header_fields()
        for rq in required:
            try:
                if rq not in getattr(vcf.header, self.csq_header_attribute):
                    raise RuntimeError("Could not find required consequence " +
                                       "annotation '{}' in VCF".format(rq))
            except KeyError:
                raise RuntimeError("Could not identify CSQ or ANN fields in " +
                                   "VCF header. Please ensure your input is " +
                                   "annotated with the relevant tool")
        self.retain_labels = defaultdict(set)
        if retain_labels:
            for lbl in retain_labels:
                split_lbl = lbl.split('=', 1)
                if len(split_lbl) != 2:
                    raise RuntimeError("Error in --retain_label - CSQ field " +
                                       "and value must be separated by a " +
                                       "'=' character")
                self.retain_labels[split_lbl[0]].add(split_lbl[1])
        self.blacklist = self._read_blacklist(blacklist)

    def filter(self, record):
        try:
            csqs = getattr(record, self.csq_attribute)
        except KeyError:
            raise RuntimeError("Could not identify CSQ or ANN fields in VCF " +
                               "header. Please ensure your input is " +
                               "annotated with the appropriate tool.")
        filter_csq = [True] * len(csqs)
        # whether each csq should be filtered or not
        filter_alleles = [True] * len(record.alts)
        # whether an ALT allele should be filtered or not
        override_alleles = [False] * len(record.alts)
        # whether an ALT allele should be filtered or not
        for i, c in enumerate(csqs):
            alt_i = c['alt_index'] - 1
            if override_alleles[alt_i]:  # always filter this allele
                continue
            try:
                f_csq, filter_alt = self.filter_csq(c)
            except Exception as e:
                self.logger.error("CSQ filtering failed at {}:{}".format(
                    record.chrom, record.pos))
                self.logger.error(e)
                raise
            if filter_alt:
                override_alleles[alt_i] = True
                filter_alleles[alt_i] = True
            elif not f_csq:
                filter_csq[i] = False
                filter_alleles[alt_i] = False
        return filter_alleles, filter_csq

    def _retain_label_matched(self, csq):
        for k, v in self.retain_labels.items():
            for lbl in csq[k].split('&'):
                if lbl in v:
                    return True
        return False

    def _read_blacklist(self, blacklist):
        '''
            Reads a file of feature IDs to ignore. Only reads first
            non-whitespace text from each line.
        '''
        if blacklist is None:
            return None
        with open(blacklist, 'rt') as bfile:
            features = set(line.split()[0] for line in bfile)
        self.logger.info("{:,} unique feature IDs blacklisted from {}.".format(
                         len(features), blacklist))
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

    def get_required_header_fields(self):
        '''
        Check which CSQ/ANN annotation fields are required given arguments
        passed to __init__
        '''
        raise NotImplementedError("Child implementation required")

    def filter_csq(csq):
        '''
        Returns two boolean values. The first indicates whether the consequence
        annotation should be filtered. The second indicates whether the ALT
        allele should be filtered irrespective of the given or any other
        consequence annotation.
        '''
        raise NotImplementedError("Child implementation required")
