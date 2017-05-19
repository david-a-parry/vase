import sys
from .parse_vcf.parse_vcf import * 
from .insilico_filter import *


class VepFilter(object):
    '''An object that filters VCF records based on annotated VEP data.'''

    def __init__(self, csq=[], canonical=False, biotypes=[], in_silico=[],
                 filter_unpredicted=False, keep_any_damaging=False):

        default_csq, valid_csq = self._read_csq_file()
        default_biotypes, valid_biotypes = self._read_biotype_file()
        self.csq = set()
        self.biotypes = set()
        if len(csq) == 0:
            csq = ['default']
        for c in csq:
            if c.lower() == 'default':
                self.csq.update(default_csq)
            else:
                if c.lower() in valid_csq:
                    self.csq.add(c.lower())
                else:
                    raise Exception("ERROR: Unrecognised VEP consequence " +  
                                    "class '{}'".format(c))
        if len(biotypes) == 0:
            biotypes = ['default']
        for b in biotypes:
            if b.lower() == 'all':
                self.biotypes.update(valid_biotypes)
            elif b.lower() == 'default':
                self.biotypes.update(default_biotypes)
            else:
                if b.lower() in valid_biotypes:
                    self.biotypes.add(b.lower())
                else:
                    raise Exception("ERROR: Unrecognised VEP biotype " +  
                                    "'{}'".format(b))

        self.canonical = canonical
        self.in_silico = False
        if in_silico:
            in_silico = set(in_silico)
            self.in_silico = InSilicoFilter(in_silico, filter_unpredicted, 
                                            keep_any_damaging)

    def filter(self, record):
        filter_alleles = [True] * (len(record.ALLELES) -1)
        #whether an ALT allele should be filtered or not
        try:
            filter_csq = [True] * len(record.CSQ) 
        except HeaderError:
            raise Exception("Could not identify CSQ or ANN fields in VCF " +
                            "header. Please ensure your input is annotated " +
                            "with Ensembl's VEP")
        #whether each csq should be filtered or not
        i = -1
        for c in record.CSQ:
            i += 1
            if self.canonical:
                try: 
                    if c['CANONICAL'] != 'YES':
                        continue
                except KeyError:
                    pass
            if c['BIOTYPE'] not in self.biotypes:
                continue
            consequence = c['Consequence'].split('&')
            for s_csq in consequence:
                if s_csq in self.csq:
                    if self.in_silico and s_csq == 'missense_variant':
                        do_filter = self.in_silico.filter(c)
                        if not do_filter:
                            filter_alleles[c['alt_index'] -1] = False
                        filter_csq[i] = do_filter
                    else:
                        filter_alleles[c['alt_index'] -1] = False
                        filter_csq[i] = False
        return filter_alleles, filter_csq

    def _read_csq_file(self):
        data_file = os.path.join(os.path.dirname(__file__), "data", "vep_classes.tsv")
        return self._get_valid_and_default(data_file)

    def _read_biotype_file(self):
        data_file = os.path.join(os.path.dirname(__file__), "data", "biotypes.tsv")
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
