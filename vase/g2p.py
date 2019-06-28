import sys
from collections import defaultdict
from .utils import csv_to_dict

allelic_req_to_label = {'biallelic'                 : ['recessive'],
                        'digenic'                   : [],
                        'hemizygous'                : ['recessive'],
                        'imprinted'                 : [],
                        'mitochondrial'             : [],
                        'monoallelic'               : ['de novo', 'dominant'],
                        'mosaic'                    : ['de novo', 'dominant'],
                        'x-linked dominant'         : ['de novo', 'dominant'],
                        'x-linked over-dominance'   : ['de novo', 'dominant'],
                        ''                          : [],
                       }

mutation_to_csq = {
    'loss of function'                   : ['frameshift_variant',
                                           'stop_gained',
                                           'start_lost',
                                           'initiator_codon_variant',
                                           'splice_acceptor_variant',
                                           'splice_donor_variant',
                                           'transcript_ablation',
                                           'feature_truncation',],
    'uncertain'                          : None,
    ''                                   : None,
    'all missense/in frame'              : ['missense_variant',
                                           'inframe_deletion',
                                           'inframe_insertion',],
    'dominant negative'                  : ['missense_variant',
                                           'inframe_deletion',
                                           'inframe_insertion',],
    'activating'                         : ['missense_variant',
                                           'inframe_deletion',
                                           'inframe_insertion',],
    'cis-regulatory or promotor mutation': ['regulatory_region_ablation',
                                           'regulatory_region_amplification',
                                           'regulatory_region_variant',
                                           'TFBS_ablation',
                                           'TFBS_amplification',
                                           'TF_binding_site_variant',],
    'part of contiguous gene duplication': [],
    'increased gene dosage'              : [],
    'gain of function'                   : ['missense_variant',
                                           'inframe_deletion',
                                           'inframe_insertion',],
    '5_prime or 3_prime UTR mutation'    : ['3_prime_UTR_variant',
                                           '5_prime_UTR_variant',],
}


class G2P(object):
    ''' Filter variants based on requirements from a G2P CSV file.'''

    def __init__(self, g2p_file):
        self.g2p = self._read_g2p_csv(g2p_file)

    def _read_g2p_csv(self, g2p):
        required_fields = ['gene symbol', 'disease name', 'DDD category',
                           'allelic requirement', 'mutation consequence',
                           'organ specificity list', 'prev symbols']
        g2p = csv_to_dict(g2p, 'gene symbol', required_fields)
        prev_symbol_d = defaultdict(list)
        for row in [x for y in g2p.values() for x in y]:
            if row['prev symbols']:
                for prev in row['prev symbols'].split(';'):
                    prev_symbol_d[prev].append(row)
        g2p.update(prev_symbol_d)
        return g2p

    def consequence_requirement_met(self, record):
        '''
            For a given VcfRecord, for each CSQ annotation from Ensembl
            VEP return True or False indicating whether the consequence
            annotation matches a 'mutation consequence' requirement from
            the gene's G2P annotations.
        '''
        try:
            return [csq_matches_requirement(x) for x in record.CSQ]
        except HeaderError:
            raise RuntimeError("Could not identify CSQ or ANN fields in VCF " +
                               "header. Please ensure your input is " +
                               "annotated with Ensembl's VEP")

    def csq_matches_requirement(self, csq, keep_uncertain=True):
        '''
            For a single CSQ annotation (from a VcfRecord) return True
            if the consequence annotation matches a 'mutation
            consequence' requirement from the gene's G2P annotations.

            Args:
                csq:    A single item from a VcfRecord's CSQ property.

                keep_uncertain:
                        If True return True for any consequence for G2P
                        genes with a value of 'uncertain' or '' in the
                        'mutation consequence' column. If False return
                        False for these consequences.
        '''
        if csq['SYMBOL'] in self.g2p:
            for req in (x['mutation consequence'] for x in
                        self.g2p[csq['SYMBOL']]):
                if mutation_to_csq[req] is None:
                    if keep_uncertain:
                        return True
                    else:
                        return False
                elif any(x in csq['Consequence'].split('&') for x in
                         mutation_to_csq[req]):
                    return True
        return False

    def consequences_from_gene(self, gene, uncertain_to_none=True):
        '''
            Args:
                gene:   gene to search for

                uncertain_to_none:
                        set to True if you want to return None if any
                        of the 'mutation consequence' values for the gene
                        is 'uncertain'. If set to False, no consequences
                        will be returned for annotations with an
                        'uncertain' 'mutation consequence'.
        '''
        csqs = set()
        for d in self.g2p[gene]:
            if (d['mutation consequence'] == 'uncertain' or
                    d['mutation consequence'] == ''):
                if uncertain_to_none:
                    return None
            else:
                csqs.update(mutation_to_csq[d['mutation consequence']])
        return csqs

    def allelic_requirement_met(self, record, inheritance):
        '''
            Return True or False for each CSQ annotation in a record
            indicating whether the consequence affects a gene in the G2P
            data associated with the given inheritance pattern.

            Args:
                record: VcfRecord object with VEP annotations.

                inheritance:
                        Inheritance patten to check - i.e. one of
                        'recessive', 'dominant' or 'de novo'.

        '''
        return (any(inheritance in allelic_req_to_label[y] for x in
                    self.g2p[csq['SYMBOL']] for y in
                    self.g2p[csq['SYMBOL']] for y in
                    x['allelic requirement'].split(',')) for csq in record.CSQ)

    def csq_and_allelic_requirement_met(self, record, inheritance,
                                        keep_uncertain=True):
        '''
            Return True or False for each CSQ annotation in a record
            indicating whether the consequence affects a gene in the G2P
            data associated with the given inheritance pattern with a
            matching consequence type.

            Args:
                record: VcfRecord object with VEP annotations.

                inheritance:
                        Inheritance patten to check - i.e. one of
                        'recessive', 'dominant' or 'de novo'.

                keep_uncertain:
                        If True return True for any consequence for G2P
                        genes with a value of 'uncertain' or '' in the
                        'mutation consequence' column. If False return
                        False for these consequences.
        '''
        met = []
        for csq in record.CSQ:
            verdict = False
            for d in self.g2p[csq['SYMBOL']]:
                req = d['allelic requirement']
                if any(inheritance in allelic_req_to_label[r] for r in
                       req.split(',')):
                    if mutation_to_csq[d['mutation consequence']] is None:
                        if keep_uncertain:
                            verdict = True
                            break
                    elif any(x in csq['Consequence'].split('&') for x in
                           mutation_to_csq[d['mutation consequence']]):
                        verdict = True
                        break
            met.append(verdict)
        return met
