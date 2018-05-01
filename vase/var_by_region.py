from .bed_parser import BedParser, BedFormatError
from collections import defaultdict
import re

ENSG = re.compile(r'''^ENS\w*G\d{11}(\.\d+)?''')
ENST = re.compile(r'''^ENS\w*T\d{11}(\.\d+)?''')
ENSP = re.compile(r'''^ENS\w*P\d{11}(\.\d+)?''')
ENSR = re.compile(r'''^ENS\w*R\d{11}(\.\d+)?''')

class VarByRegion(object):
    '''
        Iterate over variants in VcfReader that overlap regions in
        BED file.

    '''
    __slots__ = ['vcfreader', 'region_iter', 'current_region', 'current_targets',
                 'gene_targets']

    def __init__(self, vcfreader, bed=None, region_iter=None, 
                 gene_targets=False):
        '''
            Args:
                vcfreader:
                    VcfReader object from parse_vcf module

                bed:
                    Filename for BED file containing regions to retrieve
                    variants from. Either this or region_iter argument
                    is required.

                region_iter:
                    An iterator that gives GenomicInterval objects.
                    Either this or bed argument is required.

                gene_targets:
                    If True, fourth column of BED input contains gene
                    identifiers (either Ensembl ID or gene symbol) to
                    be accessible during traversal via the
                    'current_targets' property.

        '''
        self.gene_targets = gene_targets
        if region_iter:
            if bed:
                raise ValueError("bed and region_iter arguments are mutually" +
                                 "exclusive")
            self.region_iter = region_iter
        elif bed:
            min_col = 4 if gene_targets else 3
            self.region_iter = BedParser(bed, min_col=min_col)
        else:
            raise ValueError("Either bed or region_iter argument is required.")
        self.vcfreader = vcfreader
        self.current_region = None
        self.current_targets = defaultdict(list)  # keys are VEP columns,
                                                  # values are lists of IDs

    def __iter__(self):
        return self

    def __next__(self):
        '''
            For each region in region_iter return each overlapping variant
            in vcfreader.
        '''
        if self.current_region is None:
            self._next_interval()
            self.vcfreader.set_region(self.current_region.contig,
                                      self.current_region.start,
                                      self.current_region.end)
        record = next(self.vcfreader, None)
        while record is None:
            self._next_interval()
            self.vcfreader.set_region(self.current_region.contig,
                                      self.current_region.start,
                                      self.current_region.end)
            record = next(self.vcfreader, None)
        return record

    def _next_interval(self):
        '''
            Retrieve next GenomicInterval and if using gene_targets set
            current_targets.
        '''
        self.current_region = next(self.region_iter)
        if self.gene_targets:
            self._targets_from_region()

    def _targets_from_region(self):
        ''' Retrieve feature names from GenomicInterval.'''
        self.current_targets.clear()
        for reg in self.current_region.regions:
            for x in reg[3].split('/'):
                if ENST.match(x) or ENSR.match(x):
                    i,c = x,'Feature'
                elif ENSG.match(x):
                    i,c= x,'Gene'
                elif ENSP.match(x):
                    i,c = x,'ENSP'
                else:
                    i,c = x,'SYMBOL'
                self.current_targets[c].append(i)

    def target_in_csq(self, csq):
        '''
            Return True if VEP consequence contains a consequence for
            a target in self.current_targets.

            Args:
                csq:    A dict for a single VEP CSQ annotation (i.e. a
                        single item from parse_vcf's VcfRecord.CSQ
                        attribute)

        '''
        for k,v in self.current_targets.items():
            if csq[k] in v:
                return True
        return False
