from .bed_parser import BedParser, BedFormatError
from collections import defaultdict

ENSG = re.compile(r'''^ENS\w*G\d{11}(\.\d+)?''')
ENST = re.compile(r'''^ENS\w*T\d{11}(\.\d+)?''')
ENSP = re.compile(r'''^ENS\w*P\d{11}(\.\d+)?''')
ENSR = re.compile(r'''^ENS\w*R\d{11}(\.\d+)?''')

class VarByRegion(object):
    '''
        Iterate over variants in VcfReader that overlap regions in
        BED file.

    '''
    __slots__ = ['vcfreader', 'bedparser', 'current_region', 'current_targets',
                 'gene_targets']

    def __init__(self, vcfreader, bed, gene_targets=False):
        '''
            Args:
                vcfreader:
                    VcfReader object from parse_vcf module

                bed:
                    Filename for BED file containing regions to retrieve
                    variants from.

                gene_targets:
                    If True, fourth column of BED input contains gene
                    identifiers (either Ensembl ID or gene symbol) to
                    be accessible during traversal via the
                    'current_targets' property.

        '''
        self.gene_targets = gene_targets
        min_col = 4 if gene_targets else 3
        self.bedparser = BedParser(bed, min_col=min_col)
        self.vcfreader = vcfreader
        self.current_region = None
        self.current_targets = defaultdict(list)  # keys are VEP columns,
                                                  # values are lists of IDs

    def __iter__(self):
        return self

    def __next__(self):
        '''
            For each region in bedparser return each overlapping variant
            in vcfreader.
        '''
        if self.current_region is None:
            self._set_interval()
            self.vcfreader.set_region(self.current_region.contig,
                                      self.current_region.start,
                                      self.current_region.end)
        try:
            return next(self.vcfreader.parser)
        except StopIteration:
            self._set_interval()
            self.vcfreader.set_region(self.current_region.contig,
                                      self.current_region.start,
                                      self.current_region.end)
            return next(self)

    def _set_interval(self):
        '''
            Retrieve next GenomicInterval and if using gene_targets set
            current_targets.
        '''
        self.current_region = next(self.bedparser)
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
        for k,v in self.current_targets:
            if csq[k] in v:
                return True
        return False
