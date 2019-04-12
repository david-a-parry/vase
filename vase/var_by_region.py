from .bed_parser import BedParser, BedFormatError
from collections import defaultdict
import re

ENSG = re.compile(r'''^ENS\w*G\d{11}(\.\d+)?''')
ENST = re.compile(r'''^ENS\w*T\d{11}(\.\d+)?''')
ENSP = re.compile(r'''^ENS\w*P\d{11}(\.\d+)?''')
ENSR = re.compile(r'''^ENS\w*R\d{11}(\.\d+)?''')


class RegionFinder(object):
    '''
        From an IntervalIter object create an index of regions per window
        and provide methods to retrieve regions from contig, start and
        end coordinates.
    '''

    __slots__ = ['regions', 'window_size']

    def __init__(self, interval_iter, window_size=100000):
        self.regions = defaultdict(dict)
        self.window_size = window_size
        for gi in interval_iter:#these should already be coordinate sorted
            r_start = int(gi.start/window_size) * window_size
            r_end = int(gi.end/window_size) * window_size
            for i in range(r_start, r_end + window_size, window_size):
                if i not in self.regions[gi.contig]:
                    self.regions[gi.contig][i] = list()
                self.regions[gi.contig][i].append(gi)

    def fetch(self, contig, start, end):
        if contig not in self.regions:
            return []
        idx_start = int(start/self.window_size) * self.window_size
        idx_end = int(end/self.window_size) * self.window_size
        candidates = []
        for i in range(idx_start, idx_end + 1, self.window_size):
            if i in self.regions[contig]:
                candidates.extend(self.regions[contig][i])
        return self._binsearch_regions(candidates, start, end)

    def _binsearch_regions(self, regions, start, end):
        '''
            Assumes all regions are on the same chromosome. Return all
            overlapping regions.
        '''
        l = 0
        u = len(regions) - 1
        i = self._binsearch(regions, l, u, start, end)
        hits = []
        if i > -1:
            for j in range(i-1, -1, -1):
                if start <= regions[j].end and end > regions[j].start:
                    hits.append(regions[j])
                elif regions[j].end < start:
                    break
            hits.reverse()
            for j in range(i, len(regions)):
                if start <= regions[j].end and end > regions[j].start:
                    hits.append(regions[j])
                elif regions[j].start > end:
                    break
        return hits

    def _binsearch(self, regions, l, u, start, end):
        ''' Find any region overlapping start and end coordinates.'''
        if u < l:
            return -1
        i = int(l + u)//2
        if regions[i].end < start:
            return self._binsearch(regions, i + 1, u, start, end)
        elif regions[i].start > end:
            return self._binsearch(regions, l, i - 1, start, end)
        else:
            return i

class VarByRegion(object):
    '''
        Iterate over variants in VcfReader that overlap regions in
        BED file.

    '''
    __slots__ = ['vcfreader', 'region_iter', 'current_region', 'exclude',
                 'current_targets', 'gene_targets', 'region_finder']

    def __init__(self, vcfreader, bed=None, region_iter=None,
                 gene_targets=False, stream=False, exclude=False):
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

                stream:
                    If True, rather than index-jumping to region
                    locations in VCF, all variants in the VCF will be
                    read and checked for overalap against regions. This
                    allows processing of non-indexed VCFs and also
                    provides a speedup for VCFs with large structural
                    variants.

                exclude:
                    If True, variants that DO NOT overlap with regions
                    will be returned instead. This forces streaming of
                    variants.

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
        self.exclude = exclude
        self.current_region = None
        self.current_targets = defaultdict(list)  # keys are VEP columns,
                                                  # values are lists of IDs
        self.region_finder = None
        if self.exclude:
            stream = True
        if stream:
            self.region_finder = RegionFinder(self.region_iter)
            self.region_iter = None

    def __iter__(self):
        return self

    def __next__(self):
        '''
            For each region in region_iter return each overlapping variant
            in vcfreader.
        '''
        if self.region_finder is not None:
            return self._next_from_region_finder()
        else:
            return self._next_from_region_iterator()

    def _next_from_region_iterator(self):
        if self.current_region is None:
            self._next_interval()
            self.vcfreader.set_region(self.current_region.contig,
                                      self.current_region.start,
                                      self.current_region.end)
        record = self._get_record_if_no_overlap()
        while record is None:
            self._next_interval()
            self.vcfreader.set_region(self.current_region.contig,
                                      self.current_region.start,
                                      self.current_region.end)
            record = self._get_record_if_no_overlap()
        return record

    def _next_from_region_finder(self):
        for record in self.vcfreader:
            regions = self.region_finder.fetch(record.CHROM, record.POS,
                                               record.SPAN)
            if not regions and not self.exclude:
                continue
            elif regions and self.exclude:
                continue
            if regions: #i.e. not using exclude option
                self.current_region = regions[0]
                if self.gene_targets:
                    self.current_targets.clear()
                    for reg in [x for r in regions for x in r.regions]:
                        self._append_targets_from_region(reg)
            return record
        raise StopIteration

    def _get_record_if_no_overlap(self):
        '''
            Ensure we don't return the same record twice by getting the
            next record if the current record overlaps the previous
            interval.
        '''
        record = next(self.vcfreader, None)
        if record is None:
            return None
        if self.region_iter.previous_interval is not None:
            while (self._record_overlaps(record,
                                         self.region_iter.previous_interval)):
                record = next(self.vcfreader, None)
                if record is None:
                    break
        return record

    def _record_overlaps(self, record, interval):
        ''' Return True if record overlaps interval.'''
        return (record.CHROM == interval.contig and record.POS <= interval.end
                and record.SPAN > interval.start)

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
            self._append_targets_from_region(reg)

    def _append_targets_from_region(self, region):
            for x in region[3].split('/'):
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
