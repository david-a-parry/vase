import re
import gzip
from natsort import natsorted
from .genomic_interval import GenomicInterval

class BedParser(object):
    ''' 
        For a given BED file, read into memory and merge overlapping 
        intervals. Merged intervals are iterable as GenomicInterval
        objects with unmerged intervals retained in the 'regions' 
        property of the GenomicInterval object.
    '''

    __slots__ = ['bed', 'min_col', 'intervals', 'current'] 

    def __init__(self, bed, min_col=3):
        ''' 
            Opens given bed file, reads into memory. Regions are sorted 
            and merged to provide non-overlapping intervals for 
            traversal.
        '''
        self.bed = bed
        self.min_col = min_col if min_col > 3 else 3
        self.intervals = self._read_bed()
        self.current = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.current < len(self.intervals):
            self.current += 1
            return self.intervals[self.current-1]
        raise StopIteration
       
    def _read_bed(self):
        regions = []
        if self.bed.endswith((".gz", ".bgz")):
            bfile = gzip.open(self.bed, errors='replace', mode='rt')
        else:
            bfile = open (self.bed, 'rt') 
        for line in bfile:
            if line[0] == '#': continue
            s = line.rstrip().split("\t")
            if len(s) < self.min_col:
                raise BedFormatError("Not enough fields in BED line: " + line)
            try:
                s[1] = int(s[1])
                s[2] = int(s[2])
            except ValueError:
                raise BedFormatError("Columns 2 and 3 must be integers (for " +
                                     "line: "+ line + ")")
            regions.append(s)
        bfile.close()
        return self._merge_regions(regions)

    def _merge_regions(self, regions):
        ''' Return a list of merged regions as GenomicInterval objects.'''
        regions = natsorted(regions, key=operator.itemgetter(0, 1, 2, 3))
        genomic_invtervals = []
        prev_i = None
        for r in regions:
            gi = GenomicInterval(r)
            if prev_i is None:
                prev_i = gi
            elif prev_i.overlaps(gi):
                prev_i.merge_interval(gi)
            else:
                genomic_intervals.append(prev_i)
                prev_i = gi
        if prev_i is not None:
            genomic_intervals.append(prev_i)
        return genomic_intervals


class BedFormatError(ValueError):
    pass
