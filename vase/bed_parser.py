import gzip
from .interval_iter import IntervalIter

class BedParser(IntervalIter):
    '''
        For a given BED file, read into memory and merge overlapping
        intervals. Merged intervals are iterable as GenomicInterval
        objects with unmerged intervals retained in the 'regions'
        property of the GenomicInterval object.
    '''

    __slots__ = ['bed', 'min_col', 'intervals']

    def __init__(self, bed, min_col=3):
        '''
            Opens given bed file, reads into memory. Regions are sorted
            and merged to provide non-overlapping intervals for
            traversal.
        '''
        self.bed = bed
        self.min_col = min_col if min_col > 3 else 3
        intervals = self._read_bed()
        super().__init__(intervals)

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
        return regions


class BedFormatError(ValueError):
    pass
