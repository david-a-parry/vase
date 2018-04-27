import operator
from natsort import natsorted
from .genomic_interval import GenomicInterval
from .bed_parser import BedParser

class RegionIter(BedParser):
    '''
        For a list of intervals in format 'chr1:1000-2000' sort and
        merge overlapping intervals. Merged intervals are iterable as
        GenomicInterval objects with unmerged intervals retained in
        the 'regions' property of the GenomicInterval object.
    '''

    __slots__ = ['intervals', 'current']

    def __init__(self, regions):
        '''
            Parses each region provided into a GenomicInterval object,
            merging overlapping intervals.
        '''
        self.intervals = self._parse_regions(regions)
        self.current = 0

    def _parse_regions(self, regions):
        intervals = (self._interval_from_string(x) for x in regions)
        return self._merge_regions(intervals)

    def _interval_from_string(self, interval):
        try:
            chrom, pos = interval.split(":")
            pos_split = pos.split("-")
            if len(pos_split) == 2:
                start, end = pos_split
            else:
                start, end = pos,pos
            start = int(start.replace(',', '')) -1 #0-based, same as BED
            end = int(end.replace(',', ''))
            return [chrom, start, end]
        except ValueError as e:
            raise RuntimeError("Invalid --region specified: '{}'"
                               .format(interval) + "\n" + str(e))

