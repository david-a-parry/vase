import operator
from natsort import natsorted
from .genomic_interval import GenomicInterval

class IntervalIter(object):
    '''
        For a given set of GenomicInterval objects, merge overlapping
        intervals and create an interator over these merged intervals.
    '''

    __slots__ = ['intervals', 'next_index']

    def __init__(self, intervals):
        self.next_index = 0
        self.intervals = self._merge_regions(intervals)

    def __iter__(self):
        return self

    def __next__(self):
        if self.next_index < len(self.intervals):
            self.next_index += 1
            return self.intervals[self.current_index]
        raise StopIteration

    @property
    def previous_interval(self):
        '''
            Get the interval preceding the current interval or None
            if there is no previous interval.
        '''
        if self.current_index > 0:
            return self.intervals[self.current_index-1]
        return None

    @property
    def current_index(self):
        '''
            next_index is always incremented before returning the
            next region to make the iteration work. We return the index
            of the current region here, not the value of next_index.
        '''
        return self.next_index - 1

    @current_index.setter
    def current_index(self, value):
        self.next_index = value + 1

    def _merge_regions(self, regions):
        ''' Return a list of merged regions as GenomicInterval objects.'''
        if self._sort_needed(regions):
            regions = natsorted(regions, key=operator.itemgetter(0, 1, 2))
        genomic_intervals = []
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

    def _sort_needed(self, regions):
        if not isinstance(regions, list):
            return True
        prev = None
        seen_chroms = set()
        for r in regions:
            if prev is not None:
                if r[0] != prev[0]:
                    if r[0] in seen_chroms:
                        return True
                    seen_chroms.add(prev[0])
                elif r[1] < prev[1] or r[2] < prev[2]:
                    return True
            prev = r
        return False
