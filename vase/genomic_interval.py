

class GenomicInterval(object):
    '''
        Simple class for representing a genomic interval, potentially
        merged from several overlapping BED regions. '''

    __slots__ = ['contig', 'start', 'end', 'regions']

    def __init__(self, interval):
        '''
            Args:
                interval:   A list of columns from a BED file.
        '''
        self.contig = interval[0]
        self.start = int(interval[1]) #should be 0-based
        self.end = int(interval[2])
        self.regions = [interval]
        if start >= end:
            raise ValueError("Start of BED interval can not be equal to or " +
                             "greater than end (for interval {}:{}-{})"
                             .format(self.contig, self.start+1, self.end))

    def __str__(self):
        return self.contig + ':' + str(self.start+1) + '-' + str(self.end)

    def __eq__(self, other):
        return (self.contig == other.contig and self.start == other.start and
                self.end == other.end)

    def __ne__(self, other):
        return (self.contig != other.contig or self.start != other.start or
                self.end != other.end)

    def __lt__(self, other):
        return (self.contig < other.contig or
               (self.contig == other.contig and
               (self.start < other.start or self.start == other.start and
                self.end < other.end)))

    def __le__(self, other):
        return (self.contig < other.contig or
               (self.contig == other.contig and self.end <= other.end))

    def __gt__(self, other):
        return (self.contig > other.contig or
               (self.contig == other.contig and
               (self.start > other.start or self.start == other.start and
                self.end > other.end)))

    def __ge__(self, other):
        return (self.contig > other.contig or
               (self.contig == other.contig and self.start >= other.start))

    def overlaps(self, other):
        if self.contig != other.contig:
            return False
        elif self.start < other.end and self.end > other.start:
            #start is 0-based so only overlaps if < end
            return True
        return False

    def merge_interval(self, other):
        '''
            Merge an overlapping interval.

            Args:
                interval:   Another GenomicInterval object

        '''
        if not self.overlaps(other):
            raise NonOverlappingIntervalError("Can not merge non-overlapping" +
                    " intervals '{}' and '{}'".format(self, other))
        if other.start < self.start:
            self.start = other.start
        if other.end > self.end:
            self.end = other.end
        self.regions.extend(other.regions)
        self.regions.sort(key=operator.itemgetter(0, 1, 2, 3))


class NonOverlappingIntervalError(ValueError):
    pass

