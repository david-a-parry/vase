import pysam
import os
import gzip
import numpy as np
import struct
from stat import S_ISREG
from .vcf_record import VaseRecord
from .vcf_header import VcfHeader
from .utils import reg2bins

MAX_INT32 = int(2**31 - 1)
# prevent UnicodeDecodeError on latin characters often found in VEP plugin
# annotations
pysam.set_encoding_error_handler('pysam.latin1replace')


class VcfReader(object):
    '''
        A class providing some parsing niceties for working with VCFs.
        Most functionality handled via pysam.VariantFile.
    '''

    def __init__(self, filename, logger=None):
        self.filename = filename
        self.logger = logger
        self.variant_file = pysam.VariantFile(self.filename)
        if filename == '-':
            self._is_reg_file = False
        else:
            self._is_reg_file = S_ISREG(os.stat(self.filename).st_mode)
        if self.filename.endswith(".bcf"):
            self.index = self.filename + '.csi'
        elif self.filename.endswith((".gz", ".bgz")):
            self.index = self.filename + '.tbi'
        else:
            self.index = None
        self.record_iter = (VaseRecord(r, self) for r in self.variant_file)
        self.header = VcfHeader(self)
        self.set_region = self._index_and_set_region
        self.indices = None
        self.walk_chrom = None
        self.prev_walk = (-1, -1)
        self.walk_buffer = []
        self.reseek = False
        self.depth = 5
        self.min_shift = 14
        self.tbi = False

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.record_iter)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.variant_file.close()

    def _index_and_set_region(self, chrom, start=None, end=None, walk=False,
                              walk_region_limit=1000):
        """
            Retrieve records by genomic location rather than reading
            records line-by-line.

            Sets self.reader and self.parser to iterators over the
            records retrieved.

            First call to this method will attempt to create an index if
            it does not already exist.

            Args:
                chrom: chromosome or contig name. Required.

                start: start position on chromosome/contig. 0-based
                       Default = None

                end:   end position on chromosome/contig.
                       Default = None

                walk:  Use "walking" retrieval method which retains file
                       seek position and uses a buffer to retain the
                       last retrieved records. This provides a more
                       efficient method than standard index-based look-ups
                       when performing multiple look-ups in coordinate
                       order. Default = False

                walk_region_limit:
                       If the length of a region is greater than this
                       value and walk is True, the buffer used for
                       walking retrieval will not be used and a reseek will
                       be performed for the next region retrieval.
                       Default = 1000

            >>> v = VcfReader(my_vcf)
            >>> v.set_region('chr1') #get all variants on chr1
            >>> for record in v:
            ... #do something with each record

            Because start coordinates are 0-based, to retrieve a variant
            at (or overlapping) chr1:1000000 use the following:

            >>> v.set_region(chrom='chr1', start=999999 end=1000000)

        """
        self._create_index()
        self.set_region = self._set_region
        self.set_region(chrom, start, end, walk=walk,
                        walk_region_limit=walk_region_limit)

    def _create_index(self):
        if not self._is_reg_file:
            raise TypeError("Cannot run set_region() on a non-regular file")
        elif self.index is None:
            raise TypeError("Cannot run set_region() on non-bcf/vcf.gz file")
        if self.variant_file.index is None:
            if self.logger is not None:
                self.logger.info("No index found for {}".format(self.filename)
                                 + " - creating index")
            preset = 'bcf' if self.variant_file.is_bcf else 'vcf'
            pysam.tabix_index(self.filename, preset=preset)
            self.variant_file = pysam.VariantFile(self.filename)

    def _set_region(self, chrom, start=None, end=None, walk=False,
                    walk_region_limit=1000):
        """
            Retrieve records by genomic location rather than reading
            records line-by-line.

            Sets self.reader and self.parser to iterators over the
            records retrieved.

            Args:
                chrom: chromosome or contig name. Required.

                start: start position on chromosome/contig. 0-based
                       Default = None

                end:   end position on chromosome/contig.
                       Default = None

                walk:  Use "walking" retrieval method which retains file
                       seek position and uses a buffer to retain the
                       last retrieved records. This provides a more
                       efficient method than standard index-based look-ups
                       when performing multiple look-ups in coordinate
                       order. Default = False

                walk_region_limit:
                       If the length of a region is greater than this
                       value and walk is True, the buffer used for
                       walking retrieval will not be used and a reseek will
                       be performed for the next region retrieval.
                       Default = 1000

            >>> v = VcfReader(my_vcf)
            >>> v.set_region('chr1') #get all variants on chr1
            >>> for record in v:
            ... #do something with each record

            Because start coordinates are 0-based, to retrieve a variant
            at (or overlapping) chr1:1000000 use the following:

            >>> v.set_region(chrom='chr1', start=999999 end=1000000)

        """
        if walk:
            self.record_iter = self.walk(chrom, start, end, walk_region_limit)
        else:
            try:
                region_iter = self.variant_file.fetch(chrom, start, end)
                self.record_iter = (VaseRecord(r, self) for r in region_iter)
            except ValueError:
                self.record_iter = iter([])  # ignore missing contigs

    def _read_index(self):
        self._create_index()
        with gzip.open(self.index, 'rb') as f:
            magic = f.read(4)
            if magic == b'TBI\x01':
                self.tbi = True
                return self._read_tbi(f)
            elif magic == b'CSI\x01':
                return self._read_csi(f)
            else:
                raise ValueError('Invalid index - wrong magic number ' +
                                 '({}) for {}'.format(magic, self.index))

    def _read_tbi(self, f):
        '''
            Called for filehandle AFTER reading first 4 bytes
            (magic number)
        '''
        header = np.frombuffer(f.read(4 * 8), dtype=np.int32)
        names = f.read(header[7]).split(b'\x00')
        ridx = dict()
        for i in range(len(names)):
            if names[i] == b'':
                continue
            bindx = dict()
            for j in range(struct.unpack('<i', f.read(4))[0]):  # n_bins
                bin_key = struct.unpack('<I', f.read(4))[0]  # bin
                n_chunk = struct.unpack('<i', f.read(4))[0]  # n_chunk
                bindx[bin_key] = np.frombuffer(f.read(8 * 2 * n_chunk),
                                               dtype=np.uint64).reshape(
                                                   n_chunk, -1)
            d = {'bindx': bindx}
            n_intv = struct.unpack('<i', f.read(4))[0]
            d['n_intv'] = n_intv
            d['ioff'] = np.frombuffer(f.read(8 * n_intv), dtype=np.uint64)
            ridx[names[i].decode()] = d
        return ridx

    def _read_csi(self, f):
        '''
            Called for filehandle AFTER reading first 4 bytes
            (magic number)
        '''
        csindex = dict()
        self.min_shift, self.depth, l_aux  = np.frombuffer(f.read(4 * 3),
                                                           dtype=np.int32)
        aux = np.frombuffer(f.read(4 * l_aux), dtype=np.uint8)
        n_ref = np.frombuffer(f.read(4), dtype=np.int32)[0]
        for i in range(n_ref):
            chrom = self.variant_file.get_reference_name(i)
            bindx = dict()
            offsets = list()
            n_bins = struct.unpack('<i', f.read(4))[0]
            for j in range(n_bins):  # n_bins
                bin_key = struct.unpack('<i', f.read(4))[0]
                offs = np.frombuffer(f.read(8), dtype=np.int64)
                n_chunk = struct.unpack('<i', f.read(4))[0]  # n_chunk
                bindx[bin_key] = np.frombuffer(f.read(8 * 2 * n_chunk),
                                               dtype=np.uint64).reshape(
                                                  n_chunk, -1)  # chunk_beg/end
            #offsets = np.ravel([x for y in bindx.values() for x in y])
            #d = {'bindx': bindx, 'ioff': offsets}
            d = {'bindx': bindx}
            csindex[chrom] = d
        return csindex

    def walk(self, chrom, start=None, end=None, region_limit=1000):
        '''
            Retrieve records given by chromosome, start and end
            coordinates using a method which retains file seek position
            and uses a buffer to retain the last retrieved records. This
            provides a more efficient method than the 'set_region' method
            when performing multiple look-ups in coordinate order where a
            the same compressed VCF chunk may be accessed multiple times
            by consecutive look-ups (e.g. if looking up variants in this
            VCF that overlap variants in another VCF).

            Args:
                chrom: chromosome or contig name. Required.

                start: start position on chromosome/contig. 0-based
                       Default = None

                end:   end position on chromosome/contig.
                       Default = None

                region_limit:
                       If the length of a region is greater than this
                       value, buffer will not be used and a reseek will
                       be performed for the next region retrieval.
                       Default = 1000

        '''
        start = 0 if start is None else start
        end = MAX_INT32 if end is None else end
        if self.indices is None:
            self.indices = self._read_index()
        if self.walk_chrom != chrom or start < self.prev_walk[0]:
            self.walk_chrom = chrom
            self.reseek = True
        if chrom not in self.indices:
            return []
        self.prev_walk = (start, end)
        use_buffer = end - start < region_limit
        if self.tbi:
            i = start >> 14
            if i >= self.indices[chrom]['n_intv']:
                return []
            min_ioff = self.indices[chrom]['ioff'][i]
        else:
            min_ioff = 0
        # binning index: includes records in large interval
        bins = [self.indices[chrom]['bindx'][k] for k in
                reg2bins(start, end, self.min_shift, self.depth)
                if k in self.indices[chrom]['bindx']]
        if not bins:
            return []
        overlap = np.concatenate(bins)
        # coupled binning and linear indices (if tbi), remove low level bins
        cbins = np.sort(np.ravel(overlap[overlap[:, 1] >= min_ioff]))
        if cbins.size == 0:
            return []
        chunk_begin, chunk_end = cbins[0], cbins[-1]
        if self.reseek or chunk_begin > self.variant_file.tell():
            self.variant_file.seek(chunk_begin)
            self.walk_buffer = []
        elif self.walk_buffer:
            remove_i = set()
            if start < self.walk_buffer[-1].stop:
                for i in range(len(self.walk_buffer)):
                    if self.walk_buffer[i].start >= end:
                        break
                    if self.walk_buffer[i].stop > start:
                        yield VaseRecord(self.walk_buffer[i], self)
                    else:
                        remove_i.add(i)
                if remove_i:
                    self.walk_buffer = [self.walk_buffer[i] for i in
                                        range(len(self.walk_buffer)) if i not
                                        in remove_i]
            else:
                self.walk_buffer = []
        if not self.walk_buffer or self.walk_buffer[-1].start < end:
            for record in self.variant_file:
                if record.start >= end or self.variant_file.tell() > chunk_end:
                    if use_buffer:
                        if record.chrom == chrom:
                            self.walk_buffer.append(record)
                    break
                if record.stop > start:
                    yield VaseRecord(record, self)
                    if use_buffer:
                        self.walk_buffer.append(record)
        self.reseek = not use_buffer
