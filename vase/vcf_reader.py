import pysam
import os
import gzip
import numpy as np
import struct
from stat import S_ISREG
from vase.vcf_record import VaseRecord
from vase.vcf_header import VcfHeader


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
        self.region_limit = 1000

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.record_iter)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.variant_file.close()

    def _index_and_set_region(self, chrom, start=None, end=None):
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
        self.set_region(chrom, start, end)

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

    def _set_region(self, chrom, start=None, end=None):
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

            >>> v = VcfReader(my_vcf)
            >>> v.set_region('chr1') #get all variants on chr1
            >>> for record in v:
            ... #do something with each record

            Because start coordinates are 0-based, to retrieve a variant
            at (or overlapping) chr1:1000000 use the following:

            >>> v.set_region(chrom='chr1', start=999999 end=1000000)

        """
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
            bindx = dict()
            for j in range(struct.unpack('<i', f.read(4))[0]):  # n_bins
                bin_key = struct.unpack('<I', f.read(4))[0]  # bin
                n_chunk = struct.unpack('<i', f.read(4))[0]  # n_chunk
                bindx[bin_key] = np.frombuffer(f.read(8 * 2 * n_chunk),
                                               dtype=np.uint64).reshape(
                                                   n_chunk, -1)
            d = {'bindx': bindx}
            n_intv = struct.unpack('<i', f.read(4))[0]
            d['ioff'] = np.frombuffer(f.read(8 * n_intv), dtype=np.uint64)
            ridx[names[i].decode()] = d
        return ridx

    def _read_csi(self, f):
        '''
            Called for filehandle AFTER reading first 4 bytes
            (magic number)
        '''
        raise NotImplementedError("CSI index retrieval not implemented yet")

    def reg2bins(self, begin, end, n_lvls=5, min_shift=14):
        t, s = 0, min_shift + (n_lvls << 1) + n_lvls
        for l in range(n_lvls + 1):
            b, e = t + (begin >> s), t + (end >> s)
            n = e - b + 1
            for k in range(b, e + 1):
                yield k
                n += 1
            t += 1 << ((l << 1) + l)
            s -= 3

    def walk(self, chrom, start, end):
        recs = []
        if self.indices is None:
            self.indices = self._read_index()
        if self.walk_chrom != chrom:
            self.walk_chrom = chrom
            self.reseek = True
        elif start < self.prev_walk[0]:
            raise RuntimeError("Walk must be done in coordinate order")
        use_buffer = 1 + end - start < self.region_limit
        min_ioff = self.indices[chrom]['ioff'][start >> 14]
        # binning index: record cluster in large interval
        overlap = np.concatenate([self.indices[chrom]['bindx'][k]
                                  for k in self.reg2bins(start, end)
                                  if k in self.indices[chrom]['bindx']])
        # coupled binning and linear indices, filter out low level bins
        chunk_begin, *_, chunk_end = np.sort(
            np.ravel(overlap[overlap[:, 0] >= min_ioff]))
        if self.reseek or chunk_begin > self.variant_file.tell():
            self.variant_file.seek(chunk_begin)
        elif self.walk_buffer and start < self.walk_buffer[-1].stop:
            for record in self.walk_buffer:
                if record.start >= end:
                    break
                if record.stop >= start:
                    recs.append(record)
        if not self.walk_buffer or self.walk_buffer[-1].start < end:
            self.walk_buffer = []
            for record in self.variant_file:
                if record.start >= end:
                    if use_buffer:
                        self.walk_buffer.append(record)
                    break
                if record.stop >= start:
                    recs.append(record)
                    if use_buffer:
                        self.walk_buffer.append(record)
        self.reseek = not use_buffer
        self.prev_walk = (start, end)
        return [VaseRecord(x, self) for x in recs]
