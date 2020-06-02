import pysam
import os
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
        self.set_region = self._set_region
        self.set_region(chrom, start, end)

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
