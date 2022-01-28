import os
import logging
import gzip
import pysam
import numpy as np
from collections import defaultdict, namedtuple
from .utils import read_tbi, reg2bins

CaddRecord = namedtuple('CaddRecord', 'chrom pos stop ref alt raw phred')


class CaddFilter(object):
    '''
        An object that filters/annotates VCF records using CADD PHRED
        scores provided by at least one tabix indexed file of CADD
        scores.
    '''

    def __init__(self, cadd_files=[], cadd_dir=[], min_phred=None,
                 min_raw_score=None, to_score=None,
                 logging_level=logging.WARNING, no_walk=False,
                 force_walk=False, skip_svs=False):
        '''
            Either a directory containing at least one tabix indexed
            file of CADD scores or a list of such files must be
            provided. Optionally a minimum PHRED score for filtering
            records can be provided.

            Args:
                cadd_files:
                    One or more reference CADD files with CADD_raw_score
                    and CADD_PHRED_score columns for variants.

                cadd_dir:
                    One or more directories containing CADD files.
                    Files with '.gz' or '.bgz' extensions will be
                    assumed to be CADD reference files.

                min_phred:
                    Minimum CADD PHRED score for filtering variants.

                min_raw_score:
                    Minimum CADD raw score for filtering variants.

                to_score:
                    Name of file for writing variants which cannot be
                    found in the given CADD reference files. If a
                    variant passed to the annotate_or_filter method is
                    not found it will be written to this in a format
                    suitable for uploading to the CADD server for
                    scoring.

                no_walk:
                    If True, do not use walking retrieval method to find
                    matching records. Walking retrieval reduces
                    unnecessary seeks if look-ups are performed in
                    coordinate order but may prove inefficient if look-up
                    order is random.

                force_walk:
                     By default, if look-ups are performed out of order,
                     record retrieval will fall-back to a standard
                     tabix-style 'fetch' method. If this option is set to
                     True, walking style look-ups will continue to be
                     attempted even if out of order look-ups are detected.
                     Ignored if no_walk is True.

                skip_svs:
                     If True skip comparisons for structural variants
                     passed to annotate_and_filter_record. Default=True.

        '''
        self.to_score_file = None
        self.logger = self._get_logger(logging_level)
        self.walk = not no_walk
        self.skip_svs = skip_svs
        self.force_walk = force_walk
        if cadd_dir:
            cadd_files.extend([os.path.join(cadd_dir, f) for f in
                               os.listdir(cadd_dir) if
                               f.endswith(('.gz', '.bgz')) and
                               os.path.isfile(os.path.join(cadd_dir, f))])
        if not cadd_files:
            if cadd_dir:
                raise RuntimeError("No .gz or .bgz files identified in " +
                                   cadd_dir)
            else:
                raise RuntimeError("No CADD files or directory provided.")
        self.cadd_tabix = list()
        self.indices = dict() if self.walk else None
        self.bgzfs = dict() if self.walk else None
        self._get_tabix_files(cadd_files)
        self._has_chr = self._check_contigs()
        self.phred = min_phred
        self.raw = min_raw_score
        self.info_fields = {
            'CADD_PHRED_score': {'Number': 'A',
                                 'Type': 'Float',
                                 'Description': 'CADD PHRED score added from' +
                                                ' reference files by VASE'},
            'CADD_raw_score': {'Number': 'A',
                               'Type': 'Float',
                               'Description': 'CADD RawScore added from ' +
                                              'reference files by VASE'}
        }
        if to_score is not None:
            if not to_score.endswith('.gz'):
                to_score += '.gz'
            self.to_score_file = gzip.open(to_score, 'wt')
        self.walk_chrom = None
        self.prev_walk = (-1, -1)
        self.walk_buffer = defaultdict(list)
        self.reseek = False
        if self.walk:
            for fh in self.cadd_tabix:
                fh.close()

    def __del__(self):
        if self.to_score_file is not None:
            self.to_score_file.close()

    def annotate_or_filter(self, record):
        '''
            Annotates record with CADD raw and PHRED scores and returns
            a list of booleans indicating whether each allele should be
            filtered (i.e. each allele has a CADD raw or PHRED score
            below threshold).
        '''
        if self.skip_svs and record.IS_SV:
            return [False] * len(record.record.alts)
        filter_alleles = []
        scores = self.score_record(record)
        info_to_add = defaultdict(list)
        i = 0
        for s in scores:
            info_to_add['CADD_raw_score'].append(s[0])
            info_to_add['CADD_PHRED_score'].append(s[1])
            do_filter = False
            if self.raw and s[0] is not None:
                if s[0] < self.raw:
                    do_filter = True
            if self.phred and s[1] is not None:
                if s[1] < self.phred:
                    do_filter = True
            filter_alleles.append(do_filter)
            if self.to_score_file and s[0] is None:
                self._write_for_scoring(record, i)
            i += 1
        record.add_info_fields(info_to_add)
        return filter_alleles

    def score_record(self, record):
        '''
            Returns tuple of raw score and phred score for each allele.
            Returns the scores for the first matching record encountered
            in cadd files.
        '''
        hits = self.search_coordinates(record.record.chrom,
                                       record.record.start,
                                       record.record.stop)
        scores = []
        for i in range(len(record.DECOMPOSED_ALLELES)):
            s = (None, None)
            for h in (x for x in hits if x is not None):
                if self._compare_allele(record.DECOMPOSED_ALLELES[i], h):
                    s = (float(h.raw), float(h.phred))
                    break  # bail on first matching variant
            scores.append(s)
        return scores

    def _compare_allele(self, alt_allele, cadd):
        return (alt_allele.POS == cadd.pos and
                alt_allele.REF == cadd.ref and
                alt_allele.ALT == cadd.alt)

    def _simplify_cadd_record(self, cadd):
        '''
            Return position, ref allele, alt allele, raw and Phred score
            after reducing alleles to their most simple representation.
        '''
        if cadd.startswith('#'):
            return None
        cols = cadd.split("\t")
        if len(cols) < 6:
            self.logger.warn("Not enought columns for CADD record: {}"
                             .format(cadd))
            return None
        ref = cols[2]
        alt = cols[3]
        pos = int(cols[1])
        while len(ref) > 1 and len(alt) > 1:
            if ref[-1] == alt[-1]:               # remove identical suffixes
                ref = ref[:-1]
                alt = alt[:-1]
            else:
                break
        while len(ref) > 1 and len(alt) > 1:
            if ref[0] == alt[0]:                 # remove identical prefixes
                ref = ref[1:]
                alt = alt[1:]
                pos += 1
            else:
                break
        return CaddRecord(cols[0], pos, pos + len(ref) - 1, ref, alt, cols[4],
                          cols[5])

    def walk_coordinates(self, tbx, chrom, start, end, region_limit=1000):
        '''
            See vase.vcf_reader.VcfReader.walk for explanation of this
            retrieval method.
        '''
        recs = []
        idx = self.indices[tbx]
        use_buffer = 1 + end - start < region_limit
        if self.walk_chrom != chrom:
            self.walk_chrom = chrom
            self.reseek = True
        elif start < self.prev_walk[0]:
            self.reseek = True
            if not self.force_walk:
                self.logger.warn("Input is not sorted by coordinate, will  " +
                                 "fall back to slower indvidual index-based " +
                                 "look-ups.")
                self.walk = False
                use_buffer = False
        if chrom not in idx:
            return []
        self.prev_walk = (start, end)
        si = start >> 14
        if si >= idx[chrom]['n_intv']:
            return []
        min_ioff = idx[chrom]['ioff'][si]
        bins = [idx[chrom]['bindx'][k] for k in reg2bins(start, end)
                if k in idx[chrom]['bindx']]
        if not bins:
            return []
        overlap = np.concatenate(bins)
        # coupled binning and linear indices, filter out low level bins
        cbins = np.sort(np.ravel(overlap[overlap[:, 1] >= min_ioff]))
        if cbins.size == 0:
            return []
        chunk_begin, chunk_end = cbins[0], cbins[-1]
        if self.reseek or chunk_begin > tbx.tell():
            tbx.seek(chunk_begin)
            self.walk_buffer[tbx] = []
        elif self.walk_buffer[tbx]:
            if start < self.walk_buffer[tbx][-1].stop:
                remove_i = set()
                for i in range(len(self.walk_buffer[tbx])):
                    if self.walk_buffer[tbx][i].pos > end:
                        break
                    if self.walk_buffer[tbx][i].stop >= start:
                        recs.append(self.walk_buffer[tbx][i])
                    else:
                        remove_i.add(i)
                if remove_i:
                    self.walk_buffer[tbx] = [self.walk_buffer[tbx][i] for i in
                                             range(len(self.walk_buffer[tbx]))
                                             if i not in remove_i]
            else:
                self.walk_buffer[tbx] = []
        if not self.walk_buffer[tbx] or self.walk_buffer[tbx][-1].pos <= end:
            for row in tbx:
                record = self._simplify_cadd_record(row.decode())
                if record is None:
                    continue
                if record.pos > end or tbx.tell() > chunk_end:
                    if use_buffer and record.chrom == chrom:
                        self.walk_buffer[tbx].append(record)
                    break
                if record.stop >= start:
                    recs.append(record)
                    if use_buffer:
                        self.walk_buffer[tbx].append(record)
        self.reseek = not use_buffer
        return recs

    def search_coordinates(self, chrom, start, end):
        hits = []
        for tbx in self.cadd_tabix:
            contig = self.convert_chrom(chrom, tbx)
            if self.walk:
                for rec in self.walk_coordinates(self.bgzfs[tbx],
                                                 contig,
                                                 start,
                                                 end):
                    hits.append(rec)
            else:
                try:
                    for rec in tbx.fetch(contig, start, end):
                        hits.append(self._simplify_cadd_record(rec))
                except ValueError:  # presumably no matching contig
                    pass
        return hits

    def convert_chrom(self, chrom, tbx):
        if chrom.startswith("chr"):
            if not self._has_chr[tbx]:
                return chrom.replace("chr", "", 1)
        elif self._has_chr[tbx]:
            return 'chr' + chrom
        return chrom

    def _get_tabix_files(self, cadd_files):
        for fn in cadd_files:
            idx = fn + '.tbi'
            if not os.path.isfile(idx):  # create index if it doesn't exist
                self.logger.warn("No index found for {} - attempting to index."
                                 .format(fn))
                pysam.tabix_index(fn, preset="vcf")
                self.logger.warn("Finished indexing {}.".format(fn))
            tbx = pysam.TabixFile(fn)
            self.cadd_tabix.append(tbx)
            if self.walk:
                bgzf = pysam.BGZFile(fn)
                self.bgzfs[tbx] = bgzf
                self.indices[bgzf] = read_tbi(idx)

    def _write_for_scoring(self, record, alt):
        if record.DECOMPOSED_ALLELES[alt].ALT != '*':
            self.to_score_file.write("{}\t{}\t.\t{}\t{}\n".format(
                                           record.chrom,
                                           record.DECOMPOSED_ALLELES[alt].POS,
                                           record.DECOMPOSED_ALLELES[alt].REF,
                                           record.DECOMPOSED_ALLELES[alt].ALT))

    def _check_contigs(self):
        tbx_has_chr = dict()
        for tbx in self.cadd_tabix:
            has_chr = False
            no_chr = False
            for c in tbx.contigs:
                if c.startswith('chr'):
                    has_chr = True
                else:
                    no_chr = True
            if has_chr and no_chr:
                raise RuntimeError(
                    "CADD file '{}'".format(tbx.filename.decode()) +
                    "has chromosomes with and without 'chr' prefix - please " +
                    "only provide files with chromosomes in the same format.")
            tbx_has_chr[tbx] = has_chr
        return tbx_has_chr

    def _get_logger(self, logging_level):
        logger = logging.getLogger(__name__)
        if not logger.hasHandlers():
            logger.setLevel(logging_level)
            formatter = logging.Formatter(
                        '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
            ch = logging.StreamHandler()
            ch.setLevel(logger.level)
            ch.setFormatter(formatter)
            logger.addHandler(ch)
        return logger
