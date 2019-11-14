import os
import pysam
import logging
import gzip
from collections import defaultdict

class CaddFilter(object):
    '''
        An object that filters/annotates VCF records using CADD PHRED
        scores provided by at least one tabix indexed file of CADD
        scores.
    '''

    def __init__(self, cadd_files=[], cadd_dir=[], min_phred=None,
                 min_raw_score=None, to_score=None,
                 logging_level=logging.WARNING):
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

        '''
        self.to_score_file = None
        self.logger = self._get_logger(logging_level)
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
        self.cadd_tabix = self._get_tabix_files(cadd_files)
        self._has_chr = self._check_contigs()
        self.phred = min_phred
        self.raw = min_raw_score
        self.info_fields = {'CADD_PHRED_score': {'Number': 'A',
                                                 'Type': 'Float',
                                                 'Description': '"CADD PHRED ' +
                                                                'score added '+
                                                                'from ' +
                                                                'reference ' +
                                                                'files by ' + 
                                                                'VASE"'
                                                },
                            'CADD_raw_score': {'Number': 'A',
                                               'Type': 'Float',
                                               'Description': '"CADD RawScore' +
                                                              ' added from ' +
                                                              'reference '+
                                                              'files by VASE"'
                                              },
                           }
        if to_score is not None:
            if not to_score.endswith('.gz'):
                to_score += '.gz'
            self.to_score_file = gzip.open(to_score, 'wt')

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
        scores = self.score_record(record)
        info_to_add = defaultdict(list)
        filter_alleles = []
        i = 0
        for s in scores:
            info_to_add['CADD_raw_score'].append(s[0] or '.')
            info_to_add['CADD_PHRED_score'].append(s[1] or '.')
            do_filter = False
            if self.raw and s[0] is not None:
                if float(s[0]) < self.raw:
                    do_filter = True
            if self.phred and s[1] is not None:
                if float(s[1]) < self.phred:
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
        start = record.POS - 1
        end = record.SPAN
        hits = self.search_coordinates(record.CHROM, start, end)
        scores = []
        for i in range(len(record.DECOMPOSED_ALLELES)):
            s = (None, None)
            for h in (x for x in hits if x is not None):
                (pos, ref, alt, raw, phred) = h
                if (record.DECOMPOSED_ALLELES[i].POS == pos and
                    record.DECOMPOSED_ALLELES[i].REF == ref and
                    record.DECOMPOSED_ALLELES[i].ALT == alt):
                    s = (raw, phred)
                    break #bail on first matching variant
            scores.append(s)
        return scores

    def _simplify_cadd_record(self, cadd):
        '''
            Return position, ref allele, alt allele, raw and Phred score
            after reducing alleles to their most simple representation.
        '''
        cols = cadd.split("\t")
        if len(cols) < 6:
            self.logger.warn("Not enought columns for CADD record: {}"
                             .format(cadd))
            return None
        ref = cols[2]
        alt = cols[3]
        pos = int(cols[1])
        while len(ref) > 1 and len(alt) > 1:
            if ref[-1] == alt[-1]:               #remove identical suffixes
                ref = ref[:-1]
                alt = alt[:-1]
            else:
                break
        while len(ref) > 1 and len(alt) > 1:
            if ref[0] == alt[0]:                 #remove identical prefixes
                ref = ref[1:]
                alt = alt[1:]
                pos += 1
            else:
                break
        return (pos, ref, alt, cols[4], cols[5])


    def search_coordinates(self, chrom, start, end):
        hits = []
        for tbx in self.cadd_tabix:
            chrom = str(chrom)
            if chrom.startswith("chr"):
                if not self._has_chr[tbx]:
                    chrom = chrom.replace("chr", "", 1)
            elif self._has_chr[tbx]:
                chrom = 'chr' + chrom
            try:
                 for rec in tbx.fetch(chrom, start, end):
                    hits.append(self._simplify_cadd_record(rec))
            except ValueError: #presumably no matching contig
                pass
        return hits

    def _get_tabix_files(self, cadd_files):
        tabixfiles = []
        for fn in cadd_files:
            idx = fn + '.tbi'
            if not os.path.isfile(idx):   #create index if it doesn't exist
                self.logger.warn("No index found for {} - attempting to index."
                                 .format(fn))
                pysam.tabix_index(fn, preset="vcf")
                self.logger.warn("Finished indexing {}.".format(fn))
            tabixfiles.append(pysam.Tabixfile(fn))
        return tabixfiles

    def _write_for_scoring(self, record, alt):
        if record.DECOMPOSED_ALLELES[alt].ALT != '*':
            self.to_score_file.write("{}\t{}\t.\t{}\t{}\n".format(
                                           record.CHROM,
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
                raise RuntimeError("CADD file '{}'".format(
                                                        tbx.filename.decode())+
                                   "has chromosomes with and without 'chr' " +
                                   "prefix - please only provide files with " +
                                   "chromosomes in the same format.")
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
