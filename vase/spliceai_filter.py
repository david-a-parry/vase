import gzip
import logging
from collections import defaultdict
from .vcf_reader import VcfReader

pre_scored_fields = ["SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG",
                     "DP_AL", "DP_DG", "DP_DL"]
annot_order = ["ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG",
               "DP_AL", "DP_DG", "DP_DL"]


def filter_on_splice_ai(record, min_delta=None, max_delta=None,
                        check_symbol=False, canonical_csq=False):
    keep_alleles = [False] * len(record.record.alts)
    keep_csq = []
    if check_symbol:
        keep_csq = [False] * len(record.CSQ)
    if 'SpliceAI' not in record.record.info:
        return keep_alleles, keep_csq
    info_dicts = list()
    for s in record.record.info['SpliceAI']:
        scores = s.split('|')
        idict = dict((k, v) for k, v in zip(annot_order, scores))
        idict.update(dict((k, float(idict[k])) if idict[k] != '.' else
                          (k, None) for k in annot_order[2:6]))
        info_dicts.append(idict)
    for i in range(1, len(record.record.alleles)):
        if check_symbol:
            csqs = [(n,x) for n,x in enumerate(record.CSQ) if
                    x['alt_index'] == i]
            if canonical_csq:
                csqs = [x for x in csqs if x[1]['CANONICAL'] == 'YES']
        for idict in info_dicts:
            if idict['ALLELE'] != record.record.alleles[i]:
                continue
            if min_delta:
                over_threshold = [ds for ds in annot_order[2:6] if idict[ds] is
                                  not None and idict[ds] >= min_delta]
            if max_delta:
                over_threshold = [ds for ds in annot_order[2:6] if idict[ds] is
                                  not None and idict[ds] <= max_delta]
            if over_threshold:
                keep_alleles[i-1] = True
                if check_symbol:
                    for c in (x for x in csqs if idict['SYMBOL'] ==
                              x[1]['SYMBOL']):
                        keep_csq[c[0]] = True
                else:
                    continue
    return keep_alleles, keep_csq


class SpliceAiFilter(object):
    '''
        An object that annotates or filters VCF records based on SpliceAI
        annotations from SpliceAI score VCF.
    '''

    def __init__(self, vcfs, min_delta=None, max_delta=None, to_delta=None,
                 to_score=None, logging_level=logging.WARNING,
                 no_walk=False, force_walk=False, skip_svs=True):
        '''
            Initialize object with a VCF file and optional filtering
            arguments.

            Args:
                vcfs:       One or more VCFs containing variants to use
                            to filter or annotate records. SpliceAI INFO
                            fields must be present in the format produced
                            for pre-scored variants as downloaded from
                            Jaganathan et al. Cell (2018) or else as
                            generated by the SpliceAI program.

                min_delta:  Return True for each allele if any delta
                            score reaches this value.

                max_delta:  Return True for each allele if any delta
                            score is below this value.

                to_score:   Name of file for writing variants which
                            cannot be found in the given SpliceAI VCFs.
                            If a variant passed to the annotate_or_filter
                            method is not found it will be written to
                            this file in VCF format for scoring with
                            SpliceAI. The scored file can then be used as
                            an additional reference file for scoring
                            variants.

                logging_level:
                            Logging level to use. Default=logging.WARNING

                no_walk:    See VcfFilter documentation.

                force_walk: See VcfFilter documentation.

                skip_svs:   See VcfFilter documentation.

        '''
        self.vcfs = dict()
        self.logger = self._get_logger(logging_level)
        self.walk = not no_walk
        self.force_walk = force_walk
        self.skip_svs = skip_svs
        self.min_delta = min_delta
        self.max_delta = max_delta
        self.to_score = to_score
        self.to_score_file = None
        self.prev_coordinate = (None, -1)
        for vcf in vcfs:
            self.vcfs[vcf] = VcfReader(vcf)
        self.info_fields = {'SpliceAI': {'Number': '.',
                                         'Type': 'String',
                                         'Description': 'SpliceAI variant ' +
                                         'annotation. These include delta ' +
                                         'scores (DS) and delta positions ' +
                                         '(DP) for acceptor gain (AG), ' +
                                         'acceptor loss (AL), donor gain ' +
                                         '(DG), and donor loss (DL). Format:' +
                                         'ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|' +
                                         'DS_DL|DP_AG|DP_AL|DP_DG|DP_DL'}}
        self.vcf_is_prescored = dict()
        self._check_vcf_info()
        if to_score is not None:
            if not to_score.endswith('.gz'):
                to_score += '.gz'
            try:
                from Bio import bgzf
                self.to_score_file = bgzf.BgzfWriter(to_score)
            except ImportError:
                self.logger.warn("Can not import bgzf via biopython. Please " +
                                 "install biopython in order to write bgzip " +
                                 "compressed (.gz/.bgz) output. Missing " +
                                 "SpliceAI scores will be written using gzip" +
                                 " instead.")
                self.to_score_file = gzip.open(to_score, 'wt')
            self.to_score_file.write('##fileformat=VCFv4.0\n' +
                                     '#CHROM\tPOS\tID\tREF\tALT\tQUAL\t' +
                                     'FILTER\tINFO\n')

    def __del__(self):
        if self.to_score_file is not None:
            self.to_score_file.close()

    def _check_vcf_info(self):
       for vcf, vreader in self.vcfs.items():
            if 'SpliceAI' in vreader.header.info:
                if vreader.header.info['SpliceAI'].number == '.':
                   # VCF is in SpliceAI annotated format
                   self.vcf_is_prescored[vcf] = False
            else:  # not annotated by SpliceAI, check if has prescored fields
                for f in pre_scored_fields:
                    if f not in vreader.header.info:
                        raise RuntimeError("ERROR: neither SpliceAI or " +
                                           "individual delta score annotations " +
                                           "found in SpliceAI VCF. Please " +
                                           "ensure you are either using a VCF " +
                                           "annotated using the SpliceAI " +
                                           "program or in the pre-scored format " +
                                           "from the SpliceAI paper")
                self.vcf_is_prescored[vcf] = True

    def get_overlapping_records(self, record):
        '''
            For a given record, returns a list of overlapping records
            in the class's VCFs.
        '''
        overlapping = dict()
        if self.skip_svs and record.IS_SV:
            return overlapping
        if self.walk and not self.force_walk:
            if (record.record.start < self.prev_coordinate[1] and
                    record.recordchrom == self.prev_coordinate[0]):
                self.logger.warn("Input is not sorted by coordinate, will " +
                                 "fall back to slower indvidual index-based " +
                                 "look-ups.")
                self.walk = False
            self.prev_coordinate = (record.record.chrom, record.record.start)
        for vcf, vreader in self.vcfs.items():
            vreader.set_region(record.record.chrom,
                               record.record.start,
                               record.record.stop,
                               walk=self.walk)
            overlapping[vcf] = list(s for s in vreader)
        return overlapping

    def _get_annotation(self, record, alt_index, prescored=False):
        alt = record.alleles[alt_index + 1]
        if prescored:
            info_strings = "|".join("{:.2f}".format(record.info[x])
                                    if isinstance(record.info[x], float)
                                    else str(record.info[x])
                                    for x in pre_scored_fields)
            info_dict = dict([(x, record.info[x]) for x in
                              pre_scored_fields])
        else:
            info_dict = defaultdict(list)
            info_strings = []
            for s in record.info['SpliceAI']:
                scores = s.split('|')
                if scores[0] == alt:
                    info_strings.append('|'.join(scores[1:]))
                    for i in range(1, len(annot_order[:6])):
                        if annot_order[i] == 'SYMBOL':
                            info_dict[annot_order[i]].append(scores[i])
                        else:
                            if scores[i] != '.':
                                info_dict[annot_order[i]].append(
                                    float(scores[i]))
                            else:
                                info_dict[annot_order[i]].append(None)
        return info_strings, info_dict

    def _search_annotations(self, alt_allele, overlaps):
        if self.vcf_is_prescored:
            i_dict = defaultdict(list)
            i_strings = []
        for vcf, olap in overlaps.items():
            for o in olap:
                for i in range(len(o.DECOMPOSED_ALLELES)):
                    if alt_allele == o.DECOMPOSED_ALLELES[i]:
                        if self.vcf_is_prescored[vcf]:
                            s, d = self._get_annotation(
                                o, i, self.vcf_is_prescored[vcf])
                            i_strings.append(s)
                            for k, v in d.items():
                                i_dict[k].append(v)

                        else:
                            return self._get_annotation(
                                o, i, self.vcf_is_prescored[vcf])
        if self.vcf_is_prescored and i_strings:
            return i_strings, i_dict
        return None, None

    def annotate_or_filter(self, record, check_symbol=False,
                           canonical_csq=False):
        '''
            Add SpliceAI annotation for each ALT allele in record. If
            min_delta or max_delta are set return True/False for each
            ALT allele indicating whether they meet the criteria
            specified. Optionally also returns True/False for each CSQ
            annotation from by checking the gene symbol matches in the
            SpliceAI annotation.

            Args:
                record: VCF record to annotate

                check_symbol:
                        Check for annotation's SYMBOL annotation and
                        compare with CSQ annotations. CSQ annotations
                        will only be marked for retaining if these match.

                canonical_csq:
                        Only mark CSQ annotations for retention if
                        CANONICAL field is 'YES' when running with
                        check_symbol option.

        '''
        keep_alleles = [False] * len(record.record.alts)
        keep_csq = []
        if check_symbol:
            try:
                keep_csq = [False] * len(record.CSQ)
            except HeaderError:
                raise RuntimeError("Could not identify CSQ or ANN fields in " +
                                   "VCF header. Please ensure your input is " +
                                   "annotated with Ensembl's VEP")
        overlaps = self.get_overlapping_records(record)
        annotation = []
        for i in range(len(record.DECOMPOSED_ALLELES)):
            info_strings, info_dict = self._search_annotations(
                record.DECOMPOSED_ALLELES[i], overlaps)
            if info_dict is None:
                if self.to_score_file:
                    self._write_for_scoring(record, i)
                continue
            annotation.extend(record.record.alts[i] + "|" + x for x in
                              info_strings)
            over_threshold = []
            if self.min_delta or self.max_delta:
                for ds in annot_order[2:6]:
                    if self.min_delta:
                        over_threshold.extend((x for x in info_dict[ds] if
                                               x is not None and
                                               x >= self.min_delta))
                    if self.max_delta:
                        over_threshold.extend((x for x in info_dict[ds] if
                                               x is not None and
                                               x <= self.max_delta))
                    if over_threshold:
                        keep_alleles[i] = True
                        break
                if not over_threshold:
                    keep_alleles[i] = False
                    continue
            if check_symbol and (self.min_delta or self.max_delta):
                for j in range(len(record.CSQ)):
                    if canonical_csq and record.CSQ[j]['CANONICAL'] != 'YES':
                        continue
                    gene_match = []
                    alt_j = record.CSQ[j]['alt_index'] - 1
                    if alt_j == i:
                        for ds in annot_order[2:]:
                            if self.min_delta:
                                gene_match.extend((x for x in info_dict[ds] if
                                                   x >= self.min_delta and
                                                   info_dict['SYMBOL'] ==
                                                   record.CSQ[j]['SYMBOL']))
                            if self.max_delta:
                                gene_match.extend((x for x in info_dict[ds] if
                                                   x <= self.max_delta and
                                                   info_dict['SYMBOL'] ==
                                                   record.CSQ[j]['SYMBOL']))
                            if gene_match:
                                break
                    if gene_match:
                        keep_csq[j] = True
        if annotation:
            record.add_info_fields({'SpliceAI': ",".join(annotation)})
        return keep_alleles, keep_csq

    def _write_for_scoring(self, record, alt):
        if record.DECOMPOSED_ALLELES[alt].ALT != '*':
            self.to_score_file.write("{}\t{}\t.\t{}\t{}\t.\t.\t.\n".format(
                                           record.record.chrom,
                                           record.DECOMPOSED_ALLELES[alt].POS,
                                           record.DECOMPOSED_ALLELES[alt].REF,
                                           record.DECOMPOSED_ALLELES[alt].ALT))

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
