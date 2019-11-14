import sys
import re
import logging
import io
from parse_vcf import VcfReader
from .dbsnp_filter import dbSnpFilter, clinvar_path_annot
from .gnomad_filter import GnomadFilter
from .vcf_filter import VcfFilter
from .vep_filter import VepFilter
from .cadd_filter import CaddFilter
from .sample_filter import SampleFilter
from .ped_file import PedFile, Individual, PedError
from .family_filter import FamilyFilter, ControlFilter
from .family_filter import RecessiveFilter, DominantFilter, DeNovoFilter
from .burden_counter import BurdenCounter
from .var_by_region import VarByRegion
from .region_iter import RegionIter
from .gt_annotator import GtAnnotator
from .spliceai_filter import SpliceAiFilter, filter_on_splice_ai
from .info_filter import InfoFilter
from .g2p import G2P

class VaseRunner(object):

    def __init__(self, args):

        self._twirler = ['-', '\\', '|', '/', ]
        self.args = args
        self._set_logger()
        self.input = VcfReader(self.args.input)
        self.var_stream = self.input.parser
        self.keep_filters = None
        self.exclude_filters = None
        if args.keep_filters:
            self.keep_filters = set(args.keep_filters)
        if args.exclude_filters:
            self.exclude_filters = set(args.exclude_filters)
        self.var_types = self._parse_var_type_arg()
        self.info_filter = self._parse_info_filters()
        self.prev_cadd_phred = False
        self.prev_cadd_raw = False
        self.prev_splice_ai = False
        self._get_prev_annotations()
        self.out = self.get_output()
        self.vcf_filters = self.get_vcf_filter_classes()
        self.cadd_filter = self.get_cadd_filter()
        self.splice_ai_filter = self.get_splice_ai_filter()
        self.gt_annotators = self.get_gt_annotators()
        self.ped = None
        if args.ped:
            self.ped = PedFile(args.ped)
        self.csq_filter = None
        self.g2p = None
        self.gene_filter = None
        self.retrieving_by_region = False
        self.gt_args = dict(gq=args.gq, dp=args.dp, max_dp=args.max_dp,
                            het_ab=args.het_ab, hom_ab=args.hom_ab,
                            min_control_dp=args.control_dp,
                            max_control_dp=args.control_max_dp,
                            min_control_gq=args.control_gq,
                            control_het_ab=args.control_het_ab,
                            control_hom_ab=args.control_hom_ab,
                            con_ref_ab=args.control_max_ref_ab,
                            sv_gq=args.sv_gq, sv_dp=args.sv_dp,
                            sv_max_dp=args.sv_max_dp, sv_het_ab=args.sv_het_ab,
                            sv_hom_ab=args.sv_hom_ab,
                            sv_min_control_dp=args.sv_control_dp,
                            sv_max_control_dp=args.sv_control_max_dp,
                            sv_min_control_gq=args.sv_control_gq,
                            sv_control_het_ab=args.sv_control_het_ab,
                            sv_control_hom_ab=args.sv_control_hom_ab,
                            sv_con_ref_ab=args.sv_control_max_ref_ab,
                            del_dhffc=args.duphold_del_dhffc,
                            dup_dhbfc=args.duphold_dup_dhbfc,
                            control_del_dhffc=args.control_duphold_del_dhffc,
                            control_dup_dhbfc=args.control_duphold_dup_dhbfc,
                            )
        #region, bed and gene_bed args are mutually exclusive (handled by parser)
        if args.region is not None:
            self.var_stream = VarByRegion(self.input,
                                          region_iter=RegionIter(args.region),
                                          stream=args.stream,
                                          exclude=args.exclude_regions)
            self.retrieving_by_region = True
        if args.bed is not None:
            self.logger.info("Reading, sorting and merging intervals in " +
                             "{}".format(args.bed))
            self.var_stream = VarByRegion(self.input, bed=args.bed,
                                          stream=args.stream,
                                          exclude=args.exclude_regions)
            self.retrieving_by_region = True
            self.logger.info("Finished processing intervals.")
        if args.gene_bed is not None:
            self.logger.info("Reading, sorting and merging intervals in " +
                             "{}".format(args.gene_bed))
            self.gene_filter = VarByRegion(self.input, bed=args.gene_bed,
                                           gene_targets=True,
                                           stream=args.stream,
                                           exclude=args.exclude_regions)
            if args.csq is None:
                args.csq = ['all']
            if args.biotypes is None:
                args.biotypes = ['all']
            self.var_stream = self.gene_filter
            self.retrieving_by_region = True
            self.logger.info("Finished processing intervals.")
        if args.g2p is not None:
            self.g2p = G2P(args.g2p)
        if (args.csq is not None or args.impact is not None or self.args.g2p is
            not None):
            if args.no_vep_freq:
                if args.vep_af:
                    self.logger.warn("Ignoring --vep_af argument because " +
                                     "--no_vep_af argument is in use.")
                vep_freq = None
                vep_min_freq = None
                vep_af = []
            else:
                vep_freq = args.freq
                vep_min_freq = args.min_freq
                vep_af = args.vep_af
            self.csq_filter = VepFilter(
                    vcf=self.input,
                    csq=args.csq,
                    impact=args.impact,
                    canonical=args.canonical,
                    biotypes=args.biotypes,
                    in_silico=args.missense_filters,
                    filter_unpredicted=args.filter_unpredicted,
                    keep_any_damaging=args.keep_if_any_damaging,
                    loftee=args.loftee,
                    splice_in_silico=args.splice_filters,
                    splice_filter_unpredicted=args.splice_filter_unpredicted,
                    splice_keep_any_damaging=args.splice_keep_if_any_damaging,
                    retain_labels=args.retain_labels,
                    filter_flagged_features=args.flagged_features,
                    freq=vep_freq,
                    min_freq=vep_min_freq,
                    filter_known=self.args.filter_known,
                    filter_novel=self.args.filter_novel,
                    afs=vep_af,
                    gene_filter=self.gene_filter,
                    blacklist=args.feature_blacklist,
                    pathogenic=args.pathogenic,
                    no_conflicted=args.no_conflicted,
                    g2p=self.g2p,
                    check_g2p_consequence=self.args.check_g2p_consequence,
                    logging_level=self.logger.level)
        self.sample_filter = None
        self.burden_counter = None
        if args.burden_counts:
            if args.n_cases or args.n_controls:
                self.logger.warn("--n_cases and --n_controls arguments are " +
                                 "ignored when using --burden_counter.")
            self.burden_counter = BurdenCounter(self.input, args.burden_counts,
                                                is_gnomad=args.gnomad_burden,
                                                cases=args.cases,
                                                controls=args.controls,
                                                gq=args.gq, dp=args.dp,
                                                max_dp=args.max_dp,
                                                het_ab=args.het_ab,
                                                hom_ab=args.hom_ab,)
        elif args.cases or args.controls:
            self.sample_filter = SampleFilter(
                                    self.input, cases=args.cases,
                                    controls=args.controls,
                                    n_cases=args.n_cases,
                                    n_controls=args.n_controls,
                                    confirm_missing=args.confirm_control_gts,
                                    **self.gt_args)
        self.de_novo_filter = None
        self.dominant_filter = None
        self.recessive_filter = None
        self.family_filter = None
        self.control_filter = None
        self.variant_cache = VariantCache()
        self.use_cache = False
        self.prog_interval = args.prog_interval
        self.log_progress = args.log_progress
        self.strict_recessive_inheritance = args.strict_recessive
        self.report_fhs = self.get_report_filehandles()
        seg_info = list()
        if args.de_novo:
            seg_info.extend(self._get_de_novo_filter())
        if args.biallelic or args.singleton_recessive:
            seg_info.extend(self._get_recessive_filter())
        if args.dominant or args.singleton_dominant:
            seg_info.extend(self._get_dominant_filter())
        if args.check_g2p_inheritance:
            if (not self.dominant_filter and not self.recessive_filter and
                    not self.de_novo_filter):
                raise RuntimeError("--check_g2p_inheritance option requires " +
                                   "at least one segregation option (e.g. " +
                                   "--recessive, --dominant or --de_novo) " + 
                                   "to be selected.")
        self._check_got_inherit_filter()
        self._set_seg_annot_cleanup(seg_info)
        self.var_written = 0
        self.var_filtered = 0

    def run(self):
        ''' Run VCF filtering/annotation using args from bin/vase'''
        self.logger.info('Starting variant processing')
        self.print_header()
        var_count = 0
        prog_updates = 0
        prog_string = ''
        for record in self.var_stream:
            self.process_record(record)
            var_count += 1
            if (not self.args.no_progress and
                    var_count % self.prog_interval == 0):
                n_prog_string = ('{:,} variants processed, '.format(var_count) +
                               '{:,} filtered, {:,} written... at pos {}:{}'
                               .format(self.var_filtered, self.var_written,
                                       record.CHROM, record.POS))
                if self.retrieving_by_region and self.var_stream.region_iter:
                    n_prog_string += " (processing region {}/{})".format(
                        self.var_stream.region_iter.current_index + 1,
                        len(self.var_stream.region_iter.intervals))
                if self.log_progress:
                    self.logger.info(n_prog_string)
                else:
                    twirl = self._twirler[prog_updates % 4]
                    n_prog_string = '\r' + n_prog_string + ' ' + twirl
                    if len(prog_string) > len(n_prog_string):
                        sys.stderr.write('\r' + ' ' * len(prog_string) )
                    prog_string = n_prog_string
                    sys.stderr.write(prog_string)
                prog_updates += 1
        self.finish_up()
        if prog_string:
            sys.stderr.write('\r' + '-' * len(prog_string) + '\n')
        self.logger.info('Finished processing {:,} {}.'
                             .format(var_count, self._var_or_vars(var_count)))
        self.logger.info('{:,} {} filtered.'
                         .format(self.var_filtered,
                                 self._var_or_vars(self.var_filtered)))
        self.logger.info('{:,} {} written.'
                         .format(self.var_written,
                                 self._var_or_vars(self.var_written)))
        if self.out is not sys.stdout:
            self.out.close()

    def process_record(self, record):
        if self.filter_global(record):
            self.var_filtered += 1
            return
        filter_alleles = None
        if self.var_types:
            filter_alleles = [x.var_type not in self.var_types for x in
                              record.DECOMPOSED_ALLELES]
            if all(filter_alleles):
                #no ALT matches any variant type asked for
                self.var_filtered += 1
                return
        filter_alleles, filter_csq = self.filter_alleles_external(record,
                                                                filter_alleles)
        if all(filter_alleles):
            #all alleles should be filtered
            self.var_filtered += 1
            return
        if self.sample_filter:
            for i in range(1, len(record.ALLELES)):
                r = self.sample_filter.filter(record, i)
                if r:
                    filter_alleles[i-1] = True
                if all(filter_alleles):
                    #all alleles should be filtered
                    self.var_filtered += 1
                    return
        dom_filter_alleles = list(filter_alleles)
        if self.control_filter:
            for i in range(1, len(record.ALLELES)):
                if dom_filter_alleles[i-1]: #no need to filter again
                    continue
                r = self.control_filter.filter(record, i)
                if r:
                    dom_filter_alleles[i-1] = True
        self.remove_previous_inheritance_filters(record)
        denovo_hit = False
        dom_hit = False
        recessive_hit = False
        if self.dominant_filter:
            dom_hit = self.dominant_filter.process_record(record,
                                                          dom_filter_alleles,
                                                          filter_csq)
        if self.de_novo_filter:
            denovo_hit = self.de_novo_filter.process_record(record,
                                                            dom_filter_alleles,
                                                            filter_csq)
        if self.recessive_filter:
            recessive_hit = self.recessive_filter.process_record(record,
                                                    filter_alleles, filter_csq)
        if self.use_cache:
            if denovo_hit or dom_hit or recessive_hit:
                keep_record_anyway = False
                if self.args.min_families < 2:
                    keep_record_anyway = denovo_hit or dom_hit
                self.variant_cache.add_record(record, keep_record_anyway)
            else:
                self.variant_cache.check_record(record)
                self.var_filtered += 1
            if self.variant_cache.output_ready:
                self.output_cache()
        elif self.de_novo_filter or self.dominant_filter:
            if denovo_hit or dom_hit:
                self.output_record(record)
                if self.burden_counter:
                    #getting relevant alleles and feats is a bit of a fudge using
                    #annotations added by dom/denovo filter
                    b_filt_al = [True] * (len(record.ALLELES) -1)
                    b_filt_csq = [[True] * len(record.CSQ)] * len(b_filt_al)
                    if dom_hit:
                        b_filt_al, b_filt_csq = self._seg_alleles_from_record(
                                                  record,
                                                  self.dominant_filter.prefix,
                                                  b_filt_al,
                                                  b_filt_csq)
                    if denovo_hit:
                        b_filt_al, b_filt_csq = self._seg_alleles_from_record(
                                                  record,
                                                  self.de_novo_filter.prefix,
                                                  b_filt_al,
                                                  b_filt_csq)
                    for i in range(len(b_filt_al)):
                        if not b_filt_al[i]:
                            feat = (record.CSQ[j]['Feature'] for j in
                                    range(len(record.CSQ)) if not
                                    b_filt_csq[i][j])
                            self.burden_counter.count_samples(record,
                                                              feat, i, 1)
            else:
                self.var_filtered += 1
        else:
            if all(filter_alleles):
                self.var_filtered += 1
                return
            if self.burden_counter:
                self.burden_counter.count(record, filter_alleles, filter_csq)
            self.output_record(record)

    def output_record(self, record):
        for gt_anno in self.gt_annotators:
            gt_anno.annotate(record)
        self.out.write(str(record) + '\n')
        self.var_written += 1

    def output_cache(self, final=False):
        keep_ids = set()
        burden_vars = dict() #dict of inheritance model to segregants
        if self.recessive_filter:
            vid_to_seg = self.recessive_filter.process_potential_recessives(
                                                                   final=final)
            keep_ids.update(vid_to_seg.keys())
            if self.burden_counter:
                burden_vars['recessive'] = vid_to_seg
        if self.dominant_filter and self.args.min_families > 1:
            vid_to_seg = self.dominant_filter.process_dominants(final=final)
            keep_ids.update(vid_to_seg.keys())
            if self.burden_counter:
                burden_vars['dominant'] = vid_to_seg
        if self.de_novo_filter and self.args.min_families > 1:
            vid_to_seg = self.de_novo_filter.process_de_novos(final=final)
            keep_ids.update(vid_to_seg.keys())
            if self.burden_counter:
                burden_vars['de_novo'] = vid_to_seg
        if final:
            self.variant_cache.add_cache_to_output_ready()
        for var in self.variant_cache.output_ready:
            if var.can_output or var.var_id in keep_ids:
                self.output_record(var.record)
            else:
                self.var_filtered += 1
        if self.burden_counter:
            self._burden_from_cache(burden_vars)
        self.variant_cache.output_ready = []

    def _burden_from_cache(self, model_to_vars):
        for model in ('dominant', 'de_novo'):
            if model in model_to_vars:
                for v in model_to_vars[model]:
                    for seg in model_to_vars[model][v]:
                        self.burden_counter.count_samples(seg.record,
                                                          seg.features,
                                                          seg.allele - 1, 1)
        # we count recessives last as the max per allele (2) is higher than for
        # dom/de novos (otherwise max per sample could get lowered to 1)
        if 'recessive' in model_to_vars:
            for v in model_to_vars['recessive']:
                for seg in model_to_vars['recessive'][v]:
                    self.burden_counter.count_samples(seg.record, seg.features,
                                                      seg.allele - 1, 2)

    def finish_up(self):
        if self.use_cache:
            self.output_cache(final=True)
        if self.burden_counter:
            self.burden_counter.output_counts()
        for fh in self.report_fhs.values():
            if fh is not None:
                fh.close()

    def filter_alleles_external(self, record, remove_alleles=None):
        '''
            Return True or False for each allele indicating whether an
            allele should be filtered based on information from VEP,
            dbSNP, ClinVar or gnomAD.
        '''
        # remove_alleles indicates whether allele should be filtered;
        # keep_alleles indicates whether allele should be kept, overriding any
        # indications in remove_alleles (e.g. if labelled pathogenic in
        # ClinVar)
        # remove_csq indicates for each VEP CSQ whether that CSQ should be
        # ignored
        if not remove_alleles:
            remove_alleles = [False] * (len(record.ALLELES) -1)
        keep_alleles = [False] * (len(record.ALLELES) -1)
        matched_alleles = [False] * (len(record.ALLELES) -1)
        remove_csq = None
        #if allele is '*' should be set to filtered
        for i in range(1, len(record.ALLELES)):
            if record.ALLELES[i] == '*':
                remove_alleles[i-1] = True
        #filter on provided INFO field filters
        if self.info_filter:
            r_alts = self.info_filter.filter(record)
            self._set_to_true_if_true(remove_alleles, r_alts)
            if all(remove_alleles):
                # bail out now if no valid allele and not keeping clinvar
                return remove_alleles, remove_csq
        #check VCF's internal AF
        if self.args.af or self.args.min_af:
            r_alts = self.filter_on_af(record)
            self._set_to_true_if_true(remove_alleles, r_alts)
            if (not self.args.clinvar_path and all(remove_alleles)):
                # bail out now if no valid allele and not keeping clinvar
                # path variants - if using clinvar path we have to ensure we
                # haven't got a path variant with a non-qualifying allele
                return remove_alleles, remove_csq
        #check VCF's internal AC
        if self.args.ac or self.args.min_ac:
            r_alts = self.filter_on_ac(record)
            self._set_to_true_if_true(remove_alleles, r_alts)
            if (not self.args.clinvar_path and all(remove_alleles)):
                # bail out now if no valid allele and not keeping clinvar
                return remove_alleles, remove_csq
        #check functional consequences
        if self.csq_filter:
            r_alts, remove_csq = self.csq_filter.filter(record)
            self._set_to_true_if_true(remove_alleles, r_alts)
            if self.prev_splice_ai and (self.args.splice_ai_min_delta or
                                        self.args.splice_ai_max_delta):
                splice_alleles, splice_csq = filter_on_splice_ai(record,
                                    min_delta=self.args.splice_ai_min_delta,
                                    max_delta=self.args.splice_ai_max_delta,
                                    check_symbol=True,
                                    canonical_csq=self.args.canonical)
                self._set_to_false_if_true(remove_alleles, splice_alleles)
                self._set_to_false_if_true(remove_csq, splice_csq)
            if self.splice_ai_filter:
                splice_alleles, splice_csq = (
                self.splice_ai_filter.annotate_or_filter(record,
                                                          True,
                                                          self.args.canonical))
                if (self.args.splice_ai_min_delta
                    or self.args.splice_ai_max_delta):
                    #RETAIN Alleles/csq if SpliceAI scores meet threshold
                    self._set_to_false_if_true(remove_alleles, splice_alleles)
                    self._set_to_false_if_true(remove_csq, splice_csq)
            if (not self.args.clinvar_path and all(remove_alleles)):
                # bail out now if no valid consequence
                return remove_alleles, remove_csq
        elif (self.splice_ai_filter or self.args.splice_ai_min_delta or
              self.args.splice_ai_max_delta):
            if self.prev_splice_ai and (self.args.splice_ai_min_delta or
                                        self.args.splice_ai_max_delta):
                splice_alleles, splice_csq = filter_on_splice_ai(record,
                                    min_delta=self.args.splice_ai_min_delta,
                                    max_delta=self.args.splice_ai_max_delta,
                                    check_symbol=True,
                                    canonical_csq=self.args.canonical)
                self._set_to_false_if_true(remove_alleles, splice_alleles)
            else:
                splice_alleles, splice_csq = (
                              self.splice_ai_filter.annotate_or_filter(record))
                if (self.args.splice_ai_min_delta
                    or self.args.splice_ai_max_delta):
                    #use SpliceAI on alleles only - FILTER If threshold not met
                    #FILTER alleles if SpliceAI scores meet threshold
                    self._set_to_false_if_true(remove_alleles, splice_alleles)
        if self.prev_cadd_phred and self.args.cadd_phred:
            r_alts = self.filter_on_existing_cadd_phred(record)
            self._set_to_true_if_true(remove_alleles, r_alts)
        if self.prev_cadd_raw and self.args.cadd_raw:
            r_alts = self.filter_on_existing_cadd_raw(record)
            self._set_to_true_if_true(remove_alleles, r_alts)
        if self.cadd_filter:
            r_alts = self.cadd_filter.annotate_or_filter(record)
            self._set_to_true_if_true(remove_alleles, r_alts)
            if (not self.args.clinvar_path and all(remove_alleles)):
                # bail out now if no valid consequence and not keeping clinvar
                # path variants - if using clinvar path we have to ensure we
                # haven't got a path variant with a non-qualifying consequence
                return remove_alleles, remove_csq
        for f in self.vcf_filters:
            r, k, m = f.annotate_and_filter_record(record)
            # should only overwrite value of remove_alleles[i] or
            # keep_alleles[i] with True, not False (e.g. if already set to
            # be filtered because of a freq in ExAC we shouldn't set to
            # False just because it is absent from dbSNP)
            self._set_to_true_if_true(remove_alleles, r)
            self._set_to_true_if_true(matched_alleles, m)
            self._set_to_true_if_true(keep_alleles, k)
        if self.prev_freqs:
            r,m = self.filter_on_existing_freq(record)
            self._set_to_true_if_true(remove_alleles, r)
            self._set_to_true_if_true(matched_alleles, m)
        if self.prev_homs:
            r,m = self.filter_on_existing_homs(record)
            self._set_to_true_if_true(remove_alleles, r)
            self._set_to_true_if_true(matched_alleles, m)
        if self.prev_builds:
            r,m = self.filter_on_existing_build(record)
            self._set_to_true_if_true(remove_alleles, r)
            self._set_to_true_if_true(matched_alleles, m)
        if self.prev_clinvar:
            k,m = self.filter_on_existing_clnsig(record)
            self._set_to_true_if_true(matched_alleles, m)
            self._set_to_true_if_true(keep_alleles, k)
        verdict = []
        for i in range(len(remove_alleles)):
            if keep_alleles[i]:
                verdict.append(False)
            elif remove_alleles[i]:
                verdict.append(True)
            elif self.args.filter_known and matched_alleles[i]:
                verdict.append(True)
            elif self.vcf_filters:
                #only apply filter_novel if we have supplied an external vcf
                if self.args.filter_novel and not matched_alleles[i]:
                    verdict.append(True)
                else:
                    verdict.append(False)
            else:
                verdict.append(False)
        return verdict, remove_csq

    def filter_on_af(self, record):
        remove  = [False] * (len(record.ALLELES) -1)
        if self.args.filtering_an and self.an_below_threshold(record):
            return remove
        try:
            af = record.parsed_info_fields(fields=['AF'])['AF']
            for i in range(len(remove)):
                if self.args.af:
                    if af[i] is not None and af[i] > self.args.af:
                        remove[i] = True
                if self.args.min_af:
                    if af[i] is None or af[i] < self.args.min_af:
                        remove[i] = True
        except KeyError:
            self.logger.debug("No 'AF' in INFO at {}:{}".format(record.CHROM,
                                                               record.POS))
            if 'AN' not in record.parsed_info_fields(fields=['AN']):
                self.logger.warn("No 'AF' or 'AN' in INFO at {}:{}".format(
                        record.CHROM, record.POS) + " - will not filter on AF")
            elif 'AC' not in record.parsed_info_fields(fields=['AC']):
                self.logger.warn("No 'AF' or 'AC' in INFO at {}:{}".format(
                        record.CHROM, record.POS) + " - will not filter on AF")
            else:
                self.logger.debug("Trying AC/AN instead")
                return self.filter_on_ac_over_an(record)
        return remove

    def an_below_threshold(self, record):
        try:
            return (record.parsed_info_fields(fields=['AN'])['AN'] <
                    self.args.filtering_an)
        except KeyError:
            self.logger.warn("No 'AN' in INFO at {}:{}".format(record.CHROM,
                                                               record.POS))
            return True

    def an_under_minimum(self, record):
        try:
            return (record.parsed_info_fields(fields=['AN'])['AN'] <
                    self.args.min_an)
        except KeyError:
            self.logger.warn("No 'AN' in INFO at {}:{}".format(record.CHROM,
                                                               record.POS))
            return True

    def filter_on_ac_over_an(self, record):
        remove  = [False] * (len(record.ALLELES) -1)
        try:
            info = record.parsed_info_fields(fields=['AC', 'AN'])
            ac = info['AC']
            an = info['AN']
            af = [ac[i]/an if an > 0 else 0 for i in range(len(ac))]
            for i in range(len(remove)):
                if self.args.af:
                    if af[i] is not None and af[i] > self.args.af:
                        remove[i] = True
                if self.args.min_af:
                    if af[i] is None or af[i] < self.args.min_af:
                        remove[i] = True
        except KeyError:
            self.logger.warn("Missing 'AN' or 'AC' field in INFO at {}:{}"
                             .format(record.CHROM, record.POS) + " - will not"
                             " filter on AF")
        return remove

    def filter_on_ac(self, record):
        remove  = [False] * (len(record.ALLELES) -1)
        try:
            ac = record.parsed_info_fields(fields=['AC'])['AC']
            for i in range(len(remove)):
                if self.args.ac:
                    if ac[i] is not None and ac[i] > self.args.ac:
                        remove[i] = True
                if self.args.min_ac:
                    if ac[i] is None or ac[i] < self.args.min_ac:
                        remove[i] = True
        except KeyError:
            self.logger.warn("No 'AC' in INFO at {}:{}".format(record.CHROM,
                                                               record.POS))
        return remove

    def filter_on_existing_cadd_phred(self, record):
        remove  = [False] * (len(record.ALLELES) -1)
        phreds = record.parsed_info_fields(fields=['CADD_PHRED_score'])
        if 'CADD_PHRED_score' in phreds:
            for i in range(len(remove)):
                if (phreds['CADD_PHRED_score'][i] is not None and
                    phreds['CADD_PHRED_score'][i] < self.args.cadd_phred):
                    remove[i] = True
        return remove

    def filter_on_existing_cadd_raw(self, record):
        remove  = [False] * (len(record.ALLELES) -1)
        raws = record.parsed_info_fields(fields=['CADD_raw_score'])
        if 'CADD_raw_score' in raws:
            for i in range(len(remove)):
                if (raws['CADD_raw_score'][i] is not None and
                    raws['CADD_raw_score'][i] < self.args.cadd_raw):
                    remove[i] = True
        return remove

    def filter_on_existing_freq(self, record):
        remove  = [False] * (len(record.ALLELES) -1)
        matched = [False] * (len(record.ALLELES) -1)
        parsed = record.parsed_info_fields(fields=self.prev_freqs)
        for annot in parsed:
            if parsed[annot] is None:
                continue
            for i in range(len(remove)):
                if parsed[annot][i] is not None:
                    matched[i] = True
                    if self.args.freq:
                        if parsed[annot][i] >= self.args.freq:
                            remove[i] = True
                    if self.args.min_freq:
                        if parsed[annot][i] < self.args.min_freq:
                            remove[i] = True
        return remove,matched

    def filter_on_existing_homs(self, record):
        remove  = [False] * (len(record.ALLELES) -1)
        matched = [False] * (len(record.ALLELES) -1)
        parsed = record.parsed_info_fields(fields=self.prev_homs)
        for annot in parsed:
            if parsed[annot] is None:
                continue
            for i in range(len(remove)):
                if parsed[annot][i] is not None:
                    matched[i] = True
                    if parsed[annot][i] > self.args.max_gnomad_homozygotes:
                        remove[i] = True
        return remove,matched

    def filter_on_existing_build(self, record):
        remove  = [False] * (len(record.ALLELES) -1)
        matched = [False] * (len(record.ALLELES) -1)
        parsed = record.parsed_info_fields(fields=self.prev_builds)
        for annot in parsed:
            if parsed[annot] is None:
                continue
            for i in range(len(remove)):
                if parsed[annot][i] is not None:
                    matched[i] = True
                    if self.args.build:
                        if parsed[annot][i] <= self.args.build:
                            remove[i] = True
                    if self.args.max_build:
                        if parsed[annot][i] > self.args.max_build:
                            remove[i] = True
        return remove,matched

    def filter_on_existing_clnsig(self, record):
        keep    = [False] * (len(record.ALLELES) -1)
        matched = [False] * (len(record.ALLELES) -1)
        parsed = record.parsed_info_fields(fields=self.prev_clinvar)
        for annot in parsed:
            if parsed[annot] is None:
                continue
            for i in range(len(keep)):
                if parsed[annot][i] is not None:
                    matched[i] = True
                    if self.args.clinvar_path:
                        if any(x for x in clinvar_path_annot if x in
                                parsed[annot][i].split('|')):
                                keep[i] = True
        return keep,matched

    def filter_global(self, record):
        ''' Return True if record fails any global variant filters.'''
        if self.args.pass_filters:
            if record.FILTER != 'PASS':
                return True
        elif self.keep_filters: #--pass and --keep_filters are mutually excl.
            fs = record.FILTER.split(';')
            if not self.keep_filters.issuperset(fs):
                return True
        if self.exclude_filters:
            fs = record.FILTER.split(';')
            if self.exclude_filters.intersection(fs):
                return True
        if self.args.variant_quality is not None:
            if record.QUAL < self.args.variant_quality:
                return True
        if self.args.max_alt_alleles is not None:
            if (len([x for x in record.ALLELES if x != '*']) >
                self.args.max_alt_alleles + 1):
                return True
        if self.args.min_an and self.an_under_minimum(record):
            return True
        if self.args.filter_asterisk_only_calls:
            if len(record.ALLELES) == 2 and record.ALLELES[1] == '*':
                return True
        return False

    def get_cadd_filter(self):
        if self.args.cadd_directory or self.args.cadd_files:
            cadd_args = {'cadd_files': self.args.cadd_files,
                         'cadd_dir': self.args.cadd_directory,
                         'min_phred': self.args.cadd_phred,
                         'min_raw_score': self.args.cadd_raw,
                         'to_score': self.args.missing_cadd_scores}
            cf = CaddFilter(**cadd_args)
            for f,d in cf.info_fields.items():
                self.logger.debug("Adding {} annotation")
                self.input.header.add_header_field(name=f, dictionary=d,
                                          field_type='INFO')
            return cf
        else:
            if self.args.cadd_phred is not None and not self.prev_cadd_phred:
                raise RuntimeError("--cadd_phred cutoff given but VCF does " +
                                   "not have CADD_PHRED_score INFO field " +
                                   "in its header and no CADD score files " +
                                   "have been specified with --cadd_files or" +
                                   " --cadd_dir arguments.")
            if self.args.cadd_raw is not None and not self.prev_cadd_raw:
                raise RuntimeError("--cadd_raw cutoff given but VCF does " +
                                   "not have CADD_raw_score INFO field " +
                                   "in its header and no CADD score files " +
                                   "have been specified with --cadd_files " +
                                   "or --cadd_dir arguments.")
        return None

    def get_splice_ai_filter(self):
        if self.args.splice_ai_vcfs:
            splice_ai_args = {'vcfs': self.args.splice_ai_vcfs,
                              'min_delta': self.args.splice_ai_min_delta,
                              'max_delta': self.args.splice_ai_max_delta,
                              'to_score': self.args.missing_splice_ai_scores,
                              'logging_level': self.logger.level,}
            sf = SpliceAiFilter(**splice_ai_args)
            for f,d in sf.info_fields.items():
                self.logger.debug("Adding {} annotation")
                self.input.header.add_header_field(name=f, dictionary=d,
                                                   field_type='INFO')
            return sf
        else:
            if not self.prev_splice_ai:
                if self.args.splice_ai_min_delta is not None:
                    raise RuntimeError("--splice_ai_min_delta cutoff given " +
                                       "but VCF does not have SpliceAI INFO " +
                                       "field in its header and no SpliceAI " +
                                       "score files have been specified with " +
                                       "--splice_ai_vcfs argument.")
                if self.args.splice_ai_max_delta is not None:
                    raise RuntimeError("--splice_ai_max_delta cutoff given " +
                                       "but VCF does not have SpliceAI INFO " +
                                       "field in its header and no SpliceAI " +
                                       "score files have been specified with " +
                                       "--splice_ai_vcfs argument.")
        return None


    def get_vcf_filter_classes(self):
        filters = []
        uni_args = {}
        if self.args.freq is not None:
            uni_args["freq"] = self.args.freq
        if self.args.min_freq is not None:
            uni_args["min_freq"] = self.args.min_freq
        # get dbSNP filter
        for dbsnp in self.args.dbsnp:
            prefix = self.check_info_prefix('VASE_dbSNP')
            kwargs = {"vcf" : dbsnp, "prefix" : prefix,
                      "clinvar_path" : self.args.clinvar_path,}
            if self.args.build is not None:
                kwargs['build'] = self.args.build
            if self.args.max_build is not None:
                kwargs['max_build'] = self.args.max_build
            kwargs.update(uni_args)
            dbsnp_filter = dbSnpFilter(**kwargs)
            filters.append(dbsnp_filter)
            for f,d in dbsnp_filter.added_info.items():
                self.logger.debug("Adding dbSNP annotation {}" .format(f))
                self.input.header.add_header_field(name=f, dictionary=d,
                                          field_type='INFO')
        # get gnomAD/ExAC filters
        for gnomad in self.args.gnomad:
            prefix = self.check_info_prefix('VASE_gnomAD')
            kwargs = {"vcf" : gnomad,
                      "prefix" : prefix,
                      "pops": self.args.gnomad_pops,
                      "max_homozygotes": self.args.max_gnomad_homozygotes}
            kwargs.update(uni_args)
            gnomad_filter = GnomadFilter(**kwargs)
            filters.append(gnomad_filter)
            for f,d in gnomad_filter.added_info.items():
                self.logger.debug("Adding gnomAD/ExAC annotation {}"
                                  .format(f))
                self.input.header.add_header_field(name=f, dictionary=d,
                                          field_type='INFO')
        #get other VCF filters
        for var_filter in self.args.vcf_filter:
            vcf_and_id = var_filter.split(',')
            if len(vcf_and_id) < 2:
                raise RuntimeError("Expected two or more comma separated " +
                                   "values for --vcf_filter argument '{}'"
                                   .format(var_filter))
            prefix = self.check_info_prefix('VASE_' + vcf_and_id[1])
            kwargs = {"vcf" : vcf_and_id[0], "prefix" : prefix,
                      "annotations" : vcf_and_id[2:]}
            kwargs.update(uni_args)
            vcf_filter = VcfFilter(**kwargs)
            filters.append(vcf_filter)
            for f,d in vcf_filter.added_info.items():
                self.logger.debug("Adding {} annotation from {}"
                                  .format(f, vcf_and_id[0]))
                self.input.header.add_header_field(name=f, dictionary=d,
                                          field_type='INFO')
        return filters

    def get_gt_annotators(self):
        gt_annos = []
        if self.args.dng_vcf:
            for vcf in self.args.dng_vcf:
                g = (GtAnnotator(vcf, ['PP_DNM', 'PP_NULL']))
                for f,d in g.header_fields.items():
                    self.input.header.add_header_field(name=f, dictionary=d,
                                                       field_type='FORMAT')
                gt_annos.append(g)
        return gt_annos

    def _seg_alleles_from_record(self, record, prefix, filter_al, filter_csq):
        ps = prefix + '_samples'
        pf = prefix + '_features'
        info = record.parsed_info_fields(fields=[pf, ps])
        f_al = [False if info[ps][i] != '.' else filter_al[i] for i in
                range(len(filter_al))]
        f_csq = [[False if (info[pf][i] is not None and
                            record.CSQ[j]['Feature'] in info[pf][i].split('|'))
                  else filter_csq[i][j] for j in range(len(record.CSQ))] for i
                 in range(len(filter_al))]
        return f_al, f_csq

    def _get_prev_annotations(self):
        self.prev_annots = set()
        self.info_prefixes = set()
        for info in self.input.metadata['INFO']:
            match = re.search('^(VASE_\w+)_\w+(_\d+)?', info)
            if match:
                self.prev_annots.add(info)
                self.info_prefixes.add(match.group(1))
                self.logger.debug("Identified previously annotated VASE INFO" +
                                  " field '{}'" .format(info))
        self._parse_prev_vcf_filter_annotations()
        if 'CADD_PHRED_score' in self.input.metadata['INFO']:
            if (self.input.metadata['INFO']['CADD_PHRED_score'][-1]['Number'] == 'A' and
                self.input.metadata['INFO']['CADD_PHRED_score'][-1]['Type'] == 'Float'):
                self.prev_cadd_phred = True
        if 'CADD_raw_score' in self.input.metadata['INFO']:
            if (self.input.metadata['INFO']['CADD_raw_score'][-1]['Number'] == 'A' and
                self.input.metadata['INFO']['CADD_raw_score'][-1]['Type'] == 'Float'):
                self.prev_cadd_raw = True
        if 'SpliceAI' in self.input.metadata['INFO']:
            if self.input.metadata['INFO']['SpliceAI'][-1]['Number'] == '.':
                self.prev_splice_ai = True

    def _parse_prev_vcf_filter_annotations(self):
        frq_annots = []
        hom_annots = []
        bld_annots = []
        cln_annots = []
        get_matching = self.args.filter_known or self.args.filter_novel
        if not self.args.ignore_existing_annotations:
            if self.args.freq or self.args.min_freq or get_matching:
                for annot in sorted(self.prev_annots):
                    match = re.search('^VASE_dbSNP|gnomAD(_\d+)?_(CAF|AF)(_\w+)?',
                                      annot)
                    if match:
                        if (self.input.metadata['INFO'][annot][-1]['Number'] == 'A' and
                            self.input.metadata['INFO'][annot][-1]['Type'] == 'Float'):
                            self.logger.info("Found previous allele frequency " +
                                              "annotation '{}'".format(annot))
                            if len(match.groups()) == 3 and match.group(3) is not None:
                                pop = match.group(3).replace("_", "")
                                if pop.upper() not in [x.upper() for x in
                                                       self.args.gnomad_pops]:
                                    self.logger.info("Ignoring {} ".format(annot) +
                                                     "annotation as not in " +
                                                     "populations specified by " +
                                                     "--gnomad_pops")
                                    continue
                            frq_annots.append(annot)
            if self.args.max_gnomad_homozygotes is not None:
                for annot in sorted(self.prev_annots):
                    match =re.search(
                        '^VASE_gnomAD(_\d+)?_(Hom|Hemi|nhomalt)(_\w+)?',
                        annot)
                    if match:
                        if (self.input.metadata['INFO'][annot][-1]['Number'] == 'A' and
                            self.input.metadata['INFO'][annot][-1]['Type'] == 'Integer'):
                            self.logger.info("Found previous Hom/Hemi " +
                                              "annotation '{}'".format(annot))
                            if len(match.groups()) == 3 and match.group(3) is not None:
                                pop = match.group(3).replace("_", "")
                                if pop.upper() not in [x.upper() for x in
                                                       self.args.gnomad_pops]:
                                    self.logger.info("Ignoring {} ".format(annot) +
                                                     "annotation as not in " +
                                                     "populations specified by " +
                                                     "--gnomad_pops")
                                    continue
                            hom_annots.append(annot)
            if self.args.build or self.args.max_build or get_matching:
                for annot in sorted(self.prev_annots):
                    match = re.search('^VASE_dbSNP(_\d+)?_dbSNPBuildID', annot)
                    if match:
                        if (self.input.metadata['INFO'][annot][-1]['Number'] == 'A' and
                          self.input.metadata['INFO'][annot][-1]['Type'] == 'Integer'):
                            self.logger.info("Found previous dbSNP build " +
                                              "annotation '{}'".format(annot))
                            bld_annots.append(annot)
            if self.args.clinvar_path or get_matching:
                for annot in sorted(self.prev_annots):
                    match = re.search('^VASE_dbSNP(_\d+)?_CLNSIG', annot)
                    if match:
                        if (self.input.metadata['INFO'][annot][-1]['Number'] == 'A' and
                          self.input.metadata['INFO'][annot][-1]['Type'] == 'String'):
                            self.logger.info("Found previous ClinVar " +
                                              "annotation '{}'".format(annot))
                            cln_annots.append(annot)
            # if using gnomAD file for burden counting use internal annotations
            # for frequency filtering, if specified
            if self.args.gnomad_burden and (self.args.freq or
                                            self.args.min_freq):
                #check only relatively outbred pops by default
                pops = set(("POPMAX", "AFR", "AMR", "EAS", "FIN", "NFE", "SAS"))
                for info in self.input.metadata['INFO']:
                    match = re.search(r'''^AF_([A-Z]+)$''', info)
                    if match and match.group(1).upper() in pops:
                        if (self.input.metadata['INFO'][info][-1]['Number'] == 'A'
                                and
                                self.input.metadata['INFO'][info][-1]['Type'] ==
                                'Float'):
                            self.logger.info("Found gnomAD allele frequency " +
                                              "annotation '{}'".format(info))
                            frq_annots.append(info)
        self.prev_freqs = tuple(frq_annots)
        self.prev_homs = tuple(hom_annots)
        self.prev_builds = tuple(bld_annots)
        self.prev_clinvar = tuple(cln_annots)

    def check_info_prefix(self, name):
        if name in self.info_prefixes:
            self.logger.debug("INFO field {} already exists - trying another"
                              .format(name))
            match = re.search('_(\d+)$', name)
            if match:
                #already has an appended '_#' - increment and try again
                i = int(match.group(1))
                name = re.sub(match.group(1) + '$', str(i + 1), name)
                return self.check_info_prefix(name)
            else:
                #append _1
                return self.check_info_prefix(name + '_1')
        self.info_prefixes.add(name)
        return name

    def print_header(self):
        '''
            Write a VCF header for output that consists of the input VCF
            header data plus program arguments and any new INFO/FORMAT
            fields.
        '''

        vase_opts = []
        for k,v in vars(self.args).items():
            vase_opts.append('--{} {}'.format(k, v))
        self.input.header.add_header_field(name="vase",
                                   string='"' + str.join(" ", vase_opts) + '"')
        self.out.write(str(self.input.header))

    def get_output(self):
        '''
            Return an output filehandle. If no output specified return
            sys.stdout, else, if output name ends with .gz or .bgz return a
            bgzf.BgzfWriter object and otherwise return a standard
            filehandle.
        '''

        if isinstance(self.args.output, str):
            if self.args.output.endswith(('.gz', '.bgz')):
                try:
                    from Bio import bgzf
                except ImportError:
                    raise RuntimeError("Can not import bgzf via " +
                                       "biopython. Please install biopython " +
                                       "in order to write bgzip compressed " +
                                       "(.gz/.bgz) output.")
                fh = bgzf.BgzfWriter(self.args.output)
            else:
                fh = open(self.args.output, 'w')
        else:
            fh = sys.stdout
        return fh

    def get_report_filehandles(self):
        fhs = {'recessive': None,
               'dominant':  None,
               'de_novo':   None,}
        if self.args.report_prefix is not None:
            if self.args.biallelic or self.args.singleton_recessive:
                f = self.args.report_prefix + ".recessive.report.tsv"
                fhs['recessive'] = open(f, 'w')
            if self.args.dominant or self.args.singleton_dominant:
                f = self.args.report_prefix + ".dominant.report.tsv"
                fhs['dominant'] = open(f, 'w')
            if self.args.de_novo:
                f = self.args.report_prefix + ".de_novo.report.tsv"
                fhs['de_novo'] = open(f, 'w')
        return fhs

    def _set_to_true_if_true(self, alist, values):
        for i in range(len(alist)):
            if values[i]:
                alist[i] = True

    def _set_to_false_if_true(self, alist, values):
        for i in range(len(alist)):
            if values[i]:
                alist[i] = False

    def _get_family_filter(self):
        if self.family_filter is not None:
            return self.family_filter
        infer = True
        no_ped = False
        if not self.ped:
            if (not self.args.singleton_recessive and
                not  self.args.singleton_dominant):
                raise RuntimeError("Inheritance filtering options require a " +
                                   "PED file specified using --ped or else " +
                                   "samples specified using " +
                                   "--singleton_recessive or " +
                                   "--singleton_dominant arguments")
            else:
                self.ped = self._make_ped_io()
                no_ped = True
                infer = False
        if self.args.check_g2p_inheritance:
            g2p = self.g2p
        else:
            g2p = None
        self.family_filter = FamilyFilter(
                        ped=self.ped, vcf=self.input,
                        infer_inheritance=infer, g2p=g2p,
                        check_g2p_consequence=self.args.check_g2p_consequence,
                        logging_level=self.logger.level)
        for s in self.args.seg_controls:
            indv = Individual(s, s, 0, 0, 0, 1)
            try:
                self.ped.add_individual(indv)
            except PedError:
                pass
        if not no_ped:
            for s in set(self.args.singleton_recessive +
                         self.args.singleton_dominant):
                indv = Individual(s, s, 0, 0, 0, 2)
                try:
                    self.ped.add_individual(indv)
                except PedError:
                    raise RuntimeError("Sample '{}' ".format(s) + "specified" +
                                       " as either --singleton_recessive or " +
                                       "--singleton_dominant already exists " +
                                       "in PED file {}" .format(
                                                            self.ped.filename))
        for s in set(self.args.singleton_dominant):
            self.family_filter.inheritance_patterns[s].append('dominant')
        for s in set(self.args.singleton_recessive):
            self.family_filter.inheritance_patterns[s].append('recessive')

    def _parse_info_filters(self):
        if not self.args.info_filters:
            return None
        ifilters = []
        for expression in self.args.info_filters:
            exp = expression.split()
            if len(exp) != 3:
                raise ValueError("--info_filters must consist of three quoted"+ 
                                 " values separated by whitespace - for " +
                                 "example: 'QD > 4' The provided expression " +
                                 " '{}' is invalid.".format(expression))
            ifilters.append(exp)
        return InfoFilter(vcf=self.input, filters=ifilters)


    def _parse_var_type_arg(self):
        if not self.args.var_types:
            return None
        valid = {'SNV', 'INSERTION', 'DELETION', 'MNV', 'SV'}
        v_types = []
        for vt in self.args.var_types:
            if vt == 'INDEL':
                v_types.extend(['INSERTION', 'DELETION'])
            elif vt.upper() in valid:
                v_types.append(vt.upper())
            else:
                raise ValueError("Invalid --var_type provided: {}".format(vt))
        return v_types


    def _make_ped_io(self):
        p_string = ''
        for s in set(self.args.singleton_recessive +
                     self.args.singleton_dominant):
            p_string += str.join("\t", (s, s, "0", "0", "0", "2")) + "\n"
        ped = io.StringIO(p_string)
        return PedFile(ped)

    def _set_logger(self):
        self.logger = logging.getLogger("VASE")
        if self.args.silent:
            self.args.no_warnings = True
            self.args.no_progress = True
        if self.args.debug:
            self.logger.setLevel(logging.DEBUG)
        elif self.args.no_warnings:
            self.logger.setLevel(logging.ERROR)
        elif self.args.quiet:
            self.logger.setLevel(logging.WARNING)
        else:
            self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter(
                        '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setLevel(self.logger.level)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def _check_got_inherit_filter(self):
        if self.args.de_novo or self.args.biallelic or self.args.dominant:
            if (not self.recessive_filter and not self.de_novo_filter and
                not self.dominant_filter):
                raise RuntimeError("No inheritance filters could be created " +
                                   "with current settings. Please check your" +
                                   " ped/sample inputs or run without "+
                                   "--biallelic/--dominant/--de_novo options.")

    def _set_seg_annot_cleanup(self, seg_info):
        for_removal = [x for x in seg_info if x in self.prev_annots]
        if for_removal:
            self.info_to_remove = for_removal
            self.remove_previous_inheritance_filters = self._clean_info
            self.logger.warn("Replacing {:,} previous".format(len(for_removal))
                             + " segregation INFO fields.")
            for f in for_removal:
                self.logger.warn("Previous {} annotations will be removed"
                                 .format(f))
        else:
            self.remove_previous_inheritance_filters = lambda *args: None

    def _clean_info(self, record):
        ''' Remove predefined INFO fields from records. '''
        record.remove_info_fields(self.info_to_remove)

    def _get_dominant_filter(self):
        self._get_family_filter()
        self._get_control_filter()
        self.dominant_filter = DominantFilter(
                                       self.family_filter, self.gt_args,
                                       min_families=self.args.min_families,
                                       report_file=self.report_fhs['dominant'])
        added_info = list(self.dominant_filter.get_header_fields().keys())
        if not self.dominant_filter.affected:
            msg = ("No samples fit a dominant model - can not use dominant " +
                   "filtering")
            if not self.args.biallelic and not self.args.de_novo:
                raise RuntimeError("Error: " + msg)
            else:
                self.logger.warn(msg + ". Will continue with other models.")
                self.dominant_filter = None
        else:
            for f,d in self.dominant_filter.get_header_fields().items():
                self.logger.debug("Adding DominantFilter annotation {}"
                                  .format(f))
                self.input.header.add_header_field(name=f, dictionary=d,
                                                   field_type='INFO')
                if self.args.min_families > 1:
                    self.use_cache = True
        return added_info #so we know which fields to remove if necessary

    def _get_de_novo_filter(self):
        self._get_family_filter()
        self._get_control_filter()
        self.de_novo_filter = DeNovoFilter(
                                        self.family_filter, self.gt_args,
                                        min_families=self.args.min_families,
                                        report_file=self.report_fhs['de_novo'])
        added_info = list(self.de_novo_filter.get_header_fields().keys())
        if not self.de_novo_filter.affected:
            msg = ("No samples fit a de novo model - can not use de novo " +
                   "filtering")
            if not self.args.biallelic and not self.args.dominant:
                raise RuntimeError("Error: " + msg)
            else:
                self.logger.warn(msg + ". Will continue with other models.")
                self.de_novo_filter = None
        else:
            for f,d in self.de_novo_filter.get_header_fields().items():
                self.logger.debug("Adding DeNovoFilter annotation {}"
                                  .format(f))
                self.input.header.add_header_field(name=f, dictionary=d,
                                                   field_type='INFO')
                if self.args.min_families > 1:
                    self.use_cache = True
        return added_info #so we know which fields to remove if necessary

    def _get_recessive_filter(self):
        self._get_family_filter()
        self.recessive_filter = RecessiveFilter(
                              self.family_filter, self.gt_args,
                              min_families=self.args.min_families,
                              strict=self.strict_recessive_inheritance,
                              report_file=self.report_fhs['recessive'])
        added_info = list(self.recessive_filter.get_header_fields().keys())
        if not self.recessive_filter.affected:
            msg = ("No samples fit a recessive model - can not use biallelic" +
                  " filtering")
            if not self.args.de_novo and not self.args.dominant:
                raise RuntimeError("Error: " + msg)
            else:
                self.logger.warn(msg + ". Will continue with other models.")
                self.recessive_filter = None
        else:
            self.use_cache = True
            for f,d in self.recessive_filter.get_header_fields().items():
                self.logger.debug("Adding RecessiveFilter annotation {}"
                                  .format(f))
                self.input.header.add_header_field(name=f, dictionary=d,
                                                   field_type='INFO')
        return added_info #so we know which fields to remove if necessary

    def _get_control_filter(self):
        if self.control_filter:
            return
        self.control_filter = ControlFilter(vcf=self.input,
                                            family_filter=self.family_filter,
                                            n_controls=self.args.n_controls,
                                            gt_args=self.gt_args)

    def _var_or_vars(self, n):
        if n == 1:
            return 'variant'
        return 'variants'


class VariantCache(object):
    '''
        Store a collection of CachedVariants that can be outputted
        once certain conditions are met (e.g. once we've checked
        whether variants match compound heterozygous variation in a
        gene). Keeps track of VEP Features encountered while adding to
        the cache so that the cache can be released for output once the
        current record is outside of the relevant features.
    '''

    __slots__ = ['cache', 'features', 'output_ready']

    def __init__(self):
        self.cache = []
        self.features = set()
        self.output_ready = []

    def check_record(self, record):
        '''
            Check whether features in given record are disjoint with
            those in cache and if so move variants from cache to
            output_ready. The given record is NOT added to the cache.
        '''
        these_feats = set([x['Feature'] for x in record.CSQ])
        if self.features and these_feats.isdisjoint(self.features):
            self.add_cache_to_output_ready()
            self.features.clear()

    def add_record(self, record, can_output=False):
        these_feats = set([x['Feature'] for x in record.CSQ])
        if self.features and these_feats.isdisjoint(self.features):
            self.add_cache_to_output_ready()
            self.features = these_feats
        else:
            self.features.update(these_feats)
        self.cache.append(CachedVariant(record, can_output))

    def add_cache_to_output_ready(self):
        ''' Adds items in cache to output_ready and clears cache.'''
        self.output_ready.extend(self.cache)
        self.cache = []

class CachedVariant(object):
    '''
        Store a variant that will or might need outputting later
        depending on, for example, determining whether they fit a
        recessive inheritance pattern.
    '''

    __slots__ = ['record', 'can_output', 'var_id']

    def __init__(self, record, can_output=False):
        self.record = record
        self.can_output = can_output
        self.var_id = "{}:{}-{}/{}" .format(record.CHROM, record.POS,
                                            record.REF, record.ALT)
