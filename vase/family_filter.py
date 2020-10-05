from .sample_filter import SampleFilter, GtFilter
from .sv_gt_filter import SvGtFilter
import logging
from collections import OrderedDict, defaultdict


class FamilyFilter(object):
    '''
        Determine whether variants/alleles fit given inheritance
        patterns for families.

    '''

    def __init__(self, ped, vcf, infer_inheritance=True, g2p=None,
                 check_g2p_consequence=None, force_inheritance=None,
                 logging_level=logging.WARNING):
        '''
            Initialize with Family object from ped_file.py and a
            VcfReader object from vcf_reader.py. You may also specify an
            inheritance pattern (either 'recessive' or 'dominant'). If
            inheritance_pattern is not specified an attempt is made to
            infer an appropriate inheritance pattern based on the family
            structure and affecteds.

            Args:
                ped:    A PedFile object from ped_file.py. Must contain
                        at least one affected individual.

                vcf:    A VcfReader object containing data from at least
                        some of the affected individuals in the given
                        family.

                infer_inheritance:
                        If True, infer possible inheritance patterns
                        for each family in the PedFile. Inferred patterns
                        are stored in self.inheritance_patterns dict
                        (keys are families, values are lists of
                        inheritance patterns).

                g2p:    G2P object from vase.g2p for filtering on
                        presence and inheritance requirements from a G2P
                        file.

                check_g2p_consequence:
                        If using a G2P object for gene filtering, also
                        filter on consequence type as described for each
                        gene. Note that the mapping of mutation
                        consequence to consequence type is quite crude
                        and should be used with caution (see the
                        mutation_to_csq dict in vase/g2p.py for the
                        mappings used).

                force_inheritance:
                        Optionally specify an inheritance pattern to
                        test for each family - either 'dominant' or
                        'recessive' is allowed. If infer_inheritance is
                        True, these patterns will be tested in addition
                        to inferred patterns.

                logging_level:
                        The level at which logging messages are
                        displayed. Defaults to logging.WARNING

        '''
        self.logger = self._get_logger(logging_level)
        self.affected = tuple(ped.get_affected())
        self.unaffected = tuple(ped.get_unaffected())
        self.obligate_carriers = dict()
        self.ped = ped
        self.vcf = vcf
        self.g2p = g2p
        self.check_g2p_consequence = check_g2p_consequence
        if not self.affected:
            raise RuntimeError("No affected individuals found in PED file '{}'"
                               .format(ped.filename))
        self.vcf_affected = list(x for x in self.affected
                                 if x in self.vcf.header.samples)
        if not self.vcf_affected:
            raise RuntimeError("No affected individuals in PED file '{}'"
                               .format(ped.filename) + " found in VCF " +
                               "'{}'".format(vcf.filename))
        self.vcf_unaffected = list(x for x in self.unaffected
                                   if x in self.vcf.header.samples)
        self.vcf_samples = self.vcf_affected + self.vcf_unaffected
        self.inheritance_patterns = defaultdict(list)
        if infer_inheritance:
            self._infer_inheritance()
        if force_inheritance:
            if force_inheritance not in ('dominant', 'recessive'):
                raise RuntimeError("Unrecognised inheritance pattern " +
                                   "specified with 'force_inheritance' " +
                                   "argument. Valid options are 'dominant' " +
                                   "or 'recessive'.")
            for fid in self.ped.families:
                self.inheritance_patterns[fid].append(force_inheritance)

    def _infer_inheritance(self):
        '''
            Simplistic method for determining likely relevant
            inheritance pattern. For affected individuals in a family
            a check is made whether parents or grandparents are also
            affected. Currently only dominant or recessive inheritance
            is inferred, no attempt to infer X-linked or mitochondrial
            inheritance is made and it will not spot pseudodominance.

        '''
        for fid, fam in self.ped.families.items():
            n_affected = 0
            no_parents = True
            both_pars_unaffected = False
            dominant = False
            denovo = False
            recessive = False
            self.logger.info("Assessing inheritance pattern of family {}"
                             .format(fid))
            f_aff = tuple(fam.get_affected())
            obligate_carriers = set()
            if not f_aff:
                continue
            for iid in f_aff:
                self.logger.info("Checking affected individual {}".format(iid))
                n_affected += 1
                indv = fam.individuals[iid]
                if not indv.parents:
                    self.logger.info("No parents for affected individual {}"
                                     .format(iid))
                    continue
                no_parents = False
                p_unaff = 0
                for par in indv.parents:
                    # is parent affected
                    if par not in fam.individuals:
                        if par in self.vcf.header.samples:
                            self.logger.warn("Family '{}' parent '{}' ".format(
                                             fid, par) + "not specified in " +
                                             "PED, but present in VCF - " +
                                             "assuming unaffected")
                            self.vcf_samples.append(par)
                            self.vcf_unaffected.append(par)
                        p_unaff += 1
                        continue
                    parent = fam.individuals[par]
                    par_to_child = False
                    gpar_to_child = False
                    if parent.is_affected():
                        self.logger.info("Apparent vertical transmission " +
                                         "from {} -> {}" .format(par, iid))
                        par_to_child = True
                    else:
                        p_unaff += 1
                    for gpar in parent.parents:
                        if fam.individuals[gpar].is_affected():
                            gpar_to_child = True
                            msg = "Apparent vertical transmission "
                            if par_to_child:
                                msg += ("from {} -> {} -> {}"
                                        .format(gpar, par, iid))
                            else:
                                msg += ("with partial penetrance from " +
                                        "{} -> ({}) -> {}"
                                        .format(gpar, par, iid))
                                obligate_carriers.add(par)
                            self.logger.info(msg)
                    if par_to_child or gpar_to_child:
                        dominant = True
                if p_unaff == 2:
                    both_pars_unaffected = True
            if not dominant:
                recessive = True
            if no_parents or not both_pars_unaffected:
                # missing information on one/both parents - could be dominant
                dominant = True
            if recessive and n_affected == 1 and not no_parents:
                f_par = fam.individuals[f_aff[0]].parents
                if len(f_par) != 2:
                    self.logger.info("Can not analyze {} under ".format(fid) +
                                     "a de novo model due to missing parents" +
                                     " in ped")
                    dominant = True
                elif (f_par[0] not in self.vcf.header.samples or
                      f_par[1] not in self.vcf.header.samples):
                    self.logger.info("Can not analyze {} under ".format(fid) +
                                     "a de novo model due to missing parents" +
                                     " in VCF")
                else:
                    denovo = True
            elif recessive and n_affected > 1:
                # we can entertain apparent de novos due to somatic mosaicism
                # if all affecteds share a parent
                pars = fam.individuals[f_aff[0]].parents
                shared_pars = None
                if len(pars) != 2:
                    self.logger.info("Can not analyze {} under ".format(fid) +
                                     "a de novo model due to missing parents" +
                                     " in ped")
                    dominant = True
                else:
                    shared_pars = set(pars)
                    for i in range(1, len(f_aff)):
                        ipars = self.ped.individuals[f_aff[i]].parents
                        if ipars is None:
                            break
                        shared_pars = shared_pars.intersection(ipars)
                        if not shared_pars:
                            break
                if shared_pars:
                    denovo = True
                    for par in shared_pars:
                        if par not in self.vcf_samples:
                            self.logger.info("Can not analyze {}".format(fid) +
                                             "under a de novo model due to " +
                                             "missing parents in VCF")
                            denovo = False
                            break
            self.inheritance_patterns[fid] = []
            if recessive:
                self.logger.info("Family '{}' " .format(fid) + "can be " +
                                 "analysed under a recessive model")
                self.inheritance_patterns[fid].append('recessive')
            if denovo:
                dmodel = "de novo"
                if n_affected > 1:
                    dmodel += " (with germline mosaicism)"
                self.logger.info("Family '{}' " .format(fid) + "can be " +
                                 "analysed under a {} model" .format(dmodel))
                self.inheritance_patterns[fid].append('de_novo')
            if dominant:
                self.logger.info("Family '{}' " .format(fid) + "can be " +
                                 "analysed under a dominant model")
                self.inheritance_patterns[fid].append('dominant')
                self.obligate_carriers[fid] = tuple(obligate_carriers)

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


class InheritanceFilter(object):

    '''
        Parent class for RecessiveFilter/DominantFilter/DeNovoFilter
        object.
    '''

    def __init__(self, family_filter, gt_args, min_families=1,
                 report_file=None):
        '''
            Create genotype filter objects and initialise family filtering
            arguments.

            Args:
                family_filter:
                        Parent FamilyFilter object, initialized with
                        VCF and PED files.

                gt_args:
                        A dict of arguments to use for filtering
                        genotypes. These should all correspond to
                        arguments to provide to SampleFilter objects.

                min_families:
                        Require at least this many families to have
                        qualifying alleles in a feature before
                        outputting. Default=1.

                report_file:
                        Deprecated. Use vase_reporter to after
                        inheritance filtering to process VCFs instead.
        '''

        self.family_filter = family_filter
        self.min_families = min_families
        self.ped = family_filter.ped
        self.samples = family_filter.vcf_samples
        self.unaffected = family_filter.vcf_unaffected
        self.gt_filter = GtFilter(family_filter.vcf,
                                  gq=gt_args.get('gq'),
                                  dp=gt_args.get('dp'),
                                  max_dp=gt_args.get('max_dp'),
                                  het_ab=gt_args.get('het_ab'),
                                  hom_ab=gt_args.get('hom_ab'))
        self._gt_fields = set(self.gt_filter.fields)
        if gt_args.get('min_control_gq') is None:
            gt_args['min_control_gq'] = gt_args.get('gq')
        if gt_args.get('min_control_dp') is None:
            gt_args['min_control_dp'] = gt_args.get('dp')
        if gt_args.get('max_control_dp') is None:
            gt_args['max_control_dp'] = gt_args.get('max_dp')
        if gt_args.get('control_het_ab') is None:
            gt_args['control_het_ab'] = gt_args.get('het_ab')
        if gt_args.get('control_hom_ab') is None:
            gt_args['control_hom_ab'] = gt_args.get('hom_ab')
        self.con_gt_filter = GtFilter(family_filter.vcf,
                                      gq=gt_args.get('min_control_gq'),
                                      dp=gt_args.get('min_control_dp'),
                                      max_dp=gt_args.get('max_control_dp'),
                                      het_ab=gt_args.get('control_het_ab'),
                                      hom_ab=gt_args.get('control_hom_ab'),
                                      ref_ab_filter=gt_args.get('con_ref_ab'))
        self._gt_fields.update(self.con_gt_filter.fields)
        if gt_args.get('sv_min_control_gq') is None:
            gt_args['sv_min_control_gq'] = gt_args.get('sv_gq')
        if gt_args.get('sv_min_control_dp') is None:
            gt_args['sv_min_control_dp'] = gt_args.get('sv_dp')
        if gt_args.get('sv_max_control_dp') is None:
            gt_args['sv_max_control_dp'] = gt_args.get('sv_max_dp')
        if gt_args.get('sv_control_het_ab') is None:
            gt_args['sv_control_het_ab'] = gt_args.get('sv_het_ab')
        if gt_args.get('sv_control_hom_ab') is None:
            gt_args['sv_control_hom_ab'] = gt_args.get('sv_hom_ab')
        if gt_args.get('control_del_dhffc') is None:
            gt_args['control_del_dhffc'] = gt_args.get('del_dhffc')
        if gt_args.get('control_dup_dhbfc') is None:
            gt_args['control_dup_dhbfc'] = gt_args.get('dup_dhbfc')
        self.sv_gt_filter = SvGtFilter(family_filter.vcf,
                                       gq=gt_args.get('sv_gq'),
                                       dp=gt_args.get('sv_dp'),
                                       max_dp=gt_args.get('sv_max_dp'),
                                       het_ab=gt_args.get('sv_het_ab'),
                                       hom_ab=gt_args.get('sv_hom_ab'),
                                       del_dhffc=gt_args.get('del_dhffc'),
                                       dup_dhbfc=gt_args.get('dup_dhbfc'))
        self._sv_gt_fields = set(self.sv_gt_filter.fields)
        self.sv_con_gt_filter = SvGtFilter(
                                    family_filter.vcf,
                                    gq=gt_args.get('sv_min_control_gq'),
                                    dp=gt_args.get('sv_min_control_dp'),
                                    max_dp=gt_args.get('sv_max_control_dp'),
                                    het_ab=gt_args.get('sv_control_het_ab'),
                                    hom_ab=gt_args.get('sv_control_hom_ab'),
                                    ref_ab_filter=gt_args.get('sv_con_ref_ab'),
                                    del_dhffc=gt_args.get('control_del_dhffc'),
                                    dup_dhbfc=gt_args.get('control_dup_dhbfc'))
        self._sv_gt_fields.update(self.sv_con_gt_filter.fields)
        self._prev_coordinate = (None, None)  # to ensure records are processed
        self._processed_contigs = set()       # in coordinate order
        if self.report_file:
            self._write_report_header()

    def get_header_fields(self):
        '''
            Return dict of dicts with INFO header field names as keys
            and dicts of features as values. These are suitable for
            handing to VcfHeader class's add_header_field() method.

            Each INFO field must be defined in self.header_fields in
            the child class, which should be a list of tuples where
            each tuple consists of the name anddescription of the
            field.
        '''
        hf = dict()
        for f in self.header_fields:
            hf[f[0]] = {'Number': 'A',
                        'Type': 'String',
                        'Description': f[1]}
        return hf

    def confirm_heterozygous(self, record, samples):
        for s in samples:
            if len(set(record.samples[s]['GT'])) != 2:
                return False
        return True

    def _get_allele_counts(self, allele, record):
        a_counts = dict()
        gt_filter_args = dict()
        if record.IS_SV:
            gt_filter = self.sv_gt_filter
            control_filter = self.sv_con_gt_filter
            gt_filter_args['svtype'] = record.info.get('SVTYPE', '')
        else:
            gt_filter = self.gt_filter
            control_filter = self.con_gt_filter
        for samp in self.unaffected:
            if control_filter.gt_is_ok(record.samples, samp, allele,
                                       **gt_filter_args):
                a_counts[samp] = record.samples[samp]['GT'].count(allele)
            else:
                a_counts[samp] = None
            if (record.samples[samp]['GT'] == (0, 0) and
                    control_filter.ad_over_threshold is not None):
                if control_filter.ad_over_threshold(record.samples, samp,
                                                    allele):
                    a_counts[samp] = 1
        for samp in self.affected:
            if gt_filter.gt_is_ok(record.samples, samp, allele,
                                  **gt_filter_args):
                a_counts[samp] = record.samples[samp]['GT'].count(allele)
            else:
                a_counts[samp] = None
        return a_counts

    def _check_sorted(self, record):
        if self._prev_coordinate[0] != record.chrom:
            if record.chrom in self._processed_contigs:
                raise RuntimeError("Input must be sorted by chromosome and " +
                                   "position for recessive filtering. " +
                                   "Contig '{}' " .format(record.chrom) +
                                   "encountered before and after contig " +
                                   "'{}'." .format(self._prev_coordinate[0]))
            if self._prev_coordinate[0] is not None:
                self._processed_contigs.add(self._prev_coordinate[0])
        elif record.pos < self._prev_coordinate[1]:
            raise RuntimeError("Input must be sorted by chromosome and " +
                               "position for inheritance filtering. " +
                               "Encountered position {}:{} after {}:{}"
                               .format(record.chrom, record.pos,
                                       self._prev_coordinate[0],
                                       self._prev_coordinate[1]))
        self._prev_coordinate = (record.chrom, record.pos)

    def process_record(self, record):
        '''Return True if record should be printed/kept'''
        return NotImplementedError("process_record method should be " +
                                   "overriden by child class!")

    def _write_report_header(self):
        header = str.join("\t", (x for x in
                                 self.family_filter.vcf.header.csq_fields
                                 if x != 'Allele'))
        header += "\tALT_No.\t" + str.join("\t", self.annot_fields)
        header += "\tCHROM\tPOS\tID\tREF\tALT\tALLELE\tQUAL\tFILTER"
        self.report_file.write(header + "\n")

    def check_g2p(self, record, ignore_csq, inheritance):
        if self.family_filter.g2p:
            if self.family_filter.check_g2p_consequence:
                fail = (not x for x in
                        self.family_filter.g2p.csq_and_allelic_requirement_met(
                            record, inheritance))
            else:
                fail = (not x for x in
                        self.family_filter.g2p.allelic_requirement_met(
                                     record, inheritance))
            if ignore_csq:
                ignore_csq = [x or y for x, y in zip(ignore_csq, fail)]
            else:
                ignore_csq = list(fail)
        return ignore_csq


class RecessiveFilter(InheritanceFilter):
    '''
        This class assumes that each family has a shared biallelic
        genetic cause of disease. It will not cope with phenocopies,
        pseudodominance or other more complicated inheritance patterns.
    '''
    def __init__(self, family_filter, gt_args, min_families=1, strict=False,
                 exclude_denovo=False, report_file=None):
        '''
            Args:
                family_filter:
                        FamilyFilter object

                gt_args:
                        A dict of arguments to use for filtering
                        genotypes. These should all correspond to
                        arguments to provide to SampleFilter objects.

                strict: If True, for any affected sample with
                        parents, require confirmation of parental
                        genotypes. If either parent genotype is a
                        no-call for a record, then the record will
                        be ignored. Default=False.

                exclude_denovo:
                        If True, where there is data available from
                        both parents for an affected individual
                        ignore apparent de novo occuring alleles.
                        Default=False.

                min_families:
                        Require at least this many families to have a
                        qualifying biallelic combination of alleles in
                        a feature before outputting. Default=1.

                report_file:
                        Output filehandle for writing summaries of
                        segregating variants to. Default=None.

        '''
        self.prefix = "VASE_biallelic"
        self.header_fields = [
              ("VASE_biallelic_homozygous",
               'Samples that carry homozygous biallelic changes ' +
               ' parsed by {}' .format(type(self).__name__)),
              ("VASE_biallelic_compound_het",
               'Samples that carry compound heterozygous biallelic changes ' +
               'parsed by {}'.format(type(self).__name__)),
              ("VASE_biallelic_de_novo",
               'Samples that carry biallelic alleles that appear to have ' +
               'arisen de novo'),
              ('VASE_biallelic_families',
               'Family IDs for VASE_biallelic alleles'),
              ("VASE_biallelic_features",
               'Features (e.g. transcripts) that contain qualifying ' +
               'biallelic variants parsed by {}' .format(
                   type(self).__name__))]
        self.annot_fields = ('homozygous', 'compound_het', 'de_novo',
                             'families', 'features')
        self.report_file = report_file
        super().__init__(family_filter, gt_args, min_families=min_families,
                         report_file=report_file,)
        self.families = tuple(x for x in
                              self.family_filter.inheritance_patterns
                              if 'recessive' in
                              self.family_filter.inheritance_patterns[x])
        self.affected = tuple(x for x in family_filter.vcf_affected if
                              self.ped.individuals[x].fid in self.families)
        self._fam_to_aff = dict()
        for fid in self.families:
            self._fam_to_aff[fid] = set(x for x in
                                        self.ped.families[fid].get_affected()
                                        if x in self.affected)
            self.family_filter.logger.info("Analysing family {} ".format(fid) +
                                           "under a recessive model")
        self.strict = strict
        self.exclude_denovo = exclude_denovo
        self._potential_recessives = dict()
        self._current_features = set()
        self._processed_features = set()

    def process_record(self, record, ignore_alleles=[], ignore_csq=[]):
        '''
            Returns True if record should be stored for checking against
            other records overlapping the same features to see if they
            constitute biallelic variation.

            Stores potential recessive records per allele for
            segregation checking once overlapping features have been
            traversed.

            Args:
                record: VaseRecord

                ignore_alleles:
                        List of booleans indicating for each ALT in
                        order whether it should be ignored in relation
                        to possible recessive variation (e.g. if MAF is
                        too high, no likely pathogenic consequence
                        etc.). This will normally have been generated
                        by VaseRunner via VcfFilter and/or VepFilter
                        classes.

                ignore_csq:
                        List of booleans indicating for each CSQ in
                        order whether it should be ignored in relation
                        to possible recessive variation. This should
                        normally have been generated by a corresponding
                        VepFilter object.

        '''
        stored = False
        self._check_sorted(record)
        self._current_features = set(c['Feature'] for c in record.CSQ if
                                     c['Feature'] != '')
        ignore_csq = self.check_g2p(record, ignore_csq, 'recessive')
        if ignore_csq and all(ignore_csq):
            return False
        gt_filter_args = dict()
        if record.IS_SV:
            gt_filter = self.sv_gt_filter
            control_filter = self.sv_con_gt_filter
            gt_filter_args['svtype'] = record.info.get('SVTYPE', '')
        else:
            gt_filter = self.gt_filter
            control_filter = self.con_gt_filter
        skip_fam = set()
        added_prs = OrderedDict()
        for i in range(len(record.alts)):
            if ignore_alleles and ignore_alleles[i]:
                continue
            alt = i + 1
            skip_allele = False
            fams_with_allele = []
            for un in self.unaffected:
                if record.samples[un]['GT'] == (alt, alt):
                    if control_filter.gt_is_ok(record.samples, un, alt,
                                               **gt_filter_args):
                        # hom in a control - skip allele
                        skip_allele = True
                        break
            if skip_allele:
                continue
            for fid in self.families:
                if fid in skip_fam:
                    continue
                have_allele = set()  # affecteds carrying this allele
                for aff in self._fam_to_aff[fid]:
                    # check all affecteds carry this allele
                    if (alt in record.samples[aff]['GT'] and
                            gt_filter.gt_is_ok(record.samples, aff, alt,
                                               **gt_filter_args)):
                        have_allele.add(aff)
                    else:
                        break
                if have_allele == self._fam_to_aff[fid]:
                    # all affecteds in family carry allele
                    fams_with_allele.append(fid)
            if fams_with_allele:
                # store record and consequences
                try:
                    csqs = []
                    for j in range(len(record.CSQ)):
                        if ignore_csq and ignore_csq[j]:
                            continue
                        if record.CSQ[j]['alt_index'] == alt:
                            # store record and csq details
                            csqs.append(record.CSQ[j])
                    if csqs:
                        stored = True
                        alt_counts = self._get_allele_counts(alt, record)
                        pr = PotentialSegregant(record=record, allele=alt,
                                                csqs=csqs,
                                                allele_counts=alt_counts,
                                                families=fams_with_allele)
                        for feat in pr.features:
                            if feat in added_prs:
                                added_prs[feat][pr.alt_id] = pr
                            else:
                                added_prs[feat] = OrderedDict(
                                    [(pr.alt_id, pr)])
                            if feat in self._potential_recessives:
                                self._potential_recessives[feat][pr.alt_id] = pr
                            else:
                                self._potential_recessives[feat] = OrderedDict(
                                    [(pr.alt_id, pr)])
                except KeyError:
                    raise RuntimeError("Could not identify CSQ or ANN " +
                                       "fields in VCF header. Please ensure " +
                                       "your input is annotated with " +
                                       "Ensembl's VEP to perform recessive " +
                                       "filtering")
        return stored

    def process_potential_recessives(self, final=False):
        '''
            Check whether stored PotentialSegregant alleles make up
            biallelic variation in the same transcript for affected
            individuals/families. Adds labels to INFO fields of VCF
            records and returns an OrderedDict of 'var_ids' to
            lists of PotentialSegregant objects that appear to
            segregate consistent with recessive inheritance.

            Clears the cache of stored PotentialSegregant alleles.
        '''
        segregating = OrderedDict()  # key=alt_id, val=SegregatingBiallelic
        for feat, prs in self._potential_recessives.items():
            if not final and feat in self._current_features:
                continue
            feat_segregating = []  # list of tuples of values for creating SegregatingBiallelic
            un_hets = defaultdict(list)  # store het alleles carried by each unaffected
            aff_hets = defaultdict(list)  # store het alleles carried by each affected
            biallelics = defaultdict(list)  # store biallelic combinations for affecteds
            for pid, p in prs.items():
                for un in self.unaffected:
                    if p.allele_counts[un] == 1:  # already checked for homs when adding
                        # store allele carried in this unaffected
                        un_hets[un].append(pid)
                for aff in (x for x in self.affected
                            if self.ped.fid_from_iid(x) in p.families):
                    if p.allele_counts[aff] == 1:
                        aff_hets[aff].append(pid)
                    elif p.allele_counts[aff] == 2:
                        biallelics[aff].append(tuple([pid]))
            incompatibles = []  # create a list of sets of incompatible hets
            for hets in un_hets.values():
                if len(hets):
                    incompatibles.append(set(hets))
            for aff, hets in aff_hets.items():
                for i in range(len(hets)):
                    for j in range(i+1, len(hets)):
                        incomp = False
                        for iset in incompatibles:
                            if iset.issuperset([hets[i], hets[j]]):
                                incomp = True
                                break
                        if not incomp:
                            if not prs[hets[i]].record.in_cis_with(sample=aff,
                                allele=prs[hets[i]].allele,
                                other=prs[hets[j]].record,
                                other_allele=prs[hets[j]].allele):
                                # check phase groups in case alleles in cis
                                biallelics[aff].append(
                                                     tuple([hets[i], hets[j]]))
            if not biallelics:
                continue
            # see if all affecteds in the same family share the same biallelics
            for fid, affs in self._fam_to_aff.items():
                b_affs = set(x for x in affs if x in biallelics)
                if len(b_affs) == 0 or b_affs != affs:
                    continue
                affs = list(affs)
                absent_in_aff = False
                for i in range(len(affs)):
                    for bi in biallelics[affs[i]]:
                        for j in range(i+1, len(affs)):
                            if bi not in biallelics[affs[j]]:
                                absent_in_aff = True
                                break
                        if not absent_in_aff:
                            segs, de_novo = self._check_parents(feat, bi, affs)
                            if not segs:
                                continue
                            if len(bi) == 1:
                                model = 'homozygous'
                            else:
                                model = 'compound_het'
                            for bi_pr in (prs[x] for x in bi):
                                feat_segregating.append((bi_pr, affs, [fid],
                                                         model, [feat],
                                                         de_novo[bi_pr.alt_id],
                                                         self.prefix))
            fam_count = len(set([fam for tup in feat_segregating for fam in
                                 tup[2]]))
            if fam_count >= self.min_families:
                for tp in feat_segregating:
                    if tp[0] in segregating:
                        segregating[tp[0]].add_samples(*tp[1:6])
                    else:
                        segregating[tp[0]] = SegregatingVariant(*tp)
        var_to_segregants = OrderedDict()
        for sb in segregating.values():
            sb.annotate_record(self.report_file, self.annot_fields)
            if sb.segregant.var_id in var_to_segregants:
                var_to_segregants[sb.segregant.var_id].append(sb.segregant)
            else:
                var_to_segregants[sb.segregant.var_id] = [sb.segregant]
        # clear the cache except for the last entry which will be a new gene
        # self._potential_recessives = self._last_added
        self._potential_recessives = OrderedDict(
            (k, v) for k, v in self._potential_recessives.items() if k in
            self._current_features)
        return var_to_segregants

    def _check_parents(self, feat, alleles, samples):
        '''
            Check transmission of alleles (i.e. one from each parent)
            if parents available. Should have already checked that
            alleles are not present in this combination in any
            unaffected individual.

            Returns a tuple of booleans - first value is True if
            parental genotypes do not contradict recessive inheritance
            while the second value is a dict of alleles to lists of
            samples in which the allele allele appears to have arisen
            de novo.
        '''
        dns = defaultdict(list)
        counts = []
        for al in alleles:
            counts.append(self._potential_recessives[feat][al].allele_counts)
        if len(counts) == 1:  # homozygous
            counts.append(counts[0])
        for samp in samples:
            parents = self.ped.individuals[samp].parents
            par = list(x for x in parents if x in self.samples)
            if len(par) == 0:
                continue
            if self.strict:
                for p in par:
                    if None in (counts[i][p] for i in range(len(counts))):
                        # require both parental genotypes if self.strict
                        return (False, dns)
            if len(par) == 2:  # can check for de novos
                for i in range(len(counts)):
                    if counts[i][par[0]] == 0 and counts[i][par[1]] == 0:
                        # apparent de novo
                        self.family_filter.logger.debug(
                            "Apparent de novo allele " +
                            "{} for sample {} (parents = {} + {}) ".format(
                                alleles[-i], samp, par[0], par[1]) +
                            "for recessive combination {}|{}".format(
                                alleles[0], alleles[-1]))
                        dns[alleles[-i]].append(samp)
                        if self.exclude_denovo:
                            return (False, dns)
            elif len(par) == 1:
                # if only one parent and both alleles are absent it is more
                # likely that the two alleles are in cis from other parent
                if counts[0][par[0]] == 0 and counts[1][par[0]] == 0:
                    return(False, dns)
            # NOTE: we could do a check here to ensure that any non-affected
            #       parent does not carry both alleles, but this *SHOULD* have
            #       already been done earlier in process_potential_recessives
            #       function for ALL unaffecteds anyway
        return (True, dns)


class DominantFilter(InheritanceFilter):
    '''
        Identify variants that fit a dominant pattern in
        given families.
    '''

    def __init__(self, family_filter, gt_args, min_families=1,
                 report_file=None):
        '''
            Initialize with parent IDs, children IDs and VcfReader
            object.

            Args:
                family_filter:
                        FamilyFilter object

                gt_args:
                        A dict of arguments to use for filtering
                        genotypes. These should all correspond to
                        arguments to provide to SampleFilter objects.

                min_families:
                        Require at least this many families to have a
                        qualifying variant in a feature before
                        outputting. Default=1.

        '''
        self.prefix = "VASE_dominant"
        self.header_fields = [
            ("VASE_dominant_samples",
             'Sample IDs for alleles that segregate according to a ' +
             'dominant inheritance pattern in an affected sample as' +
             ' parsed by {}' .format(type(self).__name__)),
            ('VASE_dominant_unaffected_carrier',
             'Sample IDs for unaffected carriers of ' +
             'VASE_dominant alleles'),
            ('VASE_dominant_families',
             'Family IDs for VASE_dominant alleles'),
            ("VASE_dominant_features",
             'Features (e.g. transcripts) that contain qualifying ' +
             'dominant variants parsed by {}' .format(
                 type(self).__name__))]
        self.annot_fields = ('samples', 'unaffected_carrier', 'families',
                             'features')
        self.report_file = report_file
        super().__init__(family_filter, gt_args, min_families=min_families,
                         report_file=report_file,)
        self.families = tuple(x for x in
                              self.family_filter.inheritance_patterns
                              if 'dominant' in
                              self.family_filter.inheritance_patterns[x])
        self.affected = tuple(x for x in family_filter.vcf_affected if
                              self.ped.individuals[x].fid in self.families)
        self.filters = dict()
        self._potential_dominants = dict()
        self._last_added = OrderedDict()
        self._current_features = set()
        for fam in self.families:
            f_aff = tuple(x for x in self.ped.families[fam].get_affected()
                          if (x in self.affected or
                          x in self.family_filter.obligate_carriers[fam]))
            f_unaff = tuple(x for x in self.ped.families[fam].get_unaffected()
                            if (x in self.unaffected and x not in
                            self.family_filter.obligate_carriers[fam]))
            if fam in self.family_filter.obligate_carriers:
                self.obligate_carriers = tuple(
                    x for x in f_aff if x in
                    self.family_filter.obligate_carriers[fam])
            else:
                self.obligate_carriers = ()
            dom_filter = SampleFilter(family_filter.vcf, cases=f_aff,
                                      controls=f_unaff, confirm_missing=True,
                                      **gt_args)
            self.filters[fam] = dom_filter
            self.family_filter.logger.info("Analysing family {} ".format(fam) +
                                           "under a dominant model")

    def process_record(self, record, ignore_alleles=[], ignore_csq=[]):
        '''
            Returns True if an allele segregates consistent with
            dominant inheritance.

            Args:
                record: VaseRecord

                ignore_alleles:
                        List of booleans indicating for each ALT in
                        order whether it should be ignored in relation
                        to possible dominant variation (e.g. if MAF is
                        too high, no likely pathogenic consequence
                        etc.). This will normally have been generated
                        by VaseRunner via VcfFilter and/or VepFilter
                        classes.

        '''
        dom_alleles = ([[] for i in range(len(record.alts))])
        fam_alleles = ([[] for i in range(len(record.alts))])
        ignore_csq = self.check_g2p(record, ignore_csq, 'dominant')
        if ignore_csq and all(ignore_csq):
            return False
        if self.min_families > 1:
            self._check_sorted(record)
        for i in range(len(record.alts)):
            if ignore_alleles[i]:
                continue
            allele = i + 1
            for fam, dfilter in self.filters.items():
                # looking for (potentially shared) de novos in a single family
                is_dom = not dfilter.filter(record, allele)
                if is_dom:
                    if self.confirm_heterozygous(record, dfilter.cases):
                        dom_alleles[i].extend(dfilter.cases)
                        fam_alleles[i].append(fam)
                        self.family_filter.logger.debug(
                            "Apparent dominant allele {}:{}-{}/{} ".format(
                                record.chrom, record.pos, record.ref,
                                record.alleles[allele]) +
                            "present in {} ".format(dfilter.cases) +
                            "and absent in {}".format(dfilter.controls))
        segs = []
        for i in range(len(dom_alleles)):
            if not dom_alleles[i]:
                continue
            allele = i + 1
            csqs = []
            try:
                for j in range(len(record.CSQ)):
                    if ignore_csq and ignore_csq[j]:
                        continue
                    if record.CSQ[j]['alt_index'] == allele:
                        # store record and csq details
                        csqs.append(record.CSQ[j])
            except KeyError:
                if self.min_families > 1:
                    raise RuntimeError("Could not identify CSQ or ANN fields" +
                                       " in VCF header. Please ensure your " +
                                       "input is annotated with Ensembl's " +
                                       "VEP to perform dominant filtering.")
            if self.min_families <= 1 or csqs:
                a_counts = self._get_allele_counts(allele, record)
                pd = PotentialSegregant(record=record, allele=allele,
                                        csqs=csqs, allele_counts=a_counts,
                                        families=fam_alleles[i])
                segs.append(pd)
        if self.min_families > 1:
            for feat, od in self._last_added.items():
                if feat in self._potential_dominants:
                    self._potential_dominants[feat].update(od)
                else:
                    self._potential_dominants[feat] = od
            self._last_added = OrderedDict()
            for seg in segs:
                for feat in seg.features:
                    self._last_added[feat] = OrderedDict([(seg.alt_id, seg)])
        else:
            for seg in segs:
                affs = (x for x in self.affected
                        if x not in self.obligate_carriers and
                        self.ped.fid_from_iid(x) in seg.families)
                sv = SegregatingVariant(seg, affs, seg.families, 'samples',
                                        seg.features, [], self.prefix)
                obcs = tuple(x for x in self.obligate_carriers if
                             self.ped.fid_from_iid(x) in seg.families)
                if obcs:
                    obfs = set(self.ped.fid_from_iid(x) for x in obcs)
                    sv.add_samples(obcs, obfs, 'unaffected_carrier',
                                   seg.features, [])
                sv.annotate_record(self.report_file, self.annot_fields)
        return len(segs) > 0

    def process_dominants(self, final=False):
        '''
            Check whether stored PotentialSegregant alleles make up
            dominant variation in the same transcript for the minimum
            number of families. Adds labels to INFO fields of VCF
            records and returns an OrderedDict of 'var_ids' to
            lists of PotentialSegregant objects that appear to
            constitute dominant variation.

            Clears the cache of stored PotentialSegregant alleles.
        '''
        sds = OrderedDict()
        feat_processed = []
        if not self._potential_dominants:
            # if cache is empy, we never encountered the next set of features
            self._potential_dominants = self._last_added
            self._last_added = OrderedDict()
        elif final:
            for feat in self._last_added:
                if feat in self._potential_dominants:
                    self._potential_dominants[feat].update(
                                                        self._last_added[feat])
                else:
                    self._potential_dominants[feat] = self._last_added[feat]
            self._last_added = OrderedDict()
        for feat, pds in self._potential_dominants.items():
            if feat in self._current_features:  # still processing this feature
                continue
            feat_fams = set()
            feat_processed.append(feat)
            for pid, p in pds.items():
                feat_fams.update(p.families)
            if len(feat_fams) >= self.min_families:
                for p in pds.values():
                    samps = (x for x in self.affected
                             if self.ped.fid_from_iid(x) in p.families)
                    if p.alt_id in sds:
                        sds[p.alt_id].add_samples(samps, p.families,
                                                  'samples', [feat], [])
                    else:
                        sv = SegregatingVariant(p, samps, p.families,
                                                'samples', [feat], [],
                                                self.prefix)
                        sds[p.alt_id] = sv
        var_to_segregants = OrderedDict()
        for sv in sds.values():
            sv.annotate_record(self.report_file, self.annot_fields)
            if sv.segregant.var_id in var_to_segregants:
                var_to_segregants[sv.segregant.var_id].append(sv.segregant)
            else:
                var_to_segregants[sv.segregant.var_id] = [sv.segregant]
        # clear the cache of processed features
        for feat in feat_processed:
            del self._potential_dominants[feat]
        return var_to_segregants


class DeNovoFilter(InheritanceFilter):
    '''
        Identify and output variants occuring in a child and absent from
        the parents.
    '''

    def __init__(self, family_filter, gt_args, min_families=1, confirm_het=False,
                 report_file=None):
        '''
            Initialize with parent IDs, children IDs and VcfReader
            object.

            Args:
                family_filter:
                        FamilyFilter object

                gt_args:
                        A dict of arguments to use for filtering
                        genotypes. These should all correspond to
                        arguments to provide to SampleFilter objects.

                min_families:
                        Require at least this many families to have a
                        qualifying variant in a feature before
                        outputting. Default=1.

                confirm_het:
                        If True, apparent de novos are required to be
                        called as heterozygous. Default=False.

        '''
        self.prefix = "VASE_de_novo"
        self.header_fields = [("VASE_de_novo_samples",
                   'Samples that carry alleles occurring de novo parsed by ' +
                   '{}' .format(type(self).__name__)),
                    ('VASE_de_novo_families',
                    'Family IDs for VASE_de_novo alleles'),
                    ("VASE_de_novo_features",
                    'Features (e.g. transcripts) that contain qualifying ' +
                    'de novo variants parsed by {}' .format(
                    type(self).__name__)),]
        self.annot_fields = ('samples', 'families', 'features')
        self.report_file = report_file
        super().__init__(family_filter, gt_args, min_families=min_families,
                         report_file=report_file,)
        self.families = tuple(x for x in self.family_filter.inheritance_patterns
                             if 'de_novo' in
                             self.family_filter.inheritance_patterns[x])
        self.affected = tuple(x for x in family_filter.vcf_affected if
                             self.ped.individuals[x].fid in self.families)
        self._potential_denovos = dict()
        self._last_added = OrderedDict()
        self._current_features = set()
        self.confirm_het = confirm_het
        self.filters = defaultdict(list)
        self.prefix = "VASE_de_novo"
        for fam in self.families:
            f_aff = tuple(x for x in self.ped.families[fam].get_affected()
                          if x in self.affected)
            par_child_combos = defaultdict(list)
            for aff in f_aff:
                pars = tuple(x for x in
                             self.ped.families[fam].individuals[aff].parents
                             if x in self.samples)
                if len(pars) == 2:
                    par_child_combos[pars].append(aff)
            for parents, children in par_child_combos.items():
                par_filter = SampleFilter(family_filter.vcf, cases=children,
                                          controls=parents,
                                          confirm_missing=True, **gt_args)
                self.filters[fam].append(par_filter)
                self.family_filter.logger.info(
                    "Analysing family {} parents ({}) and children ({})"
                    .format(fam, str.join(", ", parents),
                            str.join(", ", children)) +
                    " combinations under a de novo dominant model")

    def process_record(self, record, ignore_alleles=[], ignore_csq=[]):
        '''
            Returns True if allele is an apparent de novo variant.

            Args:
                record: VaseRecord

                ignore_alleles:
                        List of booleans indicating for each ALT in
                        order whether it should be ignored in relation
                        to possible de novo  variation (e.g. if MAF is
                        too high, no likely pathogenic consequence
                        etc.). This will normally have been generated
                        by VaseRunner via VcfFilter and/or VepFilter
                        classes.

        '''
        if self.min_families > 1:
            self._check_sorted(record)
        ignore_csq = self.check_g2p(record, ignore_csq, 'de novo')
        if ignore_csq and all(ignore_csq):
            return False
        denovo_alleles = ([[] for i in range(len(record.alts))])
        fam_alleles = ([[] for i in range(len(record.alts))])
        for i in range(len(record.alts)):
            if ignore_alleles[i]:
                continue
            allele = i + 1
            for fam, filters in self.filters.items():
                # looking for (potentially shared) de novos in a single family
                dns = []
                for dfilter in filters:
                    is_denovo = not dfilter.filter(record, allele)
                    if is_denovo:
                        if (not self.confirm_het or self.confirm_heterozygous(
                                record, dfilter.cases)):
                            dns.append(dfilter.cases)
                            self.family_filter.logger.debug("Apparent de " +
                                "novo allele {}:{}-{}/{} ".format(
                                record.chrom, record.pos, record.ref,
                                record.alleles[allele]) + "present in {} "
                                .format(dfilter.cases) + "and absent in {}"
                                .format(dfilter.controls))
                if len(dns) == len(filters):  # all affecteds in fam have dnm
                    ([denovo_alleles[i].extend(x) for x in dns])
                    fam_alleles[i].append(fam)
        segs = []
        for i in range(len(denovo_alleles)):
            if not denovo_alleles[i]:
                continue
            allele = i + 1
            csqs = []
            try:
                for j in range(len(record.CSQ)):
                    if ignore_csq and ignore_csq[j]:
                        continue
                    if record.CSQ[j]['alt_index'] == allele:
                        # store record and csq details
                        csqs.append(record.CSQ[j])
            except KeyError:
                if self.min_families > 1:
                    raise RuntimeError("Could not identify CSQ or ANN fields" +
                                       " in VCF header. Please ensure your " +
                                       "input is annotated with Ensembl's " +
                                       "VEP to perform de novo filtering.")
            if self.min_families <= 1 or csqs:
                a_counts = self._get_allele_counts(allele, record)
                pd = PotentialSegregant(record=record, allele=allele,
                                        csqs=csqs, allele_counts=a_counts,
                                        families=fam_alleles[i])
                segs.append(pd)
        if self.min_families > 1:
            for feat, od in self._last_added.items():
                if feat in self._potential_denovos:
                    self._potential_denovos[feat].update(od)
                else:
                    self._potential_denovos[feat] = od
            self._last_added = OrderedDict()
            for seg in segs:
                for feat in seg.features:
                    self._last_added[feat] = OrderedDict([(seg.alt_id, seg)])
        else:
            for seg in segs:
                affs = (x for x in self.affected if self.ped.fid_from_iid(x)
                        in seg.families)
                sv = SegregatingVariant(seg, affs, seg.families, 'samples',
                                        seg.features, [], self.prefix)
                sv.annotate_record(self.report_file, self.annot_fields)
        return len(segs) > 0

    def process_de_novos(self, final=False):
        '''
            Check whether stored PotentialSegregant alleles make up
            de novo dominant variation in the same transcript for the
            minimum number of families. Adds labels to INFO fields of
            VCF records and returns an OrderedDict of 'var_ids' to
            lists of PotentialSegregant objects that appear to
            constitute de novo dominant variation.

            Clears the cache of stored PotentialSegregant alleles.
        '''
        sds = OrderedDict()
        feat_processed = []
        if not self._potential_denovos:
            # if cache is empy, we never encountered the next set of features
            self._potential_denovos = self._last_added
            self._last_added = OrderedDict()
        elif final:
            for feat in self._last_added:
                if feat in self._potential_denovos:
                    self._potential_denovos[feat].update(
                                                        self._last_added[feat])
                else:
                    self._potential_denovos[feat] = self._last_added[feat]
            self._last_added = OrderedDict()
        for feat, pds in self._potential_denovos.items():
            if feat in self._current_features:  # still processing this feature
                continue
            feat_fams = set()
            feat_processed.append(feat)
            for pid, p in pds.items():
                feat_fams.update(p.families)
            if len(feat_fams) >= self.min_families:
                for p in pds.values():
                    samps = (x for x in self.affected
                             if self.ped.fid_from_iid(x) in p.families)
                    if p.alt_id in sds:
                        sds[p.alt_id].add_samples(samps, p.families,
                                                  'samples', [feat], [])
                    else:
                        sv = SegregatingVariant(p, samps, p.families,
                                                'samples', [feat], [],
                                                self.prefix)
                        sds[p.alt_id] = sv
        var_to_segregants = OrderedDict()
        for sv in sds.values():
            sv.annotate_record(self.report_file, self.annot_fields)
            if sv.segregant.var_id in var_to_segregants:
                var_to_segregants[sv.segregant.var_id].append(sv.segregant)
            else:
                var_to_segregants[sv.segregant.var_id] = [sv.segregant]
        #clear the cache of processed features
        for feat in feat_processed:
            del self._potential_denovos[feat]
        return var_to_segregants


class ControlFilter(SampleFilter):
    ''' Filter variants if they are present in a control sample. '''

    def __init__(self, vcf, family_filter, gt_args, n_controls=0):
        '''
            Args:
                vcf:    Input VcfReader object.

                family_filter:
                        FamilyFilter object containing information on
                        which samples are controls in the input VCF.

                gt_args:
                        A dict of arguments to use for filtering
                        genotypes. These should all correspond to
                        arguments to provide to SampleFilter objects.

                n_controls:
                        Minimum number of controls required to carry an
                        ALT allele for it to be filtered. Alleles will
                        only be filtered if carried by this number of
                        controls or more. Default=0.

        '''
        if n_controls and n_controls > len(family_filter.vcf_unaffected):
            n_controls = len(family_filter.vcf_unaffected)
        super().__init__(vcf, controls=family_filter.vcf_unaffected,
                         n_controls=n_controls, confirm_missing=False,
                         **gt_args)


class SegregatingVariant(object):
    '''
        Stores details of alleles that segregate in a manner consistent
        with inheritance pattern.
    '''

    __slots__ = ['recessive', 'samples', 'families', 'model', 'features',
                 'segregant', 'prefix', 'de_novos']

    def __init__(self, segregant, samples, families, model, features,
                 de_novos=(), prefix='VASE_segregant'):
        '''
            Initialize with a PotentialSegregant object, an iterable of
            sample IDs carrying the PotentialSegregant a string
            indicating the model of inheritance (e.g. 'compound_het'),
            the name of the associated features (e.g. transcript IDs),
            prefix for INFO fields and a list of individuals for whom
            the allele appears to have arisen de novo.
        '''
        self.segregant = segregant
        self.samples = list(samples)
        self.families = set(families)
        self.model = [model] * len(self.samples)
        self.features = set(features)
        self.prefix = prefix
        self.de_novos = set(de_novos)

    def __eq__(self, other):
        return self.segregant == other.segregant

    def __hash__(self):
        return hash(self.segregant)

    def add_samples(self, samples, families, model, features, de_novos):
        ''' Add samples with corresponding model of inheritance '''
        self.samples.extend(samples)
        self.families.update(families)
        self.model.extend([model] * (len(self.samples) - len(self.model)))
        self.features.update(features)
        self.de_novos.update(de_novos)

    def annotate_record(self, report_file=None, annot_order=[]):
        ''' Add INFO field annotations for VcfRecords '''
        annots = defaultdict(set)
        for i in range(len(self.model)):
            k = self.prefix
            if self.model[i]:
                k += "_" + self.model[i]
            annots[k].add(self.samples[i])
        for k in annots:
            annots[k] = str.join("|", sorted(annots[k]))
        annots[self.prefix + '_families'] = str.join("|",
                                                     sorted(self.families))
        annots[self.prefix + '_features'] = str.join("|",
                                                     sorted(self.features))
        if self.de_novos:
            annots[self.prefix + '_de_novo'] = str.join("|",
                                                        sorted(self.de_novos))
        converted = self._convert_annotations(annots)
        for k, v in converted.items():
            self.segregant.record.info[k] = v
        if report_file:
            report_file.write(self._annot_to_string(annots, annot_order)
                              + "\n")

    def _annot_to_string(self, annots, annot_order):
        s = ''
        csq_to_join = []
        for k in (x for x in self.segregant.csqs[0] if x != 'Allele'):
            csq_to_join.append(str.join("|", (str(self.segregant.csqs[i][k])
                                              if self.segregant.csqs[i][k]
                                              else '.' for i in range(
                                                   len(self.segregant.csqs)))))
        s = str.join("\t", csq_to_join)
        if annot_order:
            annot_order = [self.prefix + "_" + x for x in annot_order]
            s += "\t" + str.join("\t", (annots[k] if isinstance(annots[k], str)
                                        else '.' for k in annot_order))
        else:
            s += "\t" + str.join("\t", (annots[k] if isinstance(annots[k], str)
                                 else '.' for k in sorted(annots)))
        r = self.segregant.record
        allele = r.alleles[self.segregant.allele]
        s += "\t" + str.join("\t", (str(x) for x in (r.chrom, r.pos, r.id,
                                                     r.ref, r.alt, allele,
                                                     r.qual, r.filter_string)))
        return s

    def _convert_annotations(self, annots):
        ''' Convert to per-allele (Number=A) format for INFO field '''
        converted_annots = dict()
        for k, v in annots.items():
            if k in self.segregant.record.info:
                allele_fields = list(self.segregant.record.info[k])
            else:
                allele_fields = ['.'] * len(self.segregant.record.alts)
            i = self.segregant.allele - 1
            allele_fields[i] = v
            converted_annots[k] = allele_fields
        return converted_annots


class PotentialSegregant(object):
    '''
        Class for storing variant details for records that might make up
        biallelic variants in affected samples.
    '''

    __slots__ = ['allele', 'allele_counts', 'features', 'families', 'alt_id',
                 'var_id', 'record', 'csqs']

    def __init__(self, record, allele, csqs, allele_counts, families):
        self.allele = allele
        self.allele_counts = allele_counts
        self.families = families
        self.var_id = "{}:{}-{}/{}".format(record.chrom, record.pos,
                                           record.ref, record.alt)
        self.alt_id = "{}:{}-{}/{}".format(record.chrom, record.pos,
                                           record.ref, record.alleles[allele])
        self.features = set(x['Feature'] for x in csqs if x['Feature'] != '')
        if not self.features:
            # if is intergenic and there is no Feature ID, use var ID
            # this way we can capture variants at same site if looking for n>1
            # in several families, but won't classify all intergenic variants
            # as the same "Feature"
            self.features.add(self.var_id.replace(',', '_'))
        self.csqs = csqs
        self.record = record

    def __eq__(self, other):
        return self.alt_id == other.alt_id

    def __hash__(self):
        return hash(self.alt_id)
