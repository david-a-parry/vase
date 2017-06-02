from .parse_vcf.parse_vcf import * 
from .ped_file import *
import logging
from collections import OrderedDict, defaultdict

class FamilyFilter(object):
    ''' 
        Determine whether variants/alleles fit given inheritance 
        patterns for families.
    
    '''

    def __init__(self, ped, vcf, inheritance_pattern=None, gq=0,
                 logging_level=logging.WARNING):
        ''' 
            Initialize with Family object from ped_file.py and a 
            VcfReader object from parse_vcf.py. You may also specify an
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

                inheritance_pattern:
                        Optionally specify the inheritance pattern(s) to 
                        test - either 'dominant' or 'recessive' is 
                        allowed. If not specified, the appropriate
                        patterns will be inferred from family structure.

                gq:     Minimum genotype quality score. Genotype calls 
                        with a GQ lower than this value will be treated
                        as no-calls.

                logging_level: 
                        The level at which logging messages are 
                        displayed. Defaults to logging.WARNING
    
        '''
        self.logger = logging.getLogger(__name__)
        if not self.logger.hasHandlers():
            self.logger.setLevel(logging_level)
            formatter = logging.Formatter(
                        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            ch = logging.StreamHandler()
            ch.setLevel(self.logger.level)
            ch.setFormatter(formatter)
            self.logger.addHandler(ch)
        self.affected = list(ped.get_affected())
        self.unaffected = list(ped.get_unaffected())
        self.ped = ped
        self.vcf = vcf
        if not self.affected:
            raise Exception("No affected individuals found in PED file '{}'"
                            .format(ped.filename))
        self.vcf_affected = []
        self.vcf_unaffected = []
        self.vcf_samples = []
        for aff in self.affected:
            if aff in self.vcf.header.samples:
                self.vcf_affected.append(aff)
        if not self.vcf_affected:
            raise Exception("No affected individuals in PED file '{}'"
                            .format(ped.filename) + " found in VCF '{}'")
        for un in self.unaffected:
            if un in self.vcf.header.samples:
                self.vcf_unaffected.append(un)
        self.vcf_samples = self.vcf_affected + self.vcf_unaffected
        if inheritance_pattern is None:
            self.inheritance_patterns = dict()
            self._infer_inheritance()
        else:
            #TODO - assign inheritance for each family 
            pass

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
            dominant = False
            denovo = False
            recessive = False
            self.logger.info("Assessing inheritance pattern of family {}"
                              .format(fid))
            f_aff =list(fam.get_affected())
            if not f_aff:
                continue
            for iid in f_aff:
                self.logger.info("Checking affected individual {}" .format(iid))
                n_affected += 1
                indv = fam.individuals[iid]
                if not indv.parents:
                    self.logger.info("No parents for affected individual {}"
                                      .format(iid))
                    continue
                no_parents = False
                obligate_carriers = set()
                for par in indv.parents:
                    #is parent affected
                    parent = fam.individuals[par]
                    par_to_child = False
                    gpar_to_child = False
                    if parent.is_affected():
                        self.logger.info("Apparent vertical transmission from " + 
                                          "{} -> {}" .format(par, iid))
                        par_to_child = True
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
            if not dominant:
                recessive = True
            if recessive and n_affected == 1 and not no_parents:
                denovo = True
            elif recessive and n_affected > 1:
                # we can entertain apparent de novos due to somatic mosaicism
                # if all affecteds share a parent
                pars = fam.individuals[f_aff[0]].parents
                shared_pars = None
                if len(pars) == 2:
                    shared_pars = set(pars)
                    for i in range(1, len(f_aff)):
                        ipars = self.ped.individuals[f_aff[i]].parents
                        if ipars == None: 
                            break
                        shared_pars = shared_pars.intersection(ipars)
                        if not shared_pars:
                            break
                if shared_pars:
                    denovo = True
            #TODO - Check who is present in VCF and arrange filtering appropriately
            
            self.inheritance_patterns[fid] = []
            if recessive:
                self.logger.info("Family '{}' " .format(fid) + "can be " + 
                                 "analysed  under a recessive model")
                self.inheritance_patterns[fid].append('recessive')
            if denovo:
                dmodel = "de novo"
                if n_affected > 1:
                    dmodel += " (with germline mosaicism)"
                self.logger.info("Family '{}' " .format(fid) + "can be " + 
                                 "analysed  under a {} model" .format(dmodel))
                self.inheritance_patterns[fid].append('de_novo')
            if dominant:
                self.logger.info("Family '{}' " .format(fid) + "can be " + 
                                 "analysed  under a dominant model")
                self.inheritance_patterns[fid].append('dominant')
            


class RecessiveFilter(object):
    ''' 
        This class assumes that each family has a shared biallelic 
        genetic cause of disease. It will not cope with phenocopies, 
        pseudodominance or other more complicated inheritance patterns.
    '''
    def __init__(self, vcf, family_filter, gq=0, strict=False, 
                 exclude_denovo=False):
        '''
            Arguments:
                    vcf:    VcfReader object from parse_vcf.py
                    
                    family_filter: 
                            FamilyFilter object

                    gq:     Minimum genotype quality score. Genotype 
                            calls with a GQ lower than this value will
                            be treated as no-calls.

                    strict: If True, for any affected sample with 
                            parents, require confirmation of parental 
                            genotypes. If either parent genotype is a 
                            no-call for a record, then the record will 
                            be ignored. Default = False.

                    exclude_denovo:
                            If True, where there is data available from
                            both parents for an affected individual 
                            ignore apparent de novo occuring alleles.
                            
        '''

        self.family_filter = family_filter
        self.ped = family_filter.ped
        self.families = list(x for x in family_filter.inheritance_patterns if
                          'recessive' in family_filter.inheritance_patterns[x])
        self.samples = family_filter.vcf_samples
        self.unaffected = family_filter.vcf_unaffected
        self.affected = list(x for x in family_filter.vcf_affected if 
                             self.ped.individuals[x].fid in self.families)
        self.gq = gq
        self._fam_to_aff = dict()
        for fid in self.families:
            self._fam_to_aff[fid] = set(x for x in 
                                        self.ped.families[fid].get_affected()
                                        if x in self.affected)
        self._prev_coordinate = (None, None)  # to ensure records are processed 
        self._processed_contigs = set()       # in coordinate order
        self.strict = strict
        self.exclude_denovo = exclude_denovo
        self._potential_recessives = dict()
        self._last_added = dict()
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
                record: VcfRecord from parse_vcf.py

                
                        ignore_alleles:
                        List of booleans indicating for each ALT in 
                        order whether it should be ignored in relation
                        to possible recessive variation (e.g. if MAF is
                        too high, no likely pathogenic consequence 
                        etc.). This will normally have been generated 
                        by VapeRunner via VcfFilter and/or VepFilter
                        classes.

                ignore_csq:
                        List of booleans indicating for each CSQ in 
                        order whether it should be ignored in relation
                        to possible recessive variation. This should 
                        normally have been generated by a corresponding 
                        VepFilter object.

        '''
        stored = False
        if self._prev_coordinate[0] != record.CHROM:
            if record.CHROM in self._processed_contigs:
                raise Exception("Input must be sorted by chromosome and " + 
                                "position for recessive filtering. " + 
                                "Contig '{}' " .format(record.CHROM) + 
                                "encountered before and after contig " + 
                                "'{}'." .format(self._prev_coordinate[0]))
            if self._prev_coordinate[0] is not None:
                self._processed_contigs.add(self._prev_coordinate[0])
        elif record.POS < self._prev_coordinate[1]:
            raise Exception("Input must be sorted by chromosome and position" + 
                            "for recessive filtering. Encountered position " + 
                            "{}:{} after {}:{}" .format(record.CHROM, 
                             record.POS, self._prev_coordinate[0], 
                             self._prev_coordinate[1]))
        self._prev_coordinate = (record.CHROM, record.POS)
        gts = record.parsed_gts(fields=['GT', 'GQ'], samples=self.samples)
        skip_fam = set()
        for i in range(len(record.ALLELES) -1):
            if ignore_alleles and ignore_alleles[i]:
                continue
            alt = i + 1
            skip_allele = False
            fams_with_allele = []
            for un in self.unaffected:
                if gts['GT'][un] == (alt, alt) and gts['GQ'][un] >= self.gq:
                    #hom in a control - skip allele
                    skip_allele = True
                    break
            if skip_allele:
                continue
            for fid in self.families:
                if fid in skip_fam:
                    continue
                have_allele = set() #affecteds carrying this allele
                for aff in self._fam_to_aff[fid]:
                    #check all affecteds carry this allele
                    if alt in gts['GT'][aff] and gts['GQ'][aff] >= self.gq:
                        have_allele.add(aff)
                    else: 
                        break
                if have_allele == self._fam_to_aff[fid]: 
                    #all affecteds in family carry allele
                    fams_with_allele.append(fid)
            if fams_with_allele:
                #store record and consequences 
                try:
                    csqs = []
                    for j in range(len(record.CSQ)):
                        if ignore_csq and ignore_csq[j]:
                            continue
                        if record.CSQ[j]['alt_index'] == alt:
                            #store record and csq details
                            csqs.append(record.CSQ[j])
                    if csqs:
                        stored = True
                        self._last_added = OrderedDict()
                        alt_counts = self._get_allele_counts(alt, gts)
                        pr = PotentialRecessive(record=record, allele=alt,
                                                csqs=csqs, 
                                                allele_counts=alt_counts, 
                                                families=fams_with_allele)
                        for feat in pr.features:
                            self._last_added[feat] = OrderedDict(
                                                             [(pr.alt_id, pr)])
                            if feat in self._potential_recessives:
                                self._potential_recessives[feat][pr.alt_id] = pr
                            else:
                                self._potential_recessives[feat] = OrderedDict(
                                    [(pr.alt_id, pr)])
                            
                        
                except HeaderError:
                    raise Exception("Could not identify CSQ or ANN fields " +
                                    "in VCF header. Please ensure your input" +
                                    " is annotated with Ensembl's VEP to " +
                                    "perform recessive filtering")
        return stored

    def process_potential_recessives(self):
        ''' 
            Check whether stored PotentialRecessive alleles make up 
            biallelic variation in the same transcript for affected 
            individuals/families. Adds labels to INFO fields of VCF 
            records and returns a list of 'var_ids' of variants that 
            constitute biallelic variation.

            Clears the cache of stored PotentialRecessive alleles.
        '''
        segregating = dict() #keys are alt_ids, values are SegregatingBiallelic 
        for feat, prs in self._potential_recessives.items():
            un_hets = defaultdict(list)  #store het alleles carried by each unaffected
            aff_hets = defaultdict(list) #store het alleles carried by each affected
            biallelics = defaultdict(list)  #store biallelic combinations for affecteds
            for pid,p in prs.items():
                for un in self.unaffected:
                    if p.allele_counts[un] == 1:  #already checked for homs when adding
                        #store allele carried in this unaffected 
                        un_hets[un].append(pid)
                for aff in (x for x in self.affected 
                            if self.ped.fid_from_iid(x) in p.families):
                    if p.allele_counts[aff] == 1:
                        aff_hets[aff].append(pid)
                    elif p.allele_counts[aff] == 2:
                        biallelics[aff].append(tuple([pid]))
            incompatibles = [] #create a list of sets of incompatible hets
            for hets in un_hets.values():
                if len(hets):
                    incompatibles.append(set(hets))
            for aff,hets in aff_hets.items():
                for i in range(len(hets)):
                    for j in range(i+1, len(hets)):
                        incomp = False
                        for iset in incompatibles:
                            if iset.issuperset([hets[i], hets[j]]):
                                incomp = True
                                break
                        if not incomp:
                            #TODO - check phase groups!
                            biallelics[aff].append(tuple([hets[i], hets[j]]))
            if not biallelics:
                continue
            #see if all affecteds in the same family share the same biallelics
            for fid,affs in self._fam_to_aff.items():
                b_affs = set(x for x in affs if x in biallelics)
                if len(b_affs) == 0 or b_affs != affs:
                    continue
                affs = list(affs)
                for i in range(len(affs)):
                    for bi in biallelics[affs[i]]:
                        absent_in_aff = False
                        for j in range(i+1, len(affs)):
                            if bi not in biallelics[affs[j]]:
                                absent_in_aff = True
                                break
                        if not absent_in_aff:
                            segs,de_novo = self._check_parents(feat, bi, affs)
                            if not segs:
                                continue
                            if len(bi) == 1:
                                model = 'homozygous'
                            else:
                                model = 'compound_het'
                            for bi_pr in (prs[x] for x in bi):
                                if bi_pr in segregating:
                                    segregating[bi_pr].add_samples(affs, model, 
                                                   feat, de_novo[bi_pr.alt_id])
                                else:
                                    segregating[bi_pr] = SegregatingBiallelic(
                                                      bi_pr, affs, model, feat,
                                                      de_novo[bi_pr.alt_id])
        b_var_ids = set()
        for sb in segregating.values():
            sb.annotate_records()
            b_var_ids.add(sb.recessive.var_id)
        #clear the cache except for the last entry which will be a new gene
        self._potential_recessives = self._last_added
        self._last_added = dict()
        return b_var_ids

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
        if len(counts) == 1:#homozygous
            counts.append(counts[0])
        for samp in samples:
            parents = self.ped.individuals[samp].parents
            par = list(x for x in parents if x in self.samples)
            if len(par) == 0:
                continue
            if self.strict:
                for p in par:
                    if None in (counts[i][p] for i in range(len(counts))):
                        #require both parental genotypes if self.strict
                        return (False, dns)
            if len(par) == 2: #can check for de novos
                for i in range(len(counts)):
                    if counts[i][par[0]] == 0 and counts[i][par[1]] == 0:
                        #apparent de novo
                        self.family_filter.logger.debug("Apparent de novo " + 
                            "allele {} for sample {} (parents = {} + {}) "
                            .format(alleles[-i], samp, par[0], par[1]) + 
                            "for recessive combination {}/{}"
                            .format(alleles[0], alleles[-1]))
                        dns[alleles[-i]].append(samp)
                        if self.exclude_denovo:
                            return (False, dns)
            elif len(par) == 1:
                # if only one parent and both alleles are absent it is more
                # likely that the two alleles are in cis from other parent
                if counts[0][par[0]] == 0 and counts[1][par[0]] == 0:
                    return(False, dns)
            #NOTE: we could do a check here to ensure that any non-affected
            #      parent does not carry both alleles, but this *SHOULD* have
            #      already been done earlier in process_potential_recessives
            #      function for ALL unaffecteds anyway
        return (True, dns)
 
    def _get_allele_counts(self, allele, gts):
        a_counts = dict()
        for samp in gts['GT']:
            gq = gts['GQ'][samp] 
            if gq is not None and gq >= self.gq:
                a_counts[samp] = gts['GT'][samp].count(allele)
            else:
                a_counts[samp] = None
        return a_counts

        
class SegregatingBiallelic(object):
    '''
        Stores details of alleles that segregate in a manner consistent
        with recessive inheritance.
    '''

    __slots__ = ['recessive', 'samples', 'model', 'features', 'de_novos']

    def __init__(self, recessive, samples, model, feature, de_novos):
        ''' 
            Initialize with a PotentialRecessive object, an iterable of
            sample IDs carrying the PotentialRecessives a string 
            indicating the model of inheritance (e.g. 'compound_het'),
            the name of the associated feature (e.g. a transcript 
            ID) and a list of individuals for whom the allele appears 
            to have arisen de novo.
        '''
        self.recessive = recessive
        self.samples = list(samples)
        self.model = [model] * len(self.samples)
        self.features = set([feature])
        self.de_novos = set(de_novos)

    def __eq__(self, other):
        return self.recessive == other.recessive

    def __hash__(self):
        return hash(self.recessive)

    def add_samples(self, samples, model, feature, de_novos):
        ''' Add samples with corresponding model of inheritance '''
        self.samples.extend(samples)
        self.model.extend([model] * (len(self.samples) - len(self.model)))
        self.features.add(feature)
        self.de_novos.update(de_novos)

    def annotate_records(self):
        ''' Add VapeBiallelic INFO field annotations for VcfRecords '''
        annots = defaultdict(set)
        for i in range(len(self.model)):
            k = 'VAPE_biallelic_' + self.model[i]
            annots[k].add(self.samples[i])
        for k in annots:
            annots[k] = str.join(",", sorted(annots[k]))
        annots['VAPE_biallelic_features'] = str.join(",", 
                                                     sorted(self.features))
        if self.de_novos:
            annots['VAPE_biallelic_de_novo'] = str.join(",", 
                                                        sorted(self.de_novos))
        self.recessive.record.add_info_fields(annots)


class PotentialRecessive(object):
    ''' 
        Class for storing variant details for records that might make up
        biallelic variants in affected samples.
    '''

    __slots__ = ['allele', 'allele_counts', 'features', 'families', 'alt_id',
                 'var_id', 'record']

    def __init__(self, record, allele, csqs, allele_counts, families):
        self.allele = allele
        self.allele_counts = allele_counts
        self.families = families
        self.features = set(x['Feature'] for x in csqs)
        self.record = record
        self.var_id = "{}:{}-{}/{}".format(record.CHROM, record.POS, 
                                           record.REF, record.ALT)
        self.alt_id = "{}:{}-{}/{}".format(record.CHROM, record.POS, 
                                           record.REF, record.ALLELES[allele])
    def __eq__(self, other):
        return self.alt_id == other.alt_id

    def __hash__(self):
        return hash(self.alt_id)

class DominantFilter(object):
    pass


class DeNovoFilter(object):
    pass


