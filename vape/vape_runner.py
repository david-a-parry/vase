import sys
import re
import logging
from .parse_vcf.parse_vcf import * 
from .dbsnp_filter import * 
from .gnomad_filter import * 
from .vep_filter import * 
from .sample_filter import *
from .ped_file import *

class VapeRunner(object):

    def __init__(self, args):

        self.args = args
        self._set_logger()
        self.input = VcfReader(self.args.input)
        self.prev_annots, self.info_prefixes = self._get_prev_annotations()
        self.new_annots = dict()
        self.out = self.get_output()
        self.vcf_filters = self.get_vcf_filter_classes()
        self.ped = None
        if args.ped:
            self.ped = PedFile(args.ped)
        self.csq_filter = None
        if args.csq is not None:
            self.csq_filter = VepFilter(args.csq, args.canonical,
                                        args.biotypes, 
                                        args.missense_filters,
                                        args.filter_unpredicted,
                                        args.keep_if_any_damaging)
        self.sample_filter = None
        if args.cases or args.controls:
            self.sample_filter = SampleFilter(self.input, args.cases, 
                                              args.controls, args.n_cases,
                                              args.n_controls, args.gq)
        self.de_novo_filters = None
        if args.de_novo:
            self._get_de_novo_filters()
            
            

    def run(self):
        ''' Run VCF filtering/annotation using args from bin/vape.py '''
        self.logger.info('Starting variant processing')
        self.print_header()
        var_count = 0
        for record in self.input.parser:
            self.process_record(record)
            var_count += 1
            if not self.args.quiet:
                sys.stderr.write('\r{} variants processed...\r' 
                                 .format(var_count))
        self.logger.info('Finished processing {} variants.' 
                             .format(var_count))
        self.finish_up()
        self.out.close()

    def process_record(self, record):
        if self.filter_global(record):
            return
        filter_alleles = self.filter_alleles_external(record)
        if sum(filter_alleles) == len(filter_alleles): 
            #all alleles should be filtered
            return
        if self.sample_filter:
            for i in range(1, len(record.ALLELES)):
                r = self.sample_filter.filter(record, i)
                if r:
                    filter_alleles[i-1] = True
                if sum(filter_alleles) == len(filter_alleles): 
                    #all alleles should be filtered
                    return
        filter_alleles = self.inhertance_filter(record, filter_alleles)
        if sum(filter_alleles) == len(filter_alleles): 
            return
        self.out.write(str(record) + '\n')
    
    def inhertance_filter(self, record, ignore_alleles):
        remove_alleles = [False] * (len(record.ALLELES) - 1)
        keep_alleles = [False] * (len(record.ALLELES) - 1)
        if self.de_novo_filters:
            found_de_novo = False
            de_novos = []
            for i in range(len(record.ALLELES) - 1):
                if ignore_alleles[i]:
                    de_novos.append('.')
                    continue
                for f in self.de_novo_filters:
                    r = f.filter(record, i + 1)
                    if r:
                        remove_alleles[i] = True
                    else:
                        found_de_novo = True
                        keep_alleles[i] = True
                        if len(de_novos) > i:
                            de_novos[i] += '|' + f.child
                        else:
                            de_novos.append(f.child)
                
            if found_de_novo:
                inf = {'VAPE_de_novo' : str.join(",", de_novos)}
                record.add_info_fields(info=inf, append_existing=True)
        for i in range(len(remove_alleles)):
            if remove_alleles[i] and not keep_alleles[i]:
                ignore_alleles[i] = True
        return ignore_alleles

    def finish_up(self):
        pass

    def filter_alleles_external(self, record):
        ''' 
            Return True or False for each allele indicating whether an 
            allele should be filtered based on information from VEP, 
            dbSNP, ClinVar or gnomAD.
        '''

        remove_alleles = [False] * (len(record.ALLELES) -1)
        keep_alleles = [False] * (len(record.ALLELES) -1)
        matched_alleles = [False] * (len(record.ALLELES) -1)
        # remove_alleles indicates whether allele should be filtered; 
        # keep_alleles indicates whether allele should be kept, overriding any 
        # indications in remove_alleles (e.g. if labelled pathogenic in 
        # ClinVar)
        # check functional consequences
        if self.csq_filter:
            r_alts, r_csq = self.csq_filter.filter(record)
            for i in range(len(r_alts)):
                if r_alts[i]:
                    remove_alleles[i] = True
            if (not self.args.clinvar_path and 
                sum(remove_alleles) == len(remove_alleles)):
                # bail out now if no valid consequence and not keeping clinvar
                # path variants - if using clinvar path we have to ensure we 
                # haven't got a path variant with a non-qualifying consequence
                return remove_alleles
        for f in self.vcf_filters:
            r, k, m = f.annotate_and_filter_record(record)
            for i in range(len(r)):
                # should only overwrite value of remove_alleles[i] or 
                # keep_alleles[i] with True, not False (e.g. if already set to 
                # be filtered because of a freq in ExAC we shouldn't set to 
                # False just because it is absent from dbSNP)
                if r[i]:
                    remove_alleles[i] = True
                if k[i]:
                    keep_alleles[i] = True
                if m[i]:
                    matched_alleles[i] = True
        
                verdict = []
        for i in range(len(remove_alleles)):
            if remove_alleles[i] and not keep_alleles[i]:
                verdict.append(True)
            elif self.args.filter_novel and not matched_alleles[i]:
                verdict.append(True)
            else:
                verdict.append(False)
        return verdict
 
    def filter_global(self, record):
        ''' Return True if record fails any global variant filters.'''
        if self.args.pass_filters:
            if record.FILTER != 'PASS':
                return True
        if self.args.variant_quality is not None:
            if record.QUAL < self.args.variant_quality:
                return True
        return False

    def get_vcf_filter_classes(self):
        filters = []
        uni_args = {}
        if self.args.freq is not None:
            uni_args["freq"] = self.args.freq
        if self.args.min_freq is not None:
            uni_args["min_freq"] = self.args.min_freq
        # get dbSNP filter
        for dbsnp in self.args.dbsnp:
            prefix = self.check_info_prefix('VAPE_dbSNP')
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
                self.input.header.addHeaderField(name=f, dictionary=d, 
                                          field_type='INFO')
        # get gnomAD/ExAC filters
        for gnomad in self.args.gnomad:
            prefix = self.check_info_prefix('VAPE_gnomAD')
            kwargs = {"vcf" : gnomad, "prefix" : prefix}
            kwargs.update(uni_args)
            gnomad_filter = GnomadFilter(**kwargs)
            filters.append(gnomad_filter)
            for f,d in gnomad_filter.added_info.items():
                self.logger.debug("Adding gnomAD/ExAC annotation {}" 
                                  .format(f))
                self.input.header.addHeaderField(name=f, dictionary=d, 
                                          field_type='INFO')
        #TODO get other VCF filters
        return filters

    def _get_prev_annotations(self):
        annots = set()
        prefixes = set()
        for info in self.input.metadata['INFO']:
            match = re.search('^(VAPE_\w+)_\w+(_\d+)?', info)
            if match:
                annots.add(info)
                prefixes.add(match.group(1))
                self.logger.debug("Identified previously annotated VAPE INFO" +
                                  " field '{}'" .format(info))
        return annots, prefixes

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

        vape_opts = []
        for k,v in vars(self.args).items():
            vape_opts.append('--{} {}'.format(k, v))
        self.input.header.addHeaderField(name="vape.py", 
                                   string='"' + str.join(" ", vape_opts) + '"')
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
                    raise Exception("Can not import bgzf via " + 
                                    "biopython. Please install biopython in" +
                                    " order to write bgzip compressed " + 
                                    "(.gz/.bgz) output.")
                fh = bgzf.BgzfWriter(self.args.output)
            else:
                fh = open(self.args.output, 'w')
        else:
            fh = sys.stdout
        return fh

    def _get_de_novo_filters(self):
        if not self.ped:
            raise Exception("--de_novo argument requires a PED file but no " +
                            "PED was specified. Please specify a PED file " + 
                            "with at least one parent-child trio using the " + 
                            "--ped argument.")

        fams = set()
        self.de_novo_filters = []
        num_trios = 0
        for iid in self.ped.get_affected():
            if iid in self.input.header.samples:
                self.logger.debug("Ped sample {} is in VCF" .format(iid))
                if (self.ped.individuals[iid].mother and 
                    self.ped.individuals[iid].father):
                    #must have both parents specified present in PED 
                    if (self.ped.individuals[iid].mother in 
                        self.input.header.samples and 
                        self.ped.individuals[iid].father in 
                        self.input.header.samples):
                        #and must have both parents present in VCF
                        df = DeNovoFilter(vcf=self.input, child=iid,
                                       mother=self.ped.individuals[iid].mother,
                                       father=self.ped.individuals[iid].father,
                                       gq=self.args.gq)
                        self.de_novo_filters.append(df)
                        fams.add(self.ped.individuals[iid].fid)
                        num_trios += 1
                        self.logger.debug("Ped sample {} has " .format(iid) +
                                          "both  parents in VCF and can be "+
                                          "de novo filtered")
                else:
                    self.logger.debug("Ped sample {} does not " .format(iid) + 
                                      "have both parents in VCF")
            else:
                self.logger.debug("Ped sample {} not in VCF" .format(iid))
        self.logger.info("Identified {} parent-child trios " .format(num_trios) 
                  + "in {} families present in PED and VCF " .format(len(fams))
                  + "from {} total families in PED." 
                    .format(len(self.ped.families)))
        if not num_trios:
            e = Exception("No parent-child trios identified from PED and VCF" +
                          " - please check your PED and VCF sample IDs")
            self.logger.exception(e)
            raise e
        
        self.input.header.addHeaderField(
            field_type='INFO',
            name='VAPE_de_novo', 
            dictionary={
                         'Number' : 'A',
                         'Type' : 'String',
                         'Description' : '"Sample IDs in which an ALT allele' +
                                         ' appears to have arisen de novo ' +  
                                         'based on the absence of an allele ' +
                                         'call in either parent"',
            })
                
    def _set_logger(self):
        self.logger = logging.getLogger("VAPE")
        if self.args.debug:
            self.logger.setLevel(logging.DEBUG)
        elif self.args.quiet:
            self.logger.setLevel(logging.WARNING)
        else:
            self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter(
                        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setLevel(self.logger.level)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

