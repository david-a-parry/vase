import sys
import re
from .parse_vcf.parse_vcf import * 
from .dbsnp_filter import * 
from .gnomad_filter import * 
from .vep_filter import * 
from Bio import bgzf

class VapeRunner(object):

    def __init__(self, args):

        self.args = args
        self.input = VcfReader(self.args.input)
        self.prev_annots, self.info_prefixes = self._get_prev_annotations()
        self.new_annots = dict()
        self.out = self.get_output()
        self.vcf_filters = self.get_vcf_filter_classes()
        self.csq_filter = None
        if args.csq is not None:
            self.csq_filter = VepFilter(args.csq, args.canonical,
                                        args.biotypes, 
                                        args.missense_filters,
                                        args.filter_unpredicted,
                                        args.keep_if_any_damaging)

    def run(self):
        ''' Run VCF filtering/annotation using args from bin/vape.py '''

        self.print_header()
        var_count = 0
        for record in self.input.parser:
            self.process_record(record, var_count)
            var_count += 1
            sys.stderr.write('\r{} variants processed...\r' .format(var_count))
        sys.stderr.write('\rFinished processing {} variants.\n' 
                         .format(var_count))
        self.finish_up()
        self.out.close()

    def process_record(self, record, var_count):
        remove_alleles = [False] * (len(record.ALLELES) -1)
        keep_alleles = [False] * (len(record.ALLELES) -1)
        matched_alleles = [False] * (len(record.ALLELES) -1)
        # remove_alleles indicates whether allele should be filtered; 
        # keep_alleles indicates whether allele should be kept, overriding any 
        # indications in remove_alleles (e.g. if labelled pathogenic in 
        # ClinVar)
        if self.args.pass_filters:
            if record.FILTER != 'PASS':
                return
        if self.args.variant_quality is not None:
            if record.QUAL < self.args.variant_quality:
                return
        # check functional consequences
        if self.csq_filter:
            r_alts, r_csq = self.csq_filter.filter(record)
            for i in range(len(r_alts)):
                if r_alts[i]:
                    remove_alleles[i] = True
            #bail out now if no valid consequence
            if sum(remove_alleles) == len(remove_alleles):
                return
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
        # TODO
        # if all ALTs for a record are set to be filtered and there is no 
        # override in keep_alleles, whole record can be filtered
        #
        # otherwise, mark any filtered alleles and keep record

        # the code below is a placeholder, assumes no inheritance checking etc. 
        # just a simple SNP filter
        if (sum(remove_alleles) == len(remove_alleles) and 
            sum(keep_alleles) == 0):
            #all alleles should be filtered and no overide in keep_alleles
            return
        if self.args.filter_novel and sum(matched_alleles[i]) == 0:
            return
        self.out.write(str(record) + '\n')

    def finish_up(self):
        pass

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
        return annots, prefixes

    def check_info_prefix(self, name):
        if name in self.info_prefixes:
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

        fh = None
        if isinstance(self.args.output, str):
            if self.args.output.endswith(('.gz', '.bgz')):
                fh = bgzf.BgzfWriter(self.args.output)
            else:
                fh = open(self.args.output, 'w')
        else:
            fh = sys.stdout
        return fh
