import sys
import re
from .parse_vcf.parse_vcf import * 
from .dbsnp_filter import * 
from .gnomad_filter import * 
from Bio import bgzf

class VapeRunner(object):

    def __init__(self, args):

        self.args = args
        self.info_prefixes = set()
        self.input = VcfReader(self.args.input)
        self.out = self.get_output()
        self.filters = self.get_filter_classes()

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
        # remove_alleles indicates whether allele should be filtered; 
        # keep_alleles indicates whether allele should be kept, overriding any 
        # indications in remove_alleles (e.g. if labelled pathogenic in 
        # ClinVar)
        for f in self.filters:
            r, k = f.annotate_and_filter_record(record)
            for i in range(len(r)):
                # should only overwrite value of remove_alleles[i] or 
                # keep_alleles[i] with True, not False (e.g. if already set to 
                # be filtered because of a freq in ExAC we shouldn't set to 
                # False just because it is absent from dbSNP)
                if r[i] and not remove_alleles[i]:
                    remove_alleles[i] = True
                if k[i] and not keep_alleles[i]:
                    keep_alleles[i] = True
        # TODO
        # if all ALTs for a record are set to be filtered and there is no 
        # override in keep_alleles, whole record can be filtered
        #
        # otherwise, mark any filtered alleles and keep record

        # the code below is a placeholder, assumes no inheritance checking etc. 
        # just a simple SNP filter
        for i in range(len(remove_alleles)):
            if not remove_alleles[i] or keep_alleles[i]:
                self.out.write(str(record) + '\n')
                break

    def finish_up(self):
        pass

    def get_filter_classes(self):
        filters = []
        # get dbSNP filter
        for dbsnp in self.args.dbsnp:
            prefix = self.check_info_prefix('VAPE_dbSNP')
            kwargs = {"vcf" : dbsnp, "prefix" : prefix}
            if self.args.freq is not None:
                kwargs["freq"] = self.args.freq
            dbsnp_filter = dbSnpFilter(**kwargs)
            filters.append(dbsnp_filter)
        # get gnomAD/ExAC filters
        for gnomad in self.args.gnomad:
            prefix = self.check_info_prefix('VAPE_gnomAD')
            kwargs = {"vcf" : gnomad, "prefix" : prefix}
            gnomad_filter = GnomadFilter(**kwargs)
            filters.append(gnomad_filter)
        #TODO get other VCF filters
        return filters

    def check_info_prefix(self, name):
        if name in self.info_prefixes:
            match = re.search('_(\d+)$', name) 
                #already has an appeneded '_#' - increment and try again
            if match:
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
        prog_header = '##vape.py="' + str.join(" ", vape_opts) + '"'
        new_info = self.get_vape_info_fields()
        self.out.write(str.join("\n", self.input.header.meta_header + new_info  
                                + [prog_header]) + "\n")
        self.out.write(str.join("\t", self.input.col_header) + "\n")

    def get_vape_info_fields(self):
        return []#TODO
        pass 

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
