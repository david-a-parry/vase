import sys
sys.path.insert(0, '')
from vape.parse_vcf.parse_vcf import * 
from vape.filter import * 
from Bio import bgzf


def run(vape_args):
    ''' Run VCF filtering/annotation using args from bin/vape.py '''

    vcf_input = VcfReader(vape_args.input)
    output = get_output(vape_args)
    filters = get_filter_classes(vape_args)
    print_header(vcf_input, output, vape_args)
    var_count = 0
    for record in vcf_input.parser:
        process_record(record, var_count, output, filters )
        var_count += 1
        #sys.stderr.write('\r{} variants processed...\r' .format(var_count))
    sys.stderr.write('\rFinished processing {} variants.\n' .format(var_count))
    finish_up()
    output.close()

def process_record(record, var_count, out, filters=[]):
    remove_alleles = [False] * (len(record.ALLELES) -1)
    keep_alleles = [False] * (len(record.ALLELES) -1)
    # remove_alleles indicates whether allele should be filtered; keep_alleles 
    # indicates whether allele should be kept, overriding any indications in 
    # remove_alleles (e.g. if labelled pathogenic in ClinVar) 
    for f in filters:
        r, k = f.annotate_and_filter_record(record)
        for i in range(len(r)):
            # should only overwrite value of remove_alleles[i] or 
            # keep_alleles[i] with True, not False (e.g. if already set to be 
            # filtered because of a freq in ExAC we shouldn't set to False just 
            # because it is absent from dbSNP)
            if r[i] and not remove_alleles[i]:
                remove_alleles[i] = True
            if k[i] and not keep_alleles[i]:
                keep_alleles[i] = True
    # TODO
    # if all ALTs for a record are set to be filtered and there is no override 
    # in keep_alleles, whole record can be filtered
    #
    # otherwise, mark any filtered alleles and keep record

    # the code below is a placeholder, assumes no inheritance checking etc., 
    # just a simple SNP filter
    for i in range(len(remove_alleles)):
        if not remove_alleles[i] or keep_alleles[i]:
            out.write(str(record) + '\n')
            break

def finish_up():
    pass

def get_filter_classes(vape_args):
    filters = []
    if vape_args.dbsnp:
        kwargs = {"dbsnp_vcf" : vape_args.dbsnp}
        if vape_args.freq is not None:
            kwargs["freq"] = vape_args.freq
        dbsnp_filter = dbSnpFilter(**kwargs)
        filters.append(dbsnp_filter)
    #TODO get gnomAD/ExAC filters
    #TODO get other VCF filters
    return filters

def print_header(vcf, out, vape_args):
    ''' 
        Write a VCF header for output that consists of the input VCF
        header data plus program arguments and any new INFO/FORMAT 
        fields.
    '''

    vape_opts = []
    for k,v in vars(vape_args).items():
        vape_opts.append('--{} {}'.format(k, v))
    prog_header = '##vape.py="' + str.join(" ", vape_opts) + '"'
    new_info = get_vape_info_fields(vape_args)
    out.write(str.join("\n", vcf.header.meta_header + new_info + 
                             [prog_header]) + "\n")
    out.write(str.join("\t", vcf.col_header) + "\n")

def get_vape_info_fields(vape_args):
    return []#TODO
    pass 

def get_output(vape_args):
    ''' 
        Return an output filehandle. If no output specified return 
        sys.stdout, else, if output name ends with .gz or .bgz return a 
        bgzf.BgzfWriter object and otherwise return a standard 
        filehandle.
    '''

    fh = None
    if isinstance(vape_args.output, str):
        if vape_args.output.endswith(('.gz', '.bgz')):
            fh = bgzf.BgzfWriter(vape_args.output)
        else:
            fh = open(vape_args.output, 'w')
    else:
        fh = sys.stdout
    return fh
