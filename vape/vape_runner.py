import sys
sys.path.insert(0, '')
from vape.parse_vcf.parse_vcf import * 
from Bio import bgzf

def run(vape_args):
    ''' Run VCF filtering/annotation using args from bin/vape.py '''

    vcf_input = VcfReader(vape_args.input)
    output = get_output(vape_args)
    filters = get_filter_methods(vape_args)
    print_header(vcf_input, output, vape_args)
    for record in vcf_input.parser:
        process_record(record)
    finish_up()
    output.close()

def process_record(record):
    pass

def finish_up():
    pass

def get_filter_methods(vape_args):
    pass

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
