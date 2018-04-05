import sys
import re
import logging
import xlsxwriter
from collections import namedtuple
from .ped_file import PedFile, Family, Individual, PedError
from parse_vcf import VcfReader, VcfHeader, VcfRecord 
from .ensembl_rest_queries import EnsemblRestQueries

vcf_output_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',] 
feat_annots = { 'VASE_biallelic_families': 'VASE_biallelic_features',
                'VASE_dominant_families': 'VASE_dominant_features',
                'VASE_de_novo_families': 'VASE_de_novo_features',}

class VaseReporter(object):
    ''' Read a VASE annotated VCF and output XLSX format summary of 
        segregating variants. '''

    def __init__(self, vcf, ped, out, families=[], all_features=False,
                 rest_lookups=False, grch37=False, prog_interval=None,
                 timeout=2.0, quiet=False, debug=False):
        self._set_logger(quiet, debug)
        self.vcf = VcfReader(vcf)
        self.ped = PedFile(ped)
        self.families = families
        self.all_features = all_features
        self.seg_fields = self._check_header()
        if not out.endswith(".xlsx"):
            out = out + ".xlsx"
        self.out = out
        self.sample_orders = dict() #key is fam id, value is list of samples
        self.workbook = xlsxwriter.Workbook(out)
        self.bold = self.workbook.add_format({'bold': True})
        self.worksheets = dict()
        self.rows = dict()
        self.rest_lookups = rest_lookups
        if not families:
            families = sorted(self.ped.families.keys())
        for f in families:#respect the order of families provided before converting to set
            if f not in self.ped.families:
                raise RuntimeError("Familiy '{}' not found in ped ".format(f) +
                                   "'{}'.".format(ped))
            if f in self.worksheets:
                logger.warn("Duplicate family '{}' specified".format(f))
            else:
                self.worksheets[f] = self._initialize_worksheet(f)
                self.rows[f] = 1
        self.families = set(families)
        if self.rest_lookups:
            self.ensembl_rest = EnsemblRestQueries(use_grch37_server=grch37,
                                                   timeout=timeout)
        self.rest_cache = dict()
        if prog_interval is not None:
            self.prog_interval = prog_interval
        elif self.rest_lookups:
            self.prog_interval = 100
        else:
            self.prog_interval = 1000

    def _set_logger(self, quiet=False, debug=False):
        self.logger = logging.getLogger("VASE Reporter")
        if debug:
            self.logger.setLevel(logging.DEBUG)
        elif quiet:
            self.logger.setLevel(logging.WARNING)
        else:
            self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter(
                        '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setLevel(self.logger.level)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)


    def _check_header(self):
        inheritance_fields = dict()
        if 'VASE_biallelic_families' in self.vcf.metadata['INFO']:
            inheritance_fields['VASE_biallelic_families'] = 'recessive'
            self.logger.info("Found VASE biallelic annotations.")
            seg = True
        if 'VASE_dominant_families' in self.vcf.metadata['INFO']:
            inheritance_fields['VASE_dominant_families'] = 'dominant'
            self.logger.info("Found VASE dominant annotations.")
            seg = True
        if 'VASE_de_novo_families' in self.vcf.metadata['INFO']:
            inheritance_fields['VASE_de_novo_families'] = 'de novo'
            self.logger.info("Found VASE de_novo annotations.")
            seg = True
        if not seg:
            raise RuntimeError("no vase recessive/dominant/de novo " + 
                               "annotations found in vcf - please run vase " +
                               "with --biallelic, --dominant or --de_novo " + 
                               " options first.")
        return inheritance_fields

    def _initialize_worksheet(self, family):
        worksheet = self.workbook.add_worksheet(family)
        header = ['INHERITANCE'] + vcf_output_columns + ['ALLELE', 'AC', 'AN']
        header.extend(self._get_sample_order(family))
        header.extend(x for x in self.vcf.header.csq_fields if x != 'Allele')
        if self.rest_lookups:
            header.extend(["ENTREZ", "GO", "REACTOME", "MOUSE_TRAITS", 
                           "MIM_MORBID"])
        for i in range(len(header)):
            worksheet.write(0, i, header[i], self.bold)
        return worksheet

    def _get_sample_order(self, family):
        if family in self.sample_orders:
            return self.sample_orders[family]
        samples = []
        children = []
        for par in sorted(self.ped.families[family].parents):
            if par in self.vcf.header.samples:
                samples.append(par)
            for child in self.ped.families[family].parents[par]:
                if child not in children:
                    children.append(child)
        for child in children:
            if child not in samples and child in self.vcf.header.samples:
                samples.append(child)
        self.sample_orders[family] = samples
        return samples

    def get_ensembl_rest_data(self, csq):
        if csq['Feature'] not in self.rest_cache:
            try:
                self.rest_cache[csq['Feature']] = self._get_rest_data(csq)
            except Exception as err:
                self.logger.warn(err)
                self.logger.warn("REST lookups for {} failed".format(
                                                               csq['Feature']))
                self.rest_cache[csq['Feature']] = [''] * 5
        return self.rest_cache[csq['Feature']]

    def _get_rest_data(self, csq):
        entrez = ''
        go = ''
        reactome = ''
        traits = ''
        mim = ''
        go_data = self.ensembl_rest.get_xref(csq['Feature'], external_db='GO')
        if go_data:
            go = str.join("|", (x['description'] for x in go_data if 
                                x['description'] is not None))
        if csq['Gene']:
            xref_data = self.ensembl_rest.get_xref(csq['Gene'])
            
            entrez = str.join("|", (x['primary_id'] for x in xref_data 
                                    if x['dbname'] == 'EntrezGene'))
            reactome = str.join("|", (x['description'] for x in xref_data 
                                      if x['dbname'] == 'Reactome_gene'))
            mim = str.join("|", (x['description'] for x in xref_data 
                                 if x['dbname'] == 'MIM_MORBID'))
            orth = self.ensembl_rest.lookup_ortholog(csq['Gene'])
            if orth is not None:
                traits = str.join("|", self.ensembl_rest.get_traits(orth))
        return [entrez, go, reactome, traits, mim]

    def write_row(self, worksheet, row, values):
        ''' Write a list of values to given worksheet and row '''
        col = 0
        for x in values:
            worksheet.write(row, col, x)
            col += 1
        return col

    def write_rest_data(self, worksheet, row, col, csq):
        entrez, *rest = (self.get_ensembl_rest_data(csq))
        if entrez:
            worksheet.write_url(row, col, 'https://www.ncbi.nlm.nih.gov/gene/'+
                                entrez, string=entrez)
        else:
            worksheet.write(row, col, entrez)
        col += 1
        for x in rest:
            worksheet.write(row, col, x)
            col += 1

    def write_records(self, record, family, inheritance, allele, features):
        for csq in (x for x in record.CSQ if x['Feature'] in features and 
                    x['alt_index'] == allele):
            #column order is: Inheritance, vcf_output_columns, allele, AC, AN, 
            #                 GTS, VEP fields
            values = [inheritance]
            values.extend(getattr(record, f) for f in vcf_output_columns)
            values.append(allele)
            values.append(record.INFO_FIELDS['AC'])
            values.append(record.INFO_FIELDS['AN'])
            values.extend(record.CALLS[x] for x in 
                          self._get_sample_order(family))
            values.extend(csq[x] for x in self.vcf.header.csq_fields if 
                              x != 'Allele')
            col = self.write_row(self.worksheets[family], self.rows[family], 
                                 values)
            if self.rest_lookups:
                self.write_rest_data(self.worksheets[family],  
                                     self.rows[family],
                                     col, 
                                     csq)
            self.rows[family] += 1

    def write_report(self):
        ''' Write a line for every segregating allele in VCF to XLSX'''
        self.logger.info("Reading variants and writing report")
        n = 0
        w = 0
        for record in self.vcf.parser:
            n += 1
            info = record.parsed_info_fields(list(self.seg_fields.keys()) + 
                                             list(feat_annots.values()))
            for annot,pattern in self.seg_fields.items():
               if annot in info:
                    for i in range(len(info[annot])): 
                        if info[annot][i] is None: 
                            continue
                        if self.all_features:
                            feat = list(x['Feature'] for x in record.CSQ if 
                                        x['Feature'] != '' and 
                                        x['alt_index'] == i +1)
                        else:
                            feat = info[feat_annots[annot]][i].split("|")
                        for fam in info[annot][i].split("|"):
                            if fam in self.families:
                                self.write_records(record, fam, pattern, i+1, 
                                                   feat)
            if n % self.prog_interval == 0:
                self.logger.info("Parsed {:,} records".format(n))
        self.logger.info("Finished parsing {:,} records".format(n))
        self.workbook.close()
