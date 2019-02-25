import sys
import re
import logging
import xlsxwriter
import os
import gzip
import csv
from collections import namedtuple
from .ped_file import PedFile, Family, Individual, PedError
from parse_vcf import VcfReader, VcfHeader, VcfRecord
from .ensembl_rest_queries import EnsemblRestQueries

ENST = re.compile(r'''^ENS\w*T\d{11}(\.\d+)?''')
vcf_output_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',]

feat_annots = { 'VASE_biallelic_families': 'VASE_biallelic_features',
                'VASE_dominant_families': 'VASE_dominant_features',
                'VASE_de_novo_families': 'VASE_de_novo_features',}

allelic_req_to_label = {'biallelic'                 : ['recessive'],
                        'digenic'                   : None,
                        'hemizygous'                : ['recessive'],
                        'imprinted'                 : None,
                        'mitochondrial'             : None,
                        'monoallelic'               : ['de novo', 'dominant'],
                        'mosaic'                    : ['de novo', 'dominant'],
                        'x-linked dominant'         : None,
                        'x-linked over-dominance'   : None,
                       }
impact_order = dict((k,n) for n,k in enumerate(['HIGH', 'MODERATE', 'LOW',
                                              'MODIFIER']))

class VaseReporter(object):
    ''' Read a VASE annotated VCF and output XLSX format summary of
        segregating variants. '''

    def __init__(self, vcf, ped, out, families=[], all_features=False,
                 rest_lookups=False, grch37=False, ddg2p=None, blacklist=None,
                 recessive_only=False, dominant_only=False, de_novo_only=False,
                 filter_non_ddg2p=False, allelic_requirement=False,
                 choose_transcript=False, prog_interval=None, timeout=2.0,
                 max_retries=2, quiet=False, debug=False, force=False,
                 hide_empty=False):
        self._set_logger(quiet, debug)
        self.vcf = VcfReader(vcf)
        self.ped = PedFile(ped)
        self.families = families
        self.all_features = all_features
        self.recessive_only = recessive_only
        self.dominant_only = dominant_only
        self.de_novo_only = de_novo_only
        self.choose_transcript = choose_transcript
        self.seg_fields = self._get_seg_fields()
        if not out.endswith(".xlsx"):
            out = out + ".xlsx"
        if os.path.exists(out) and not force:
            sys.exit("Output file '{}' already exists - ".format(out) +
                     "choose another name or use --force to overwrite.")
        self.out = out
        self.sample_orders = dict() #key is fam id, value is list of samples
        self.workbook = xlsxwriter.Workbook(out)
        self.bold = self.workbook.add_format({'bold': True})
        self.worksheets = dict()
        self.rows = dict()
        self.hide_empty = hide_empty
        self.rest_lookups = rest_lookups
        self.require_ddg2p = filter_non_ddg2p
        self.ddg2p = None
        self.allelic_requirement = allelic_requirement
        if ddg2p:
            self.ddg2p = self._read_ddg2p_csv(ddg2p)
        elif self.require_ddg2p:
            raise RuntimeError("--filter_non_ddg2p option requires a DDG2P " +
                               "CSV to be supplied with the --ddg2p argument")
        self.biotype_order = self._get_biotype_order()
        self.blacklist = None
        if blacklist:
            self.blacklist = self._read_blacklist(blacklist)
        if not families:
            families = sorted(self.ped.families.keys())
        fams_with_samples = []
        for f in families:#respect the order of families provided before converting to set
            if f not in self.ped.families:
                raise RuntimeError("Family '{}' not found in ped ".format(f) +
                                   "'{}'.".format(ped))
            if f in self.worksheets:
                logger.warn("Duplicate family '{}' specified".format(f))
            else:
                w = self._initialize_worksheet(f)
                if w is not None:
                    self.worksheets[f] = w
                    self.rows[f] = 1
                    fams_with_samples.append(f)
        self.families = set(fams_with_samples)
        if self.rest_lookups:
            self.ensembl_rest = EnsemblRestQueries(use_grch37_server=grch37,
                                                   timeout=timeout,
                                                   max_retries=max_retries,
                                                   log_level=self.logger.level)
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


    def _get_seg_fields(self):
        inheritance_fields = dict()
        selected_fields = dict()
        seg = False
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
            self.logger.info("Found VASE de novo annotations.")
            seg = True
        if not seg:
            raise RuntimeError("no vase recessive/dominant/de novo " +
                               "annotations found in vcf - please run vase " +
                               "with --biallelic, --dominant or --de_novo " +
                               " options first.")
        return self._select_seg_fields(inheritance_fields)

    def _select_seg_fields(self, seg_fields):
        if (not self.recessive_only and not self.dominant_only and not
            self.de_novo_only):
            return seg_fields
        selected = dict()
        if self.dominant_only:
            self.logger.info("Selecting dominant annotations")
            if 'VASE_dominant_families' in seg_fields:
                selected['VASE_dominant_families'] = 'dominant'
            else:
                raise RuntimeError("no vase dominant annotations found in " +
                                   "vcf - please run vase with --dominant " +
                                   "option first.")
        if self.de_novo_only:
            self.logger.info("Selecting de novo annotations")
            if 'VASE_de_novo_families' in seg_fields:
                selected['VASE_de_novo_families'] = 'de_novo'
            else:
                raise RuntimeError("no vase de_novo annotations found in " +
                                   "vcf - please run vase with --de_novo " +
                                   "option first.")
        if self.recessive_only:
            self.logger.info("Selecting recessive annotations")
            if 'VASE_biallelic_families' in seg_fields:
                selected['VASE_biallelic_families'] = 'recessive'
            else:
                raise RuntimeError("no vase recessive annotations found in " +
                                   "vcf - please run vase with --biallelic " +
                                   "option first.")
        return selected

    def _initialize_worksheet(self, family):
        header = self._get_header_columns(family)
        if header is None:
            return None
        sheet_name = re.sub(r'[\[\]\:\*\?\/]', '_', family)
        worksheet = self.workbook.add_worksheet(sheet_name)
        for i in range(len(header)):
            worksheet.write(0, i, header[i], self.bold)
        return worksheet

    def _get_header_columns(self, family=None):
        header = ['INHERITANCE'] + vcf_output_columns + ['ALLELE', 'AC', 'AN',
                                                         'FORMAT']
        if family is not None:
            samples = self._get_sample_order(family)
            if not samples:
                self.logger.warn("No samples in VCF for family {}".format(
                                                                        family))
                return None
            self.sample_orders[family] = samples
            header.extend(samples)
        if 'CADD_PHRED_score' in self.vcf.header.metadata['INFO']:
            header.append("CADD_PHRED_score")
        header.extend(x for x in self.vcf.header.csq_fields if x != 'Allele')
        if self.rest_lookups:
            header.extend(["ENTREZ", "GO", "REACTOME", "MOUSE_TRAITS",
                           "MIM_MORBID"])
        if self.ddg2p:
            header.extend(["DDG2P_disease", "DDG2P_Category",
                           "DDG2P_Allelic_Requirement", "DDG2P_consequences",
                           "DDG2P_organs"])
        return header


    def _get_sample_order(self, family):
        if family in self.sample_orders:
            return self.sample_orders[family]
        samples = []
        children = []
        for par in sorted(self.ped.families[family].parents):
            if par in self.vcf.header.samples:
                samples.append(par)
            for child in self.ped.families[family].parents[par]:
                if child not in children and child in self.vcf.header.samples:
                    children.append(child)
        for child in children:
            if child not in samples and child in self.vcf.header.samples:
                samples.append(child)
        for indv in self.ped.families[family].individuals:
            if indv not in samples and indv in self.vcf.header.samples:
                samples.append(indv)
        return samples

    def _get_biotype_order(self):
        data_file = os.path.join(os.path.dirname(__file__),
                                 "data",
                                 "biotypes.tsv")
        biotype_order = dict()
        with open(data_file, 'rt') as bfile:
            for line in bfile:
                if line.startswith('#'):
                    continue
                cols = line.rstrip().split('\t')
                if len(cols) < 3:
                    continue
                biotype_order[cols[0]] = int(cols[2])
        return biotype_order

    def _read_ddg2p_csv(self, ddg2p):
        g2p = dict()
        method = open
        if ddg2p.endswith(".gz"):
            method = gzip.open
        with method(ddg2p, 'rt', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            required_fields = ['gene symbol', 'disease name', 'DDD category',
                               'allelic requirement', 'mutation consequence',
                               'organ specificity list', 'prev symbols']
            for f in required_fields:
                if f not in reader.fieldnames:
                    raise RuntimeError("Missing '{}' ".format(f) + "field in" +
                                       "DDG2P file '{}'".format(f, ddg2p))
                for row in reader:
                    g2p[row['gene symbol']] = row
                    if row['prev symbols']:
                        for prev in row['prev symbols'].split(';'):
                            g2p[prev] = row
        return g2p

    def _read_blacklist(self, blacklist):
        '''
            Reads a file of feature IDs to ignore. Only reads first
            non-whitespace text from each line.
        '''
        if blacklist is None:
            return None
        with open (blacklist, 'rt') as bfile:
            features = set(line.split()[0] for line in bfile)
        self.logger.info("{:,} unique feature IDs blacklisted from {}.".format(
                         len(features), blacklist))
        return features

    def get_ensembl_rest_data(self, csq):
        if not csq['Feature']:
            return [''] * 5
        if csq['Feature'] not in self.rest_cache:
            if not ENST.match(csq['Feature']):
                self.logger.warn("Skipping REST lookup of non-Ensembl " +
                                 "transcript feature" +
                                 "'{}'".format(csq['Feature']))
                self.rest_cache[csq['Feature']] = [''] * 5
            else:
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
        try:
            go_data = self.ensembl_rest.get_xref(csq['Feature'],
                                                 external_db='GO')
            if go_data:
                go = str.join("|", (x['description'] for x in go_data if
                                    x['description'] is not None))
        except Exception as err:
            self.logger.warn(err)
            self.logger.warn("GO lookup for {} failed".format(csq['Feature']))
            go = 'LOOKUP FAILED'
        if csq['Gene']:
            try:
                xref_data = self.ensembl_rest.get_xref(csq['Gene'])
                entrez = str.join("|", (x['primary_id'] for x in xref_data
                                        if x['dbname'] == 'EntrezGene'))
                reactome = str.join("|", (x['description'] for x in xref_data
                                          if x['dbname'] == 'Reactome_gene'))
                mim = str.join("|", (x['description'] for x in xref_data
                                     if x['dbname'] == 'MIM_MORBID'))
            except Exception as err:
                self.logger.warn(err)
                self.logger.warn("XREF lookups for {} failed".format(
                                                               csq['Feature']))
                entrez = 'LOOKUP FAILED'
                reactome = 'LOOKUP FAILED'
                mim = 'LOOKUP FAILED'
            try:
                orth = self.ensembl_rest.lookup_ortholog(csq['Gene'])
                if orth is not None:
                    traits = str.join("|", self.ensembl_rest.get_traits(orth))
            except Exception as err:
                self.logger.warn(err)
                self.logger.warn("Orthology lookup for {} failed".format(
                                                               csq['Feature']))
                traits = 'LOOKUP FAILED'
        return [entrez, go, reactome, traits, mim]

    def write_row(self, worksheet, row, values):
        ''' Write a list of values to given worksheet and row '''
        col = 0
        for x in values:
            worksheet.write(row, col, x)
            col += 1
        return col

    def write_ddg2p_data(self, worksheet, row, col, csq):
        d = [''] * 5
        if csq['SYMBOL'] in self.ddg2p:
            d = [self.ddg2p[csq['SYMBOL']][f] for f in
                 ['disease name', 'DDD category', 'allelic requirement',
                  'mutation consequence', 'organ specificity list']]
        for x in d:
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
        return col

    def write_records(self, record, family, inheritance, allele, features):
        for csq in (x for x in record.CSQ if x['Feature'] in features and
                    x['alt_index'] == allele):
            #column order is: Inheritance, vcf_output_columns, allele, AC, AN,
            #                 GTS, VEP fields
            if self.require_ddg2p:
                if csq['SYMBOL'] not in self.ddg2p:
                    continue
                if self.allelic_requirement:
                    req = self.ddg2p[csq['SYMBOL']]['allelic requirement']
                    if req and allelic_req_to_label[req] is not None:
                        if inheritance not in allelic_req_to_label[req]:
                            continue
            if self.blacklist:
                if csq['Feature'] in self.blacklist:
                    continue
            values = [inheritance]
            values.extend(getattr(record, f) for f in vcf_output_columns)
            values.append(allele)
            for i_field in (('AC', 'AN')):
                if i_field in record.INFO_FIELDS:
                    values.append(record.INFO_FIELDS[i_field])
                else:
                    values.append('.')
            values.append(record.FORMAT)
            values.extend(record.CALLS[x] for x in
                          self._get_sample_order(family))
            if 'CADD_PHRED_score' in self.vcf.header.metadata['INFO']:
                try:
                    cinf = record.parsed_info_fields(['CADD_PHRED_score'])
                    values.append(cinf['CADD_PHRED_score'][allele-1])
                except KeyError:
                    values.append('.')
            values.extend(csq[x] for x in self.vcf.header.csq_fields if
                              x != 'Allele')
            col = self.write_row(self.worksheets[family], self.rows[family],
                                 values)
            if self.rest_lookups:
                col = self.write_rest_data(self.worksheets[family],
                                           self.rows[family],
                                           col,
                                           csq)
            if self.ddg2p:
                col = self.write_ddg2p_data(self.worksheets[family],
                                            self.rows[family],
                                            col,
                                            csq)
            self.rows[family] += 1

    def pick_transcript(self, features, allele, csq):
        '''
            Return a single feature from the list provided, picking
            transcripts with the highest impact and preferring biotypes
            as ordered in vase/data/biotypes.tsv and canonical
            transcripts where impacts are the same.
        '''
        feat_csq = (x for x in csq if x['Feature'] in features and
                    x['alt_index'] == allele)
        sorted_csq = sorted(feat_csq, key=lambda x: x['CANONICAL'],
                            reverse=True)
        sorted_csq = sorted(sorted_csq, key=lambda x:
                            self.biotype_order[x['BIOTYPE']])
        sorted_csq = sorted(sorted_csq, key=lambda x:
                            impact_order[x['IMPACT']])
        return sorted_csq[0]['Feature']

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
                        if self.choose_transcript:
                            feat = [self.pick_transcript(feat, i+1,
                                                         record.CSQ)]
                        for fam in info[annot][i].split("|"):
                            if fam in self.families:
                                self.write_records(record, fam, pattern, i+1,
                                                   feat)
            if n % self.prog_interval == 0:
                self.logger.info("Parsed {:,} records".format(n))
        self.logger.info("Finished parsing {:,} records".format(n))
        if self.hide_empty:
            for family in self.worksheets.keys():
                if self.rows[family] < 2:
                    self.worksheets[family].hide()
        self.workbook.close()
