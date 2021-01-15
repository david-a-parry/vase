import sys
import io
import re
import logging
import pysam
import xlsxwriter
import os
import json
from collections import OrderedDict, defaultdict
from .ped_file import PedFile, Individual, PedError
from .vcf_reader import VcfReader
from .ensembl_rest_queries import EnsemblRestQueries
from .utils import csv_to_dict
from .g2p import G2P, allelic_req_to_label

ENST = re.compile(r'''^ENS\w*T\d{11}(\.\d+)?''')
ENTREZ_RE = re.compile(r'''(\d+)(\|(\d+))*''')
SUPPORTED_OUTPUT = ['json', 'xlsx']

vcf_output_columns = ['chrom', 'pos', 'id', 'ref', 'alts', 'qual', 'filter']

MG_FIELDS = ['entrezgene', 'name', 'summary', 'go', 'MIM', 'generif']

feat_annots = {'VASE_biallelic_families': 'VASE_biallelic_features',
               'VASE_dominant_families': 'VASE_dominant_features',
               'VASE_de_novo_families': 'VASE_de_novo_features'}

impact_order = dict((k, n) for n, k in enumerate(['HIGH', 'MODERATE', 'LOW',
                                                  'MODIFIER']))

vep_keys = {'biotype': 'BIOTYPE',
            'impact': 'IMPACT',
            'symbol': 'SYMBOL',
            'gene': 'Gene',
            'feature': 'Feature'}

snpeff_keys = {'biotype': 'Transcript_Biotype',
               'impact': 'Annotation_Impact',
               'symbol': 'Gene_Name',
               'gene': 'Gene_ID',
               'feature': 'Feature_ID'}


class VaseReporter(object):
    '''
        Read a VASE annotated VCF and output XLSX or JSON summary of
        segregating variants.
    '''

    def __init__(self, vcf, out, ped=None, singletons=[], families=[],
                 output_type='xlsx', all_features=False, rest_lookups=False,
                 grch37=False, g2p=None, blacklist=None,
                 recessive_only=False, dominant_only=False, de_novo_only=False,
                 filter_non_g2p=False, allelic_requirement=False,
                 mutation_requirement=False, mygene_lookups=False,
                 info_fields=[], gnomad_constraint=None,
                 choose_transcript=False, prog_interval=None, timeout=2.0,
                 max_retries=2, quiet=False, debug=False, force=False,
                 hide_empty=False, custom_feat_annots=None):
        self._set_logger(quiet, debug)
        self.output_type = output_type.lower()
        if self.output_type not in SUPPORTED_OUTPUT:
            raise ValueError("Unsupported output type: {}\n".format(
                self.output_type) + "Supported types are:\n" + "\n\t".join(
                    SUPPORTED_OUTPUT))
        if type(vcf) == VcfReader:
            self.vcf = vcf
        else:
            self.vcf = VcfReader(vcf)
        if ped:
            self.ped = PedFile(ped)
            if singletons:
                self._add_singletons(singletons)
        elif singletons:
            p_string = ''
            for s in singletons:
                p_string += str.join("\t", (s, s, "0", "0", "0", "2")) + "\n"
            self.ped = PedFile(io.StringIO(p_string))
        self.families = families
        self.all_features = all_features
        self.recessive_only = recessive_only
        self.dominant_only = dominant_only
        self.de_novo_only = de_novo_only
        self.choose_transcript = choose_transcript
        self.seg_fields = self._get_seg_fields()
        self.info_annotations = self._get_info_fields(info_fields)
        self.csq_fields, self._snpeff_mode = self._get_csq_fields()
        self._csq_keys = snpeff_keys if self._snpeff_mode else vep_keys
        if not out.endswith(".{}".format(output_type)):
            out = out + ".{}".format(output_type)
        if os.path.exists(out) and not force:
            sys.exit("Output file '{}' already exists - ".format(out) +
                     "choose another name or use --force to overwrite.")
        self.out = out
        self.sample_orders = dict()  # key is fam id, value = list of samples
        self._header_columns = dict()  # key is fam id, value = list of columns
        self.out_fh = self._get_output_handle()
        self.rest_lookups = rest_lookups
        self.require_g2p = filter_non_g2p
        self.g2p = None
        self.allelic_requirement = allelic_requirement
        self.mutation_requirement = mutation_requirement
        if g2p:
            self.g2p = G2P(g2p, snpeff_mode=self._snpeff_mode)
        elif self.require_g2p:
            raise RuntimeError("--filter_non_g2p option requires a G2P " +
                               "CSV to be supplied with the --g2p argument")
        self.constraint = None
        if gnomad_constraint:
            self.constraint = self._read_gnomad_constraint(gnomad_constraint)
        self.biotype_order = self._get_biotype_order()
        self.blacklist = None
        if blacklist:
            self.blacklist = self._read_blacklist(blacklist)
        if not families:
            families = sorted(self.ped.families.keys())
        self._fam_order = list()
        for f in families:  # respect fam order provided b4 converting to set
            if f not in self.ped.families:
                raise RuntimeError("Family '{}' not found in ped ".format(f) +
                                   "'{}'.".format(ped))
            if f in self._fam_order:
                self.logger.warn("Duplicate family '{}' specified".format(f))
            else:
                self._fam_order.append(f)
        self.families = set(self._fam_order)
        if self.rest_lookups:
            self.ensembl_rest = EnsemblRestQueries(use_grch37_server=grch37,
                                                   timeout=timeout,
                                                   max_retries=max_retries,
                                                   log_level=self.logger.level)
        self.rest_cache = dict()
        self.mygene_cache = dict()
        self.mygene_lookups = False
        if mygene_lookups:
            try:
                import mygene
                self.mg = mygene.MyGeneInfo()
                self.mygene_lookups = True
            except ModuleNotFoundError:
                self.logger.warn("Error importing mygene - please install " +
                                 "mygene  (e.g. pip3 install mygene) to use " +
                                 "the --mygene_lookups option.")
                self.logger.warn("Continuing without mygene lookups")
        if custom_feat_annots:
            self.feat_annots = custom_feat_annots
        else:
            self.feat_annots = feat_annots
        if prog_interval is not None:
            self.prog_interval = prog_interval
        elif self.rest_lookups or self.mygene_lookups:
            self.prog_interval = 100
        else:
            self.prog_interval = 1000
        self.hide_empty = hide_empty
        if self.output_type == 'xlsx':
            self._intialize_workbook()
        elif self.output_type == 'json':
            self.json_dict = defaultdict(list)
        else:
            raise RuntimeError("Unsupported output: {}".format(
                self.output_type))

    def _finish_up(self):
        if self.output_type == 'xlsx':
            if self.hide_empty:
                for family in self.worksheets.keys():
                    if self.rows[family] < 2:
                        self.worksheets[family].hide()
        elif self.output_type == 'json':
            json.dump(self.json_dict, self.out_fh, indent=2,)
        self.out_fh.close()

    def _intialize_workbook(self):
        self.bold = self.out_fh.add_format({'bold': True})
        self.worksheets = dict()
        self.rows = dict()
        for f in self._fam_order:
            w = self._initialize_worksheet(f)
            if w is not None:
                self.worksheets[f] = w
                self.rows[f] = 1

    def _get_output_handle(self):
        '''
            Returns xlsxwriter.Workbook or plain filehandle depending on
            value of self.output_type.
        '''
        if self.output_type == 'xlsx':
            return xlsxwriter.Workbook(self.out)
        else:
            return open(self.out, 'wt')

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

    def _get_info_fields(self, info_fields):
        idict = OrderedDict()
        for x in info_fields:
            if x not in self.vcf.header.info:
                raise ValueError("Requested INFO annotation {} ".format(x) +
                                 "not found in VCF header.")
            idict[x] = self.vcf.header.info[x]
        return idict

    def _add_info_annotations(self, record, allele):
        annots = []
        for inf, d in self.info_annotations.items():
            if d.number == 'A' or d.number == 'R':
                i = allele
                if d.number == 'A':
                    i -= 1
                try:
                    annots.append(record.info[inf][i])
                except KeyError:
                    annots.append('.')
            else:
                try:
                    if d.type == 'Flag':
                        annots.append(inf in record.info)
                    elif d.number == 1:
                        annots.append(record.info[inf])
                    else:
                        annots.append(','.join(str(x) for x in
                                               record.info[inf]))
                except KeyError:
                    annots.append('.')
        return annots

    def _add_singletons(self, singletons):
        for s in singletons:
            try:
                indv = Individual(s, s, 0, 0, 0, 2)
                self.ped.add_individual(indv)
            except PedError:
                raise RuntimeError("Sample '{}' ".format(s) + "specified" +
                                   " by --singletons already exists " +
                                   "in PED file {}" .format(self.ped.filename))

    def _get_csq_fields(self):
        fields = None
        snpeff_mode = False
        try:
            fields = self.vcf.header.csq_fields
            self.logger.info("Found VEP annotations")
        except KeyError:
            try:
                fields = self.vcf.header.ann_fields
                self.logger.info("Found SnpEff annotations")
                snpeff_mode = True
            except KeyError:
                raise KeyError("Neither VEP nor SnpEff annotations found in" +
                               "input file {}".format(self.vcf.filename))
        return fields, snpeff_mode

    def _get_seg_fields(self):
        inheritance_fields = dict()
        seg = False
        if 'VASE_biallelic_families' in self.vcf.header.info:
            inheritance_fields['VASE_biallelic_families'] = 'recessive'
            self.logger.info("Found VASE biallelic annotations.")
            seg = True
        if 'VASE_dominant_families' in self.vcf.header.info:
            inheritance_fields['VASE_dominant_families'] = 'dominant'
            self.logger.info("Found VASE dominant annotations.")
            seg = True
        if 'VASE_de_novo_families' in self.vcf.header.info:
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
        worksheet = self.out_fh.add_worksheet(sheet_name)
        for i in range(len(header)):
            worksheet.write(0, i, header[i], self.bold)
        return worksheet

    def _get_header_columns(self, family=None):
        if family in self._header_columns:
            return self._header_columns[family]
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
        if 'CADD_PHRED_score' in self.vcf.header.header.info:
            header.append("CADD_PHRED_score")
        for inf in self.info_annotations.keys():
            header.append(inf)
        header.extend(x for x in self.csq_fields if x != 'Allele')
        if self.rest_lookups:
            header.extend(["ENTREZ", "Full_Name", "GO", "REACTOME",
                           "MOUSE_TRAITS", "MIM_MORBID"])
        if self.mygene_lookups:
            header.extend(["ENTREZ_ID", "Name", "Summary", "GO_BP", "GO_CC",
                           "GO_MF", "MIM", "GeneRIFs"])
        if self.g2p:
            header.extend(["G2P_disease", "G2P_Category",
                           "G2P_Allelic_Requirement", "G2P_consequences",
                           "G2P_organs"])
        if self.constraint:
            header.extend(["pLI", "pRec", "pNull", "mis_z", "syn_z",
                           "constraint_issues"])
        self._header_columns[family] = header  # None is key for default header
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

    def _read_gnomad_constraint(self, constraint_file):
        self.constraint_cols = ["pLI", "pRec", "pNull", "mis_z", "syn_z",
                                "gene_issues"]
        required_fields = ["gene", "transcript", "canonical", "mis_z", "syn_z",
                           "pLI", "pRec", "pNull", "gene_issues"]
        try:
            d = csv_to_dict(constraint_file, 'transcript', required_fields,
                            delimiter='\t', keys_are_unique=True)
        except ValueError:
            required_fields.pop()
            required_fields.append("constraint_flag")
            d = csv_to_dict(constraint_file, 'transcript', required_fields,
                            delimiter='\t', keys_are_unique=True)
            self.constraint_cols.pop()
            self.constraint_cols.append("constraint_flag")
        # in case our transcript ref differs, also index on gene name
        gene_d = dict()
        for v in d.values():
            if v["gene"] not in d and v["canonical"] == 'true':
                gene_d[v["gene"]] = v
        d.update(gene_d)
        return d

    def _read_blacklist(self, blacklist):
        '''
            Reads a file of feature IDs to ignore. Only reads first
            non-whitespace text from each line.
        '''
        if blacklist is None:
            return None
        with open(blacklist, 'rt') as bfile:
            features = set(line.split()[0] for line in bfile)
        self.logger.info("{:,} unique feature IDs blacklisted from {}.".format(
                         len(features), blacklist))
        return features

    def get_ensembl_rest_data(self, csq):
        if not csq[self._csq_keys['feature']]:
            return [''] * 5
        if csq[self._csq_keys['feature']] not in self.rest_cache:
            if not ENST.match(csq[self._csq_keys['feature']]):
                self.logger.warn("Skipping REST lookup of non-Ensembl " +
                                 "transcript feature " +
                                 "'{}'".format(csq[self._csq_keys['feature']]))
                self.rest_cache[csq[self._csq_keys['feature']]] = [''] * 5
            else:
                try:
                    self.rest_cache[csq[self._csq_keys['feature']]] = \
                        self._get_rest_data(csq)
                except Exception as err:
                    self.logger.warn(err)
                    self.logger.warn("REST lookups for {} failed".format(
                        csq[self._csq_keys['feature']]))
                    self.rest_cache[csq[self._csq_keys['feature']]] = [''] * 5
        return self.rest_cache[csq[self._csq_keys['feature']]]

    def _get_rest_data(self, csq):
        entrez = ''
        full_name = ''
        go = ''
        reactome = ''
        traits = ''
        mim = ''
        try:
            go_data = self.ensembl_rest.get_xref(
                csq[self._csq_keys['feature']], external_db='GO')
            if go_data:
                go = str.join("|", (x['description'] for x in go_data if
                                    x['description'] is not None))
        except Exception as err:
            self.logger.warn(err)
            self.logger.warn("GO lookup for {} failed".format(
                csq[self._csq_keys['feature']]))
            go = 'LOOKUP FAILED'
        if csq[self._csq_keys['gene']]:
            try:
                xref_data = self.ensembl_rest.get_xref(
                    csq[self._csq_keys['gene']])
                entrez = str.join("|", (x['primary_id'] for x in xref_data
                                        if x['dbname'] == 'EntrezGene'))
                full_name = str.join("|", (x['description'] for x in xref_data
                                           if x['dbname'] == 'EntrezGene'))
                reactome = str.join("|", (x['description'] for x in xref_data
                                          if x['dbname'] == 'Reactome_gene'))
                mim = str.join("|", (x['description'] for x in xref_data
                                     if x['dbname'] == 'MIM_MORBID'))
            except Exception as err:
                self.logger.warn(err)
                self.logger.warn("XREF lookups for {} failed".format(
                    csq[self._csq_keys['feature']]))
                entrez = 'LOOKUP FAILED'
                full_name = 'LOOKUP FAILED'
                reactome = 'LOOKUP FAILED'
                mim = 'LOOKUP FAILED'
            try:
                orth = self.ensembl_rest.lookup_ortholog(
                    csq[self._csq_keys['gene']])
                if orth is not None:
                    traits = str.join("|", self.ensembl_rest.get_traits(orth))
            except Exception as err:
                self.logger.warn(err)
                self.logger.warn("Orthology lookup for {} failed".format(
                    csq[self._csq_keys['feature']]))
                traits = 'LOOKUP FAILED'
        return [entrez, full_name, go, reactome, traits, mim]

    def get_mygene_data(self, csq):
        if csq[self._csq_keys['gene']]:
            if csq[self._csq_keys['gene']] in self.mygene_cache:
                return self.mygene_cache[csq[self._csq_keys['gene']]]
            results = self.mg.query(csq[self._csq_keys['gene']],
                                    scopes='ensemblgene,entrezgene',
                                    species=9606, fields=",".join(MG_FIELDS))
            if len(results['hits']) == 0:
                results = self.mg.query(csq[self._csq_keys['symbol']],
                                        scopes='symbol',
                                        species=9606,
                                        fields=",".join(MG_FIELDS))
            if len(results['hits']) == 0:
                self.logger.warn("No MyGene hits for gene {}/{}".format(
                    csq[self._csq_keys['gene']],
                    csq[self._csq_keys['symbol']]))
                self.mygene_cache[csq[self._csq_keys['gene']]] = [''] * 8
                return [''] * 8
            elif len(results['hits']) > 1:
                self.logger.warn("Multiple ({}) ".format(
                                 len(results['hits'])) +
                                 "MyGene hits for gene {}/{}".format(
                                     csq[self._csq_keys['gene']],
                                     csq[self._csq_keys['symbol']]) +
                                 " - will use first hit only.")
            data = []
            for field in MG_FIELDS:
                if field not in results['hits'][0]:
                    if field == 'go':
                        data.extend(["Not found"] * 3)
                    else:
                        data.append("Not found")
                elif field == 'go':
                    for subgo in ['BP', 'CC', 'MF']:
                        if subgo not in results['hits'][0]['go']:
                            data.append("Not found")
                        elif isinstance(results['hits'][0]['go'][subgo], dict):
                            data.append(
                                results['hits'][0]['go'][subgo]['term'])
                        elif isinstance(results['hits'][0]['go'][subgo], list):
                            data.append("|".join([x['term'] for x in
                                        results['hits'][0]['go'][subgo]]))
                        else:
                            data.append("Not found")
                elif field == 'generif':
                    data.append("|".join([x['text'] for x in
                                          results['hits'][0]['generif']]))
                else:
                    data.append(results['hits'][0][field])
                self.mygene_cache[csq[self._csq_keys['gene']]] = data
            return data
        return [''] * 8

    def write_row(self, worksheet, row, values, family):
        ''' Write a list of values to given worksheet and row '''
        entrez_cols = []
        mim_col = -1
        if "ENTREZ" in self._header_columns[family]:  # from Ensembl REST
            entrez_cols.append(self._header_columns[family].index("ENTREZ"))
        if "ENTREZ_ID" in self._header_columns[family]:  # from MyGene
            entrez_cols.append(self._header_columns[family].index("ENTREZ_ID"))
        if "MIM" in self._header_columns[family]:  # from MyGene
            mim_col = self._header_columns[family].index("MIM")
        for col in range(len(values)):
            if col in (entrez_cols):  # provide link to Entrez Gene
                m = ENTREZ_RE.match(values[col])
                if m:  # if multiple ENTREZ ids pick the most recent/largest
                    eid = max([int(x) for x in m.groups() if x and
                               "|" not in x])
                    worksheet.write_url(row, col,
                                        'https://www.ncbi.nlm.nih.gov/gene/' +
                                        str(eid), string=values[col])
                else:
                    worksheet.write(row, col, values[col])
            elif col == mim_col:  # provide link to OMIM
                worksheet.write_url(row, col,
                                    'https://omim.org/entry/' +
                                    values[col], string=values[col])
            else:
                worksheet.write(row, col, values[col])
        return col

    def get_g2p_data(self, csq):
        d = [''] * 5
        if csq[self._csq_keys['symbol']] in self.g2p.g2p:
            d = ["|".join(x[f] for x in
                 self.g2p.g2p[csq[self._csq_keys['symbol']]]) for f in
                 ['disease name', 'DDD category', 'allelic requirement',
                 'mutation consequence', 'organ specificity list']]
        return d

    def get_constraint_data(self, csq):
        cons = [''] * 5
        k = None
        if csq[self._csq_keys['feature']] in self.constraint:
            k = csq[self._csq_keys['feature']]
        elif csq[self._csq_keys['symbol']] in self.constraint:
            k = csq[self._csq_keys['symbol']]
        if k is not None:
            cons = [self.constraint[k][f] for f in self.constraint_cols]
        return cons

    def write_records(self, record, csqs, family, inheritance, allele,
                      features):
        for csq in (x for x in csqs if x[self._csq_keys['feature']] in features
                    and x['alt_index'] == allele):
            # column order is: Inheritance, vcf_output_columns, allele, AC, AN,
            #                  GTS, VEP fields
            if self.require_g2p:
                if csq[self._csq_keys['symbol']] not in self.g2p.g2p:
                    continue
                if self.allelic_requirement:
                    inh_ok = False
                    for req in (x['allelic requirement'] for x in
                                self.g2p.g2p[csq[self._csq_keys['symbol']]]):
                        for r in req.split(","):
                            if inheritance in allelic_req_to_label[r]:
                                inh_ok = True
                                break
                    if not inh_ok:
                        continue
                if self.mutation_requirement:
                    if not self.g2p.csq_matches_requirement(
                            csq, keep_uncertain=True):
                        continue
            if self.blacklist:
                if csq[self._csq_keys['feature']] in self.blacklist:
                    continue
            values = [inheritance]
            values.extend(",".join(str(y) for y in x) if isinstance(x, tuple)
                          or isinstance(x, pysam.libcbcf.VariantRecordFilter)
                          else str(x) if x is not None else '.' for x in
                          (getattr(record, f) for f in vcf_output_columns))
            values.append(allele)
            for i_field in (('AC', 'AN')):
                if i_field in record.info:
                    if isinstance(record.info[i_field], tuple):
                        values.append(",".join(str(x) for x in
                                               record.info[i_field]))
                    else:
                        values.append(record.info[i_field])
                else:
                    values.append('.')
            values.append(":".join(record.format))
            for s in self._get_sample_order(family):
                values.append(":".join(
                    ("/".join(str(y) for y in v) if k == 'GT'
                     else (",".join(str(y) for y in v) if isinstance(v, tuple)
                           else str(v)) if v is not None else '.'
                     for k, v in record.samples[s].items())))
            if 'CADD_PHRED_score' in self.vcf.header.header.info:
                try:
                    values.append(record.info['CADD_PHRED_score'][allele-1])
                except KeyError:
                    values.append('.')
            values.extend(self._add_info_annotations(record, allele))
            values.extend(csq[x] for x in self.csq_fields if x != 'Allele')
            if self.rest_lookups:
                values.extend(self.get_ensembl_rest_data(csq))
            if self.mygene_lookups:
                values.extend(self.get_mygene_data(csq))
            if self.g2p:
                values.extend(self.get_g2p_data(csq))
            if self.constraint:
                values.extend(self.get_constraint_data(csq))
            if self.output_type == 'xlsx':
                _ = self.write_row(self.worksheets[family],
                                   self.rows[family],
                                   values,
                                   family)
                self.rows[family] += 1
            elif self.output_type == 'json':
                jrow = dict((k, v) for k, v in zip(
                    self._get_header_columns(family), values))
                self.json_dict[family].append(jrow)
            else:
                self.out_fh.write("\t".join([family] + values) + "\n")

    def pick_transcript(self, features, allele, csq):
        '''
            Return a single feature from the list provided, picking
            transcripts with the highest impact and preferring biotypes
            as ordered in vase/data/biotypes.tsv and canonical
            transcripts where impacts are the same.
        '''
        feat_csq = (x for x in csq if x[self._csq_keys['feature']] in features
                    and x['alt_index'] == allele)
        if 'CANONICAL' in csq[0]:
            sorted_csq = sorted(feat_csq, key=lambda x: x['CANONICAL'],
                                reverse=True)
        sorted_csq = sorted(sorted_csq, key=lambda x:
                            self.biotype_order[x[self._csq_keys['biotype']]])
        sorted_csq = sorted(sorted_csq, key=lambda x:
                            impact_order[x[self._csq_keys['impact']]])
        if sorted_csq:
            return sorted_csq[0][self._csq_keys['feature']]
        # none of the features are in the CSQ annotations (e.g. feature is a
        # genomic coordinate)
        return features[0]

    def write_report(self):
        '''Write an entry for every segregating allele in VCF'''
        self.logger.info("Reading variants and writing report")
        n = 0
        for record in self.vcf:
            n += 1
            if self._snpeff_mode:
                csqs = record.ANN
            else:
                csqs = record.CSQ
            for annot, pattern in self.seg_fields.items():
                if annot in record.info:
                    for i in range(len(record.info[annot])):
                        if record.info[annot][i] is None:
                            continue
                        if self.all_features:
                            feat = list(x[self._feature_key] for x in csqs if
                                        x[self._feature_key] != '' and
                                        x['alt_index'] == i + 1)
                        else:
                            feat = \
                             record.info[self.feat_annots[annot]][i].split("|")
                        if self.choose_transcript:
                            feat = [self.pick_transcript(feat, i+1,
                                                         csqs)]
                        for fam in record.info[annot][i].split("|"):
                            if fam in self.families:
                                self.write_records(record, csqs, fam, pattern,
                                                   i+1, feat)
            if n % self.prog_interval == 0:
                self.logger.info("Parsed {:,} records".format(n))
        self.logger.info("Finished parsing {:,} records".format(n))
        self._finish_up()
