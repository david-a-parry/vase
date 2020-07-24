import sys
import os
import tempfile
import gzip
import pysam
import numpy as np
from nose.tools import *
from vase.vase_runner import VaseRunner
from argparse import Namespace

dir_path = os.path.dirname(os.path.realpath(__file__))
input_prefix = os.path.join(dir_path, 'test_data', 'ex1')
vcf_filter = os.path.join(dir_path,
                          "test_data",
                          "vcf_filter_test.vcf.gz,test_vcf")
gnomad = os.path.join(dir_path, "test_data", "gnomadTest.vcf.gz")
dbsnp = os.path.join(dir_path, "test_data", "dbSnpTest.vcf.gz")

all_args = {
    'input': input_prefix + '.vcf',
    'output': None,
    'report_prefix': None,
    'burden_counts': None,
    'gnomad_burden': False,
    'variant_quality': None,
    'pass_filters': False,
    'keep_filters': None,
    'exclude_filters': None,
    'var_types': None,
    'max_alt_alleles': None,
    'filter_asterisk_only_calls': False,
    'af': None,
    'min_af': None,
    'filtering_an': 0,
    'min_an': 0,
    'ac': None,
    'min_ac': None,
    'info_filters': None,
    'csq': None,
    'impact': None,
    'canonical': False,
    'flagged_features': False,
    'biotypes': [],
    'feature_blacklist': None,
    'loftee': False,
    'missense_filters': [],
    'filter_unpredicted': False,
    'keep_if_any_damaging': False,
    'splice_filters': None,
    'splice_filter_unpredicted': False,
    'splice_keep_if_any_damaging': False,
    'retain_labels': None,
    'no_vep_freq': False,
    'vep_af': [],
    'pathogenic': False,
    'no_conflicted': False,
    'g2p': None,
    'check_g2p_consequence': False,
    'check_g2p_inheritance': False,
    'region': None,
    'bed': None,
    'gene_bed': None,
    'stream': False,
    'exclude_regions': False,
    'cadd_files': [],
    'cadd_directory': None,
    'missing_cadd_scores': None,
    'cadd_phred': None,
    'cadd_raw': None,
    'dbsnp': [],
    'gnomad': [],
    'gnomad_pops': ['AFR', 'AMR', 'EAS', 'FIN', 'NFE', 'SAS'],
    'vcf_filter': [],
    'dng_vcf': None,
    'freq': None,
    'min_freq': None,
    'max_gnomad_homozygotes': None,
    'build': None,
    'max_build': None,
    'filter_known': False,
    'filter_novel': False,
    'clinvar_path': False,
    'ignore_existing_annotations': False,
    'splice_ai_vcfs': [],
    'splice_ai_min_delta': None,
    'splice_ai_max_delta': None,
    'missing_splice_ai_scores': None,
    'cases': [],
    'controls': [],
    'ped': None,
    'gq': 20,
    'dp': 0,
    'max_dp': 0,
    'het_ab': 0.0,
    'hom_ab': 0.0,
    'min_support': None,
    'control_gq': None,
    'control_dp': None,
    'control_max_dp': None,
    'control_het_ab': None,
    'control_hom_ab': None,
    'control_max_ref_ab': None,
    'mosaic_ab': 0.01,
    'mosaic_control_ab': 0.001,
    'sv_gq': 20,
    'sv_dp': 0,
    'sv_max_dp': 0,
    'sv_het_ab': 0.0,
    'sv_hom_ab': 0.0,
    'sv_control_gq': None,
    'sv_control_dp': None,
    'sv_control_max_dp': None,
    'sv_control_het_ab': None,
    'sv_control_hom_ab': None,
    'sv_control_max_ref_ab': None,
    'duphold_del_dhffc': None,
    'duphold_dup_dhbfc': None,
    'control_duphold_del_dhffc': None,
    'control_duphold_dup_dhbfc': None,
    'n_cases': None,
    'n_controls': None,
    'confirm_control_gts': False,
    'biallelic': False,
    'de_novo': False,
    'mosaic': False,
    'dominant': False,
    'min_families': 1,
    'singleton_recessive': [],
    'singleton_dominant': [],
    'seg_controls': [],
    'strict_recessive': False,
    'prog_interval': 1000,
    'log_progress': False,
    'no_progress': True,
    'quiet': True,
    'debug': False,
    'no_warnings': False,
    'silent': False
    }


def get_tmp_out(suffix='.vcf'):
    f, fname = tempfile.mkstemp(suffix=suffix)
    return fname


def get_args(test_args):
    args = all_args.copy()
    args.update(test_args)
    return Namespace(**args)


def get_expected_out(f):
    f += ".txt"
    f = os.path.join(dir_path, "test_data", "expected_outputs", f)
    with open(f, 'rt') as infile:
        lines = [x for x in infile.read().split("\n") if x != '']
    return lines


def var_string_from_record(record):
    return "{}:{}-{}/{}".format(record.chrom,
                                record.pos,
                                record.ref,
                                ','.join(record.alts))


def convert_results(f):
    with pysam.VariantFile(f) as vcf:
        return [var_string_from_record(x) for x in vcf]


def info_fields_from_vcf(f, i):
    results = []
    with pysam.VariantFile(f) as vcf:
        for record in vcf:
            results.extend(record.info[i] if i in record.info else None)
    return results


def run_args(args, output=None, func_name=None):
    args = get_args(args)
    runner = VaseRunner(args)
    runner.run()
    if output is not None and func_name is not None:
        results = convert_results(output)
        expected = get_expected_out(func_name)
        return (results, expected)
