import sys
import os
import tempfile
import pysam
from nose.tools import *
from vase.vase_runner import VaseRunner
from argparse import Namespace

dir_path = os.path.dirname(os.path.realpath(__file__))
input_prefix = os.path.join(dir_path, 'test_data', 'ex1')
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
    'control_gq': None,
    'control_dp': None,
    'control_max_dp': None,
    'control_het_ab': None,
    'control_hom_ab': None,
    'control_max_ref_ab': None,
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


def get_tmp_out():
    f, fname = tempfile.mkstemp(suffix='.vcf')
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


def run_test(args, output, func_name):
    args = get_args(args)
    runner = VaseRunner(args)
    runner.run()
    results = convert_results(output)
    expected = get_expected_out(func_name)
    return (results, expected)


def test_alts_filters():
    output = get_tmp_out()
    test_args = dict(
        max_alt_alleles=1,
        keep_filters=['VQSRTrancheSNP99.90to99.95'],
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_info_filters():
    output = get_tmp_out()
    test_args = dict(
        af=0.4,
        ac=5,
        info_filters=['FS < 1', 'DB = True'],
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_case_control():
    output = get_tmp_out()
    test_args = dict(
        cases=['Sample3', 'Sample2'],
        controls=['Sample1'],
        het_ab=0.005,
        gq=20,
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_csq():
    output = get_tmp_out()
    test_args = dict(
        csq=[],
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_impact():
    output = get_tmp_out()
    test_args = dict(
        impact=["HIGH"],
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_vartype():
    output = get_tmp_out()
    test_args = dict(
        var_types=["INDEL"],
        max_alt_alleles=1,
        output=output,
        pass_filters=True,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_de_novo():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_name, "test_data", "test.ped"),
        de_novo=True,
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_de_novo2():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_name, "test_data", "test.ped"),
        de_novo=True,
        max_alt_alleles=1,
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_de_novo3():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_name, "test_data", "test.ped"),
        de_novo=True,
        het_ab=0.25,
        max_alt_alleles=1,
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_biallelic():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_name, "test_data", "test.ped"),
        biallelic=True,
        csq=True,
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_biallelic2():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_name, "test_data", "test.ped"),
        biallelic=True,
        impact="HIGH",
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_biallelic3():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_name, "test_data", "test.ped"),
        biallelic=True,
        impact="HIGH",
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_dominant():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_name, "test_data", "test2.ped"),
        dominant=True,
        csq=True,
        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest = __name__)
