from .utils import *

bed = os.path.join(dir_path, 'test_data', 'test_regions.bed')


def test_fail_on_uncompressed():
    test_args = dict(
        region=['1:1060742-1061726', '1:1083580-1084363'],
        output='/dev/null'
    )
    assert_raises(TypeError, run_args, test_args)


def test_var_by_region_bcf():
    var_by_region('.bcf')


def test_var_by_region_vcf():
    var_by_region('.vcf.gz')


def test_var_by_region_stream():
    var_by_region('.vcf', True)


def test_var_by_region_stream_compressed():
    var_by_region('.vcf.gz', True)


def test_var_from_bed():
    output = get_tmp_out()
    test_args = dict(
        input=input_prefix + '.vcf.gz',
        bed=bed,
        output=output,
    )
    results, expected = run_args(test_args, output, 'test_var_by_region')
    assert_equal(results, expected)
    os.remove(output)


def test_var_from_gene_bed():
    output = get_tmp_out()
    test_args = dict(
        input=input_prefix + '.vcf.gz',
        gene_bed=bed,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    not_expected = get_expected_out('test_var_by_region')
    assert_not_equal(results, not_expected)
    os.remove(output)


def var_by_region(suffix, stream=False):
    output = get_tmp_out()
    test_args = dict(
        input=input_prefix + suffix,
        region=['1:1060742-1061726', '1:1083580-1084363'],
        output=output,
        stream=stream,
    )
    results, expected = run_args(test_args, output, 'test_var_by_region')
    assert_equal(results, expected)
    os.remove(output)
