from .utils import *


def test_read_bcf():
    output = get_tmp_out()
    test_args = dict(
        input=input_prefix + '.bcf',
        output=output,
    )
    run_args(test_args)
    expected = convert_results(input_prefix + '.bcf')
    results = convert_results(output)
    assert_equal(results, expected)
    os.remove(output)


def test_read_bgz():
    output = get_tmp_out()
    test_args = dict(
        input=input_prefix + '.vcf.gz',
        output=output,
    )
    run_args(test_args)
    expected = convert_results(input_prefix + '.vcf.gz')
    results = convert_results(output)
    assert_equal(results, expected)
    os.remove(output)


def test_read_vcf():
    output = get_tmp_out()
    test_args = dict(
        input=input_prefix + '.vcf',
        output=output,
    )
    run_args(test_args)
    expected = convert_results(input_prefix + '.vcf')
    results = convert_results(output)
    assert_equal(results, expected)
    os.remove(output)


def test_write_bcf():
    output = get_tmp_out(suffix='.bcf')
    test_args = dict(
        output=output,
    )
    run_args(test_args)
    expected = convert_results(input_prefix + '.vcf')
    results = convert_results(output)
    assert_equal(results, expected)
    with pysam.VariantFile(output) as vcf:
        assert(vcf.is_bcf)
        assert(vcf.compression == 'BGZF')
    os.remove(output)


def test_write_bgz():
    output = get_tmp_out(suffix='.vcf.gz')
    test_args = dict(
        output=output,
    )
    run_args(test_args)
    expected = convert_results(input_prefix + '.vcf')
    results = convert_results(output)
    assert_equal(results, expected)
    with pysam.VariantFile(output) as vcf:
        assert(vcf.is_vcf)
        assert(vcf.compression == 'BGZF')
    os.remove(output)


def test_write_vcf():
    output = get_tmp_out(suffix='.vcf')
    test_args = dict(
        output=output,
    )
    run_args(test_args)
    expected = convert_results(input_prefix + '.vcf')
    results = convert_results(output)
    assert_equal(results, expected)
    with pysam.VariantFile(output) as vcf:
        assert(vcf.is_vcf)
        assert(vcf.compression == 'NONE')
    os.remove(output)



def test_pass_filters():
    output = get_tmp_out()
    test_args = dict(
        pass_filters=True,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_alts_filters():
    output = get_tmp_out()
    test_args = dict(
        max_alt_alleles=1,
        keep_filters=['VQSRTrancheSNP99.90to99.95'],
        output=output,
    )
    results, expected = run_args(test_args, output,
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
    results, expected = run_args(test_args, output,
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
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
