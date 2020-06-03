from .utils import *


def test_vcf_filter_known():
    output = get_tmp_out()
    test_args = dict(
        vcf_filter=[vcf_filter],
        filter_known=True,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_vcf_filter_novel():
    output = get_tmp_out()
    test_args = dict(
        vcf_filter=[vcf_filter],
        filter_novel=True,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_vcf_filter_freq():
    output = get_tmp_out()
    test_args = dict(
        vcf_filter=[vcf_filter],
        freq=0.1,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_dbsnp_known():
    output = get_tmp_out()
    test_args = dict(
        dbsnp=[dbsnp],
        filter_known=True,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_dbsnp_novel():
    output = get_tmp_out()
    test_args = dict(
        dbsnp=[dbsnp],
        filter_novel=True,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_dbsnp_freq():
    output = get_tmp_out()
    test_args = dict(
        dbsnp=[dbsnp],
        freq=0.005,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_gnomad_novel():
    output = get_tmp_out()
    test_args = dict(
        gnomad=[gnomad],
        filter_novel=True,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_gnomad_known():
    output = get_tmp_out()
    test_args = dict(
        gnomad=[gnomad],
        filter_known=True,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_gnomad_freq():
    output = get_tmp_out()
    test_args = dict(
        gnomad=[gnomad],
        freq=0.0005,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
