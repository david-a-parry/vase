from .utils import *

vcf_filter = os.path.join(dir_path,
                          "test_data",
                          "vcf_filter_test.vcf.gz")
gnomad = os.path.join(dir_path, "test_data", "gnomadTest.vcf.gz")
dbsnp = os.path.join(dir_path, "test_data", "dbSnpTest.vcf.gz")


def teardown_module():
    for f in [vcf_filter, gnomad, dbsnp]:
        idx = f + '.tbi'
        if os.path.exists(idx):
            os.remove(idx)


def test_vcf_filter_known():
    output = get_tmp_out()
    test_args = dict(
        vcf_filter=[vcf_filter + ',test_vcf'],
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
        vcf_filter=[vcf_filter + ',test_vcf'],
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
        vcf_filter=[vcf_filter + ',test_vcf'],
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
