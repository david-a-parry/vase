from .utils import *

vcf_filter = os.path.join(dir_path,
                          "test_data",
                          "vcf_filter_test.vcf.gz")
gnomad = os.path.join(dir_path, "test_data", "gnomadTest.vcf.gz")
dbsnp = os.path.join(dir_path, "test_data", "dbSnpTest.vcf.gz")


def teardown_module():
    for f in [vcf_filter, gnomad, dbsnp]:
        tbi = f + '.tbi'
        csi = f.replace('.vcf.gz', '.bcf.csi')
        for idx in [tbi, csi]:
            if os.path.exists(idx):
                os.remove(idx)


def test_vcf_filter_known():
    for f in [vcf_filter, vcf_filter.replace('.vcf.gz', '.bcf')]:
        output = get_tmp_out()
        test_args = dict(
            vcf_filter=[f + ',test_vcf'],
            filter_known=True,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                     sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


def test_vcf_filter_novel():
    for f in [vcf_filter, vcf_filter.replace('.vcf.gz', '.bcf')]:
        output = get_tmp_out()
        test_args = dict(
            vcf_filter=[f + ',test_vcf'],
            filter_novel=True,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                    sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


def test_vcf_filter_freq():
    for f in [vcf_filter, vcf_filter.replace('.vcf.gz', '.bcf')]:
        output = get_tmp_out()
        test_args = dict(
            vcf_filter=[f + ',test_vcf'],
            freq=0.1,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                    sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


def test_dbsnp_known():
    for f in [dbsnp, dbsnp.replace('.vcf.gz', '.bcf')]:
        output = get_tmp_out()
        test_args = dict(
            dbsnp=[f],
            filter_known=True,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                     sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


def test_dbsnp_novel():
    for f in [dbsnp, dbsnp.replace('.vcf.gz', '.bcf')]:
        output = get_tmp_out()
        test_args = dict(
            dbsnp=[f],
            filter_novel=True,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                    sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


def test_dbsnp_freq():
    for f in [dbsnp, dbsnp.replace('.vcf.gz', '.bcf')]:
        output = get_tmp_out()
        test_args = dict(
            dbsnp=[f],
            freq=0.005,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                    sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


def test_gnomad_novel():
    for f in [gnomad, gnomad.replace('.vcf.gz', '.bcf')]:
        output = get_tmp_out()
        test_args = dict(
            gnomad=[f],
            filter_novel=True,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                    sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


def test_gnomad_known():
    for f in [gnomad, gnomad.replace('.vcf.gz', '.bcf')]:
        output = get_tmp_out()
        test_args = dict(
            gnomad=[f],
            filter_known=True,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                    sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


def test_gnomad_freq():
    for f in [gnomad, gnomad.replace('.vcf.gz', '.bcf')]:
        output = get_tmp_out()
        test_args = dict(
            gnomad=[f],
            freq=0.0005,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                    sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


def test_gnomad_homozygotes():
    for f in [gnomad, gnomad.replace('.vcf.gz', '.bcf')]:
        output = get_tmp_out()
        test_args = dict(
            gnomad=[f],
            freq=0.0005,
            max_gnomad_homozygotes=0,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                    sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
