from .utils import *

cadd_file = os.path.join(dir_path, "test_data", "test_cadd_scores.tsv.gz")


def teardown_module():
    idx = cadd_file + '.tbi'
    if os.path.exists(idx):
        os.remove(idx)


def test_cadd_annot():
    output = get_tmp_out()
    missing_cadd = get_tmp_out(suffix='.vcf')
    test_args = dict(
        output=output,
        missing_cadd_scores=missing_cadd,
        cadd_files=[cadd_file],
    )
    run_args(test_args)
    expected_results = os.path.join(dir_path,
                                    "test_data",
                                    "expected_outputs",
                                    "test_cadd_annot.vcf.gz")
    for i in ['CADD_PHRED_score', 'CADD_raw_score']:
        results = np.array(info_fields_from_vcf(output, i), dtype=float)
        expected = np.array(info_fields_from_vcf(expected_results, i),
                            dtype=float)
        np.testing.assert_array_almost_equal(results, expected, decimal=5)
    expected_missing = os.path.join(dir_path,
                                    "test_data",
                                    "expected_outputs",
                                    "test_cadd_missing.vcf.gz")
    with gzip.open(missing_cadd + '.gz', 'rt') as infile:
        results = [x for x in infile.read().split("\n") if x != '']
    with gzip.open(expected_missing, 'rt') as infile:
        expected = [x for x in infile.read().split("\n") if x != '']
    assert_equal(results, expected)
    for f in (output, missing_cadd):
        os.remove(f)


def test_cadd_phred_filters():
    output = get_tmp_out()
    test_args = dict(
        output=output,
        cadd_files=[cadd_file],
        cadd_phred=30
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)


def test_cadd_raw_filters():
    output = get_tmp_out()
    test_args = dict(
        output=output,
        cadd_files=[cadd_file],
        cadd_raw=2.25
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)

