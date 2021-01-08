from .utils import *
from vase.vep_filter import VepFilter
from vase.vcf_reader import VcfReader

is_input = os.path.join(dir_path, 'test_data', 'ex5.bcf')


def test_csq():
    output = get_tmp_out()
    test_args = dict(
        csq=[],
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_impact():
    output = get_tmp_out()
    test_args = dict(
        impact=["HIGH"],
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_insilico_freq():
    output = get_tmp_out()
    test_args = dict(
        input=is_input,
        output=output,
        csq=[],
        freq=0.01,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_insilico_pred():
    output = get_tmp_out()
    test_args = dict(
        input=is_input,
        output=output,
        csq=[],
        biotypes=['default', 'transcribed_processed_pseudogene'],
        missense_filters=['sift=deleterious'],
    )
    results, expected = run_args(test_args, output,
                                 'test_insilico_pred_sift')
    assert_equal(results, expected)
    os.remove(output)
    output = get_tmp_out()
    test_args = dict(
        input=is_input,
        output=output,
        csq=[],
        missense_filters=['Polyphen'])
    results, expected = run_args(test_args, output,
                                 'test_insilico_pred_polyphen')
    assert_equal(results, expected)
    output = get_tmp_out()
    test_args = dict(
        input=is_input,
        output=output,
        csq=[],
        missense_filters=['Polyphen', 'sift=deleterious'],
    )
    results, expected = run_args(test_args, output,
                                 'test_insilico_pred_both')
    assert_equal(results, expected)
    os.remove(output)


def test_insilico_freq_pred():
    output = get_tmp_out()
    test_args = dict(
        input=is_input,
        output=output,
        csq=[],
        freq=0.01,
        missense_filters=['sift=deleterious'],
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_insilico_score():
    output = get_tmp_out()
    test_args = dict(
        input=is_input,
        output=output,
        csq=[],
        missense_filters=['SIFT_score=0.1'])
    results, expected = run_args(test_args, output,
                                 'test_insilico_score')
    assert_equal(results, expected)
    os.remove(output)


def test_insilico_unpredicted():
    output = get_tmp_out()
    test_args = dict(
        input=is_input,
        output=output,
        csq=[],
        filter_unpredicted=True,
        missense_filters=['SIFT_score=0.1'])
    results, expected = run_args(test_args, output,
                                 'test_insilico_unpredicted')
    assert_equal(results, expected)
    os.remove(output)


def test_insilico_keep_any_damaging():
    output = get_tmp_out()
    test_args = dict(
        input=is_input,
        output=output,
        csq=[],
        keep_if_any_damaging=True,
        missense_filters=['Polyphen', 'sift=deleterious'],
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_canonical():
    vcf = VcfReader(is_input)
    csq_filter = VepFilter(
        vcf=vcf,
        csq=[],
        canonical=True)
    results_alts = []
    results_csq = []
    expected_alts = [[False]] * 7
    expected_csq = [[False, True, True, True],
                    [False, True, True, True],
                    [False, True, True, True],
                    [False, False, True, False],
                    [False, True, True, True],
                    [False, True, True, True],
                    [False, True, True, True]]
    for record in vcf:
        r_alts, r_csq = csq_filter.filter(record)
        results_csq.append(r_csq)
        results_alts.append(r_alts)
    assert_equal(expected_csq, results_csq)
    assert_equal(expected_alts, results_alts)


def test_flags():
    vcf = VcfReader(is_input)
    csq_filter = VepFilter(
        vcf=vcf,
        csq=['all'],
        biotypes=['all'],
        filter_flagged_features=True)
    results_alts = []
    results_csq = []
    expected_alts = [[False]] * 7
    expected_csq = [[False, False, False, False],
                    [False, False, False, False],
                    [False, False, False, False],
                    [False, False, False, False],
                    [False, False, False, False],
                    [False, True, False, False],
                    [False, True, False, False]]
    for record in vcf:
        r_alts, r_csq = csq_filter.filter(record)
        results_csq.append(r_csq)
        results_alts.append(r_alts)
    assert_equal(expected_csq, results_csq)
    assert_equal(expected_alts, results_alts)


def test_canonical_stop_gained():
    vcf = VcfReader(is_input)
    csq_filter = VepFilter(
        vcf=vcf,
        csq=['stop_gained'],
        canonical=True)
    results_alts = []
    results_csq = []
    expected_alts = [[True]] * 5 + [[False]] * 2
    expected_csq = [[True, True, True, True],
                    [True, True, True, True],
                    [True, True, True, True],
                    [True, True, True, True],
                    [True, True, True, True],
                    [False, True, True, True],
                    [False, True, True, True]]
    for record in vcf:
        r_alts, r_csq = csq_filter.filter(record)
        results_csq.append(r_csq)
        results_alts.append(r_alts)
    assert_equal(expected_csq, results_csq)
    assert_equal(expected_alts, results_alts)


def test_canonical_lof():
    vcf = VcfReader(is_input)
    csq_filter = VepFilter(
        vcf=vcf,
        csq=['stop_gained'],
        loftee=True,
        canonical=True)
    results_alts = []
    results_csq = []
    expected_alts = [[True]] * 5 + [[False]] + [[True]]
    expected_csq = [[True, True, True, True],
                    [True, True, True, True],
                    [True, True, True, True],
                    [True, True, True, True],
                    [True, True, True, True],
                    [False, True, True, True],
                    [True, True, True, True]]
    for record in vcf:
        r_alts, r_csq = csq_filter.filter(record)
        results_csq.append(r_csq)
        results_alts.append(r_alts)
    assert_equal(expected_csq, results_csq)
    assert_equal(expected_alts, results_alts)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
