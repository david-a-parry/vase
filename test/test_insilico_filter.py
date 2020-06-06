from .utils import *

input = os.path.join(dir_path, 'test_data', 'ex5.bcf')


def test_insilico_freq():
    output = get_tmp_out()
    test_args = dict(
        input=input,
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
        input=input,
        output=output,
        csq=[],
        missense_filters=['sift=deleterious'],
    )
    results, expected = run_args(test_args, output,
                                 'test_insilico_pred_sift')
    assert_equal(results, expected)
    os.remove(output)
    output = get_tmp_out()
    test_args = dict(
        input=input,
        output=output,
        csq=[],
        missense_filters=['Polyphen'])
    results, expected = run_args(test_args, output,
                                 'test_insilico_pred_polyphen')
    assert_equal(results, expected)
    output = get_tmp_out()
    test_args = dict(
        input=input,
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
        input=input,
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
        input=input,
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
        input=input,
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
        input=input,
        output=output,
        csq=[],
        keep_if_any_damaging=True,
        missense_filters=['Polyphen', 'sift=deleterious'],
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
