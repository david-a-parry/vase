from .utils import *


def test_g2p():
    output = get_tmp_out()
    input = os.path.join(dir_path, 'test_data', 'ex2.bcf')
    test_args = dict(
        no_warnings=True,
        input=input,
        output=output,
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        de_novo=True,
        biallelic=True,
        csq=['default'],
        check_g2p_consequence=True,
        check_g2p_inheritance=True,
        g2p=os.path.join(dir_path, "test_data", "test_g2p.csv")
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


if __name__ == '__main__':
    test_g2p()
    import nose
    nose.run(defaultTest=__name__)
