from .utils import *

input = os.path.join(dir_path, 'test_data', 'ex4.bcf')


def test_sv_biallelic_csq():
    output = get_tmp_out()
    test_args = dict(
        input=input,
        max_alt_alleles=1,
        output=output,
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        biallelic=True,
        csq=[],
        sv_gq=99,
        sv_het_ab=0.3,
        duphold_del_dhffc=0.7,
        duphold_dup_dhbfc=1.3,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_sv_biallelic_lof():
    output = get_tmp_out()
    test_args = dict(
        input=input,
        max_alt_alleles=1,
        output=output,
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        biallelic=True,
        impact=['HIGH'],
        sv_gq=99,
        sv_het_ab=0.3,
        duphold_del_dhffc=0.7,
        duphold_dup_dhbfc=1.3,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_sv_de_novo():
    output = get_tmp_out()
    test_args = dict(
        input=input,
        max_alt_alleles=1,
        output=output,
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        de_novo=True
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_sv_de_novo_ad():
    output = get_tmp_out()
    test_args = dict(
        input=input,
        max_alt_alleles=1,
        output=output,
        ad=12,
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        de_novo=True
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_sv_de_novo_filters():
    output = get_tmp_out()
    test_args = dict(
        input=input,
        max_alt_alleles=1,
        output=output,
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        csq=[],
        sv_gq=99,
        sv_het_ab=0.3,
        duphold_del_dhffc=0.7,
        duphold_dup_dhbfc=1.3,
        de_novo=True
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
