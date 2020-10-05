from .utils import *

fb_input = os.path.join(dir_path, 'test_data', 'ex1.fb.vcf')
mosaic_input = os.path.join(dir_path, 'test_data', 'ex9.bcf')

def test_case_control():
    output = get_tmp_out()
    test_args = dict(
        cases=['Sample3', 'Sample2'],
        controls=['Sample1'],
        het_ab=0.005,
        gq=20,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_case_control_fb():
    output = get_tmp_out()
    test_args = dict(
        input=fb_input,
        cases=['Sample3', 'Sample2'],
        controls=['Sample1'],
        het_ab=0.005,
        gq=20,
        output=output,
    )
    results, expected = run_args(test_args, output, "test_case_control")
    assert_equal(results, expected)
    os.remove(output)


def test_de_novo():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        de_novo=True,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_de_novo2():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        de_novo=True,
        max_alt_alleles=1,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_de_novo3():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        de_novo=True,
        het_ab=0.25,
        control_max_ref_ab=0.05,
        max_alt_alleles=1,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_de_novo_fb():
    output = get_tmp_out()
    test_args = dict(
        input=fb_input,
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        de_novo=True,
        het_ab=0.25,
        control_max_ref_ab=0.05,
        max_alt_alleles=1,
        output=output,
    )
    results, expected = run_args(test_args, output, "test_de_novo3")
    assert_equal(results, expected)
    os.remove(output)


def test_de_novo_no_csq():
    output = get_tmp_out()
    test_args = dict(
        input=os.path.join(dir_path, 'test_data', 'ex9.vcf.gz'),
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        de_novo=True,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_biallelic():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        biallelic=True,
        csq=[],
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_biallelic2():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        biallelic=True,
        impact=['HIGH', 'MODERATE'],
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_biallelic3():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        biallelic=True,
        impact=['HIGH'],
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_biallelic_no_ped():
    output = get_tmp_out()
    test_args = dict(
        singleton_recessive=['Sample1'],
        seg_controls=['Sample2', 'Sample3'],
        csq=[],
        output=output,
    )
    results, expected = run_args(test_args, output, "test_biallelic")
    assert_equal(results, expected)
    os.remove(output)


def test_biallelic_seg_control():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_path, "test_data", "test3.ped"),
        singleton_recessive=['Sample1'],
        seg_controls=['Sample2', 'Sample3'],
        csq=[],
        output=output,
    )
    assert_raises(ValueError, run_args, test_args)
    test_args = dict(
        ped=os.path.join(dir_path, "test_data", "test3.ped"),
        biallelic=True,
        seg_controls=['Sample2', 'Sample3'],
        csq=[],
        output=output,
    )
    results, expected = run_args(test_args, output, "test_biallelic")
    assert_equal(results, expected)
    os.remove(output)


def test_dominant():
    output = get_tmp_out()
    test_args = dict(
        ped=os.path.join(dir_path, "test_data", "test2.ped"),
        dominant=True,
        csq=[],
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_mosaic():
    output = get_tmp_out()
    test_args = dict(
        input=mosaic_input,
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        mosaic=True,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_mosaic2():
    output = get_tmp_out()
    test_args = dict(
        input=mosaic_input,
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        mosaic=True,
        mosaic_control_ab=8e-4,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_mosaic3():
    output = get_tmp_out()
    test_args = dict(
        input=mosaic_input,
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        mosaic=True,
        min_support=2,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
