from .utils import *

vep_and_snpeff_inputs = [(input_prefix + '.vcf.gz', False),
                         (input_prefix + '.snpeff.vcf.gz', True)]


def test_case_control():
    for vcf, snpeff in vep_and_snpeff_inputs:
        output = get_tmp_out()
        test_args = dict(
            input=vcf,
            snpeff=snpeff,
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


def test_de_novo():
    for vcf, snpeff in vep_and_snpeff_inputs:
        output = get_tmp_out()
        test_args = dict(
            input=vcf,
            snpeff=snpeff,
            ped=os.path.join(dir_path, "test_data", "test.ped"),
            de_novo=True,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                    sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


def test_de_novo2():
    for vcf, snpeff in vep_and_snpeff_inputs:
        output = get_tmp_out()
        test_args = dict(
            input=vcf,
            snpeff=snpeff,
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
    for vcf, snpeff in vep_and_snpeff_inputs:
        output = get_tmp_out()
        test_args = dict(
            input=vcf,
            snpeff=snpeff,
            ped=os.path.join(dir_path, "test_data", "test.ped"),
            de_novo=True,
            het_ab=0.25,
            max_alt_alleles=1,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                    sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


def test_de_novo4():
    for vcf, snpeff in vep_and_snpeff_inputs:
        output = get_tmp_out()
        test_args = dict(
            input=vcf,
            snpeff=snpeff,
            ped=os.path.join(dir_path, "test_data", "test.ped"),
            de_novo=True,
            ad=4,
            output=output,
        )
        results, expected = run_args(test_args, output,
                                     sys._getframe().f_code.co_name)
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



def test_de_novo_ref_ad():
    output = get_tmp_out()
    test_args = dict(
        input=os.path.join(dir_path, 'test_data', 'ex9.vcf.gz'),
        ped=os.path.join(dir_path, "test_data", "test.ped"),
        de_novo=True,
        control_max_ref_ad=2,
        output=output,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_biallelic():
    for vcf, snpeff in vep_and_snpeff_inputs:
        output = get_tmp_out()
        test_args = dict(
            input=vcf,
            snpeff=snpeff,
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
    for vcf, snpeff in vep_and_snpeff_inputs:
        output = get_tmp_out()
        test_args = dict(
            input=vcf,
            snpeff=snpeff,
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
    for vcf, snpeff in vep_and_snpeff_inputs:
        output = get_tmp_out()
        test_args = dict(
            input=vcf,
            snpeff=snpeff,
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
    for vcf, snpeff in vep_and_snpeff_inputs:
        output = get_tmp_out()
        test_args = dict(
            input=vcf,
            snpeff=snpeff,
            singleton_recessive=['Sample1'],
            seg_controls=['Sample2', 'Sample3'],
            csq=[],
            output=output,
        )
        results, expected = run_args(test_args, output, "test_biallelic")
        assert_equal(results, expected)
        os.remove(output)


def test_biallelic_seg_control():
    for vcf, snpeff in vep_and_snpeff_inputs:
        output = get_tmp_out()
        test_args = dict(
            input=vcf,
            snpeff=snpeff,
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
    for vcf, snpeff in vep_and_snpeff_inputs:
        output = get_tmp_out()
        test_args = dict(
            input=vcf,
            snpeff=snpeff,
            ped=os.path.join(dir_path, "test_data", "test2.ped"),
            dominant=True,
            csq=[],
            output=output,
        )
        results, expected = run_args(test_args, output,
                                    sys._getframe().f_code.co_name)
        assert_equal(results, expected)
        os.remove(output)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
