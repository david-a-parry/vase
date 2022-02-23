from .utils import *

splice_ai_vcf = os.path.join(dir_path, "test_data", "splice_ai_scores.vcf.gz")
prescored_vcf = os.path.join(dir_path, "test_data", "splice_ai_prescored.vcf.gz")
snpeff_vcf = input_prefix + '.snpeff.vcf.gz'
cadd_file = os.path.join(dir_path, "test_data", "test_cadd_scores.tsv.gz")


def teardown_module():
    for vcf in [splice_ai_vcf, prescored_vcf]:
        idx = vcf + '.tbi'
        if os.path.exists(idx):
            os.remove(idx)


def get_info_annotations(anno_vcf, annot):
    expected_annots = dict()
    with pysam.VariantFile(anno_vcf) as vcf:
        for record in vcf:
            rid = var_string_from_record(record)
            expected_annots[rid] = record.info[annot]
    return expected_annots


def test_annotate_prescored():
    output = get_tmp_out()
    test_args = dict(
        splice_ai_vcfs=[prescored_vcf],
        output=output,
    )
    run_args(test_args)
    annot = 'SpliceAI'
    expected = get_info_annotations(splice_ai_vcf, annot)
    hits = 0
    with pysam.VariantFile(output) as vcf:
        for record in vcf:
            rid = var_string_from_record(record)
            if rid in expected:
                assert_equal(record.info[annot], expected[rid])
                hits += 1
            else:
                assert(annot not in record.info)
    assert_equal(hits, len(expected))
    os.remove(output)


def test_annotate_splice_ai():
    output = get_tmp_out()
    test_args = dict(
        splice_ai_vcfs=[splice_ai_vcf],
        output=output,
    )
    run_args(test_args)
    annot = 'SpliceAI'
    expected = get_info_annotations(splice_ai_vcf, annot)
    hits = 0
    with pysam.VariantFile(output) as vcf:
        for record in vcf:
            rid = var_string_from_record(record)
            if rid in expected:
                assert_equal(record.info[annot], expected[rid])
                hits += 1
            else:
                assert(annot not in record.info)
    assert_equal(hits, len(expected))
    os.remove(output)


def test_filter_splice_ai():
    output = get_tmp_out()
    test_args = dict(
        splice_ai_vcfs=[splice_ai_vcf],
        output=output,
        splice_ai_min_delta=0.5,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_filter_splice_ai_prescored():
    output = get_tmp_out()
    test_args = dict(
        splice_ai_vcfs=[prescored_vcf],
        output=output,
        splice_ai_min_delta=0.5,
    )
    results, expected = run_args(test_args, output, "test_filter_splice_ai")
    assert_equal(results, expected)
    os.remove(output)


def test_filter_splice_ai_multi():
    output = get_tmp_out()
    test_args = dict(
        splice_ai_vcfs=[splice_ai_vcf, prescored_vcf],
        output=output,
        splice_ai_min_delta=0.5,
    )
    results, expected = run_args(test_args, output, "test_filter_splice_ai")
    assert_equal(results, expected)
    os.remove(output)


def test_filter_splice_ai_vep():
    output = get_tmp_out()
    test_args = dict(
        splice_ai_vcfs=[splice_ai_vcf],
        output=output,
        splice_ai_min_delta=0.5,
        impact=["HIGH"],
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_filter_splice_ai_snpeff():
    output = get_tmp_out()
    test_args = dict(
        input=snpeff_vcf,
        snpeff=True,
        splice_ai_vcfs=[splice_ai_vcf],
        output=output,
        splice_ai_min_delta=0.5,
        impact=["HIGH"],
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_filter_splice_ai_cadd():
    output = get_tmp_out()
    test_args = dict(
        splice_ai_vcfs=[splice_ai_vcf],
        output=output,
        splice_ai_min_delta=0.5,
        cadd_phred=20,
        cadd_files=[cadd_file],
        prioritise_cadd=True,
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


def test_filter_splice_ai_vep_cadd():
    output = get_tmp_out()
    test_args = dict(
        splice_ai_vcfs=[splice_ai_vcf],
        output=output,
        splice_ai_min_delta=0.5,
        impact=["HIGH"],
        cadd_phred=20,
        cadd_files=[cadd_file],
    )
    results, expected = run_args(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
