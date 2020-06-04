from .utils import *
from vase.gt_annotator import GtAnnotator
from collections import defaultdict

def get_annotated_gts(anno_vcf, annot):
    expected_annots = defaultdict(dict)
    with pysam.VariantFile(anno_vcf) as vcf:
        for record in vcf:
            rid = var_string_from_record(record)
            for s in record.samples:
                expected_annots[rid][s] = record.samples[s][annot]
    return expected_annots


def check_annots(test_vcf, annot, expected_annots):
    with pysam.VariantFile(test_vcf) as vcf:
        n = vcf.header.formats[annot].number
        for record in vcf:
            rid = var_string_from_record(record)
            for s in record.samples:
                if rid in expected_annots:
                    sys.stderr.write("Found {} = {} for {} at {}\n".format(
                        annot,
                        expected_annots[rid][s],
                        s,
                        rid)
                    )
                    if n == 1:
                        assert_almost_equal(record.samples[s][annot],
                                            expected_annots[rid][s],
                                            5)
                    else:
                        r = np.array(record.samples[s][annot], dtype=float)
                        e = np.array(expected_annots[rid][s], dtype=float)
                        np.testing.assert_array_almost_equal(r, e, decimal=4)
                else:
                    assert(annot not in record.samples[s])


def test_gt_annot():
    anno_vcf = os.path.join(dir_path, 'test_data', 'ex3.bcf')
    gtanno = (GtAnnotator(anno_vcf, ['RND']))
    with pysam.VariantFile(input_prefix + '.bcf') as vcf:
        vcf.header.formats.add("RND", "G", "Float", "Random numbers")
        output = pysam.VariantFile(get_tmp_out(), mode='w', header=vcf.header)
        for record in vcf:
            gtanno.annotate(record)
            output.write(record)
    output.close()
    expected_annots = get_annotated_gts(anno_vcf, 'RND')
    check_annots(output.filename, 'RND', expected_annots)
    os.remove(output.filename)


def test_dng_annot():
    output = get_tmp_out()
    anno_vcf = os.path.join(dir_path, 'test_data', 'dng_test.bcf')
    test_args = dict(
        output=output,
        dng_vcf=[anno_vcf]
    )
    run_args(test_args)
    for pp in ['PP_DNM', 'PP_NULL']:
        expected_annots = get_annotated_gts(anno_vcf, pp)
        check_annots(output, pp, expected_annots)
    os.remove(output)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
