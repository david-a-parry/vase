from .utils import *
from vase.vase_reporter import VaseReporter
from openpyxl import load_workbook
import json
from collections import defaultdict

rep_input = os.path.join(dir_path, 'test_data', 'ex7.bcf')
rep_snpeff_input = os.path.join(dir_path, 'test_data', 'ex7.snpeff.bcf')
ped = os.path.join(dir_path, "test_data", "test4.ped")
g2p = os.path.join(dir_path, "test_data", "test_g2p.csv")


def _get_xlsx_output(xlsx, fam):
    wb = load_workbook(xlsx)
    sheet = wb[fam]
    values = []
    for row in sheet.rows:
        values.append([x.value for x in row])
    return values


def test_reporter_no_overwrite():
    output = get_tmp_out(suffix='.json')
    kwargs = dict(
        ped=ped,
        output_type='json',
        quiet=True,
    )
    assert_raises(SystemExit, VaseReporter, rep_input, output, **kwargs)


def test_reporter_json():
    output = get_tmp_out(suffix='.json')
    vr = VaseReporter(
        vcf=rep_input,
        out=output,
        ped=ped,
        force=True,
        quiet=True,
        output_type='json')
    vr.write_report()
    with open(output, 'rt') as rfile:
        results = json.load(rfile)
    expect_json = os.path.join(
        dir_path, "test_data", "expected_outputs", "test_reporter.json")
    with open(expect_json, 'rt') as efile:
        expected = json.load(efile)
    assert_equal(results, expected)


def test_reporter_xlsx():
    output = get_tmp_out(suffix='.xlsx')
    vr = VaseReporter(
        vcf=rep_input,
        out=output,
        ped=ped,
        force=True,
        quiet=True,
        )
    vr.write_report()
    expect_xlsx = os.path.join(
        dir_path, "test_data", "expected_outputs", "test_reporter.xlsx")
    for f in ['Fam1', 'Fam2']:
        results = _get_xlsx_output(output, f)
        expected = _get_xlsx_output(expect_xlsx, f)
        assert_equal(results, expected)


def test_reporter_info_annotations():
    output = get_tmp_out(suffix='.xlsx')
    vr = VaseReporter(
        vcf=rep_input,
        out=output,
        ped=ped,
        info_fields='MQRankSum MLEAC NEGATIVE_TRAIN_SITE'.split(),
        force=True,
        quiet=True,
        )
    vr.write_report()
    expect_xlsx = os.path.join(dir_path,
                               "test_data",
                               "expected_outputs",
                               sys._getframe().f_code.co_name + ".xlsx")
    for f in ['Fam1', 'Fam2']:
        results = _get_xlsx_output(output, f)
        expected = _get_xlsx_output(expect_xlsx, f)
        assert_equal(results, expected)


def test_reporter_choose_transcript():
    output = get_tmp_out(suffix='.xlsx')
    vr = VaseReporter(
        vcf=rep_input,
        out=output,
        ped=ped,
        choose_transcript=True,
        force=True,
        quiet=True,
        )
    vr.write_report()
    expect_xlsx = os.path.join(dir_path,
                               "test_data",
                               "expected_outputs",
                               sys._getframe().f_code.co_name + ".xlsx")
    for f in ['Fam1', 'Fam2']:
        results = _get_xlsx_output(output, f)
        expected = _get_xlsx_output(expect_xlsx, f)
        assert_equal(results, expected)


def test_reporter_choose_transcript_snpeff():
    output = get_tmp_out(suffix='.xlsx')
    vr = VaseReporter(
        vcf=rep_snpeff_input,
        out=output,
        ped=ped,
        choose_transcript=True,
        force=True,
        quiet=True,
        )
    vr.write_report()
    expect_xlsx = os.path.join(dir_path,
                               "test_data",
                               "expected_outputs",
                               sys._getframe().f_code.co_name + ".xlsx")
    for f in ['Fam1', 'Fam2']:
        results = _get_xlsx_output(output, f)
        expected = _get_xlsx_output(expect_xlsx, f)
        assert_equal(results, expected)


def test_reporter_g2p():
    output = get_tmp_out(suffix='.xlsx')
    vr = VaseReporter(
        vcf=rep_input,
        out=output,
        ped=ped,
        g2p=g2p,
        force=True,
        quiet=True,
        )
    vr.write_report()
    expect_xlsx = os.path.join(dir_path,
                               "test_data",
                               "expected_outputs",
                               sys._getframe().f_code.co_name + ".xlsx")
    for f in ['Fam1', 'Fam2']:
        results = _get_xlsx_output(output, f)
        expected = _get_xlsx_output(expect_xlsx, f)
        assert_equal(results, expected)


def test_reporter_g2p_allelic():
    output = get_tmp_out(suffix='.xlsx')
    vr = VaseReporter(
        vcf=rep_input,
        out=output,
        ped=ped,
        g2p=g2p,
        allelic_requirement=True,
        force=True,
        quiet=True,
        )
    vr.write_report()
    expect_xlsx = os.path.join(dir_path,
                               "test_data",
                               "expected_outputs",
                               sys._getframe().f_code.co_name + ".xlsx")
    for f in ['Fam1', 'Fam2']:
        results = _get_xlsx_output(output, f)
        expected = _get_xlsx_output(expect_xlsx, f)
        assert_equal(results, expected)


def test_reporter_g2p_allelic_snpeff():
    output = get_tmp_out(suffix='.xlsx')
    vr = VaseReporter(
        vcf=rep_snpeff_input,
        out=output,
        ped=ped,
        g2p=g2p,
        allelic_requirement=True,
        force=True,
        quiet=True,
        )
    vr.write_report()
    expect_xlsx = os.path.join(dir_path,
                               "test_data",
                               "expected_outputs",
                               sys._getframe().f_code.co_name + ".xlsx")
    for f in ['Fam1', 'Fam2']:
        results = _get_xlsx_output(output, f)
        expected = _get_xlsx_output(expect_xlsx, f)
        assert_equal(results, expected)


def test_reporter_filter_g2p():
    output = get_tmp_out(suffix='.xlsx')
    vr = VaseReporter(
        vcf=rep_input,
        out=output,
        ped=ped,
        g2p=g2p,
        allelic_requirement=True,
        mutation_requirement=True,
        filter_non_g2p=True,
        force=True,
        quiet=True,
        )
    vr.write_report()
    expect_xlsx = os.path.join(dir_path,
                               "test_data",
                               "expected_outputs",
                               sys._getframe().f_code.co_name + ".xlsx")
    for f in ['Fam1', 'Fam2']:
        results = _get_xlsx_output(output, f)
        expected = _get_xlsx_output(expect_xlsx, f)
        assert_equal(results, expected)
 

def test_reporter_filter_g2p_snpeff():
    output = get_tmp_out(suffix='.xlsx')
    vr = VaseReporter(
        vcf=rep_snpeff_input,
        out=output,
        ped=ped,
        g2p=g2p,
        allelic_requirement=True,
        mutation_requirement=True,
        filter_non_g2p=True,
        force=True,
        quiet=True,
        )
    vr.write_report()
    expect_xlsx = os.path.join(dir_path,
                               "test_data",
                               "expected_outputs",
                               sys._getframe().f_code.co_name + ".xlsx")
    for f in ['Fam1', 'Fam2']:
        results = _get_xlsx_output(output, f)
        expected = _get_xlsx_output(expect_xlsx, f)
        assert_equal(results, expected)


def test_reporter_singleton():
    output = get_tmp_out(suffix='.xlsx')
    vr = VaseReporter(
        vcf=os.path.join(dir_path, 'test_data', 'ex8.bcf'),
        out=output,
        singletons=['Sample4'],
        force=True,
        quiet=True,
        )
    vr.write_report()
    expect_xlsx = os.path.join(dir_path,
                               "test_data",
                               "expected_outputs",
                               sys._getframe().f_code.co_name + ".xlsx")
    results = _get_xlsx_output(output, "Sample4")
    expected = _get_xlsx_output(expect_xlsx, "Sample4")
    assert_equal(results, expected)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
