from .utils import *
from vase.vase_reporter import VaseReporter
import xlrd
import json
from collections import defaultdict

rep_input = os.path.join(dir_path, 'test_data', 'ex7.bcf')
ped = os.path.join(dir_path, "test_data", "test4.ped")


def _get_xlsx_output(xlsx, fam):
    wb = xlrd.open_workbook(xlsx)
    sheet = wb.sheet_by_name(fam)
    values = []
    for i in range(sheet.nrows):
        values.append([sheet.cell_value(i, j) for j in range(sheet.ncols)])
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
    vr  = VaseReporter(
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
    vr  = VaseReporter(
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


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
