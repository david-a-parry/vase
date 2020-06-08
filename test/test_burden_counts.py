from .utils import *


def test_burden_counts():
    output = get_tmp_out(suffix='.txt')
    test_args = dict(
        cases=["Sample1", "Sample2"],
        controls=["Sample3"],
        burden_counts=output,
        csq=["default"],
        output='/dev/null',
    )
    run_args(test_args)
    expected_results = os.path.join(dir_path,
                                    "test_data",
                                    "expected_outputs",
                                    "test_burden_counts.txt")
    with open(output, 'rt') as infile:
        results = [x for x in infile.read().split("\n") if x != '']
    with open(expected_results, 'rt') as infile:
        expected = [x for x in infile.read().split("\n") if x != '']
    assert_equal(len(results), len(expected))
    assert_equal(set(results), set(expected))
    os.remove(output)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
