#!/usr/bin/env python3
import sys
import re


arg_re = re.compile(r'--\S+')


def output_args(args, func):
    print('''
def {}():
    output = get_tmp_out()
    test_args = dict('''.format(func))
    for x in args:
        print("        {}={},".format(x[0].replace("-", ""), x[1]))
    print('''        output=output,
    )
    results, expected = run_test(test_args, output,
                                 sys._getframe().f_code.co_name)
    assert_equal(results, expected)
    os.remove(output)
''')


def main(f):
    with open(f, 'rt') as infile:
        for line in infile:
            if line.rstrip() == '':
                continue
            args = []
            s = line.split()
            for i in range(len(s) - 1):
                if arg_re.match(s[i]):
                    if arg_re.match(s[i + 1]) or s[i + 1] == '|':
                        args.append((s[i], True))
                    else:
                        try:
                            x = int(s[i + 1])
                        except ValueError:
                            try:
                                x = float(s[i + 1])
                            except ValueError:
                                x = '"{}"'.format(s[i + 1])
                        args.append((s[i], x))
            func = s[-1].split('/')[-1].replace('.txt', '')
            output_args(args, func)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} tests.sh\n".format(sys.argv[0]))
    main(sys.argv[1])
