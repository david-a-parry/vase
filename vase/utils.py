import csv
from collections import defaultdict
import gzip

def csv_to_dict(f, index, fieldnames, delimiter=',', keys_are_unique=False):
    if keys_are_unique:
        d = dict()
    else:
        d = defaultdict(list)
    method = open
    if f.endswith(".gz") or f.endswith(".bgz"):
        method = gzip.open
    with method(f, 'rt', newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=delimiter)
        for x in fieldnames:
            if x not in reader.fieldnames:
                raise ValueError("Missing '{}' ".format(x) + "field in" +
                                   " file '{}'".format(f))
        for row in reader:
            if keys_are_unique:
                if row[index] in d:
                    raise ValueError("Duplicate value ('{}') ".format(
                        row[index]) + "for '{}' field ".format(index) +
                        "in file '{}'".format(f))
                d[row[index]] = row
            else:
                d[row[index]].append(row)
    return d


