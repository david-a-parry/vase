import csv
import gzip
import numpy as np
import struct
from collections import defaultdict


def csv_to_dict(f, index, fieldnames, delimiter=',', keys_are_unique=False):
    if keys_are_unique:
        d = dict()
    else:
        d = defaultdict(list)
    method = open
    if f.endswith(".gz") or f.endswith(".bgz"):
        method = gzip.open
    with method(f, 'rt', newline='', errors='replace') as csvfile:
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


def read_tbi(tbi):
    '''Return a dict with .tbi index binning information.'''
    ridx = dict()
    with gzip.open(tbi, 'rb') as f:
        magic = f.read(4)
        if magic != b'TBI\x01':
            raise ValueError('Invalid index - wrong magic number ' +
                             '({}) for {}'.format(magic, tbi))
        header = np.frombuffer(f.read(4 * 8), dtype=np.int32)
        names = f.read(header[7]).split(b'\x00')
        for i in range(len(names)):
            bindx = dict()
            for j in range(struct.unpack('<i', f.read(4))[0]):  # n_bins
                bin_key = struct.unpack('<I', f.read(4))[0]  # bin
                n_chunk = struct.unpack('<i', f.read(4))[0]  # n_chunk
                bindx[bin_key] = np.frombuffer(f.read(8 * 2 * n_chunk),
                                               dtype=np.uint64).reshape(
                                                   n_chunk, -1)
            d = {'bindx': bindx}
            n_intv = struct.unpack('<i', f.read(4))[0]
            d['n_intv'] = n_intv
            d['ioff'] = np.frombuffer(f.read(8 * n_intv), dtype=np.uint64)
            ridx[names[i].decode()] = d
    return ridx


def reg2bins(begin, end, min_shift=14, depth=5):
    '''Calculate possible bins for tbi/csi index based retrieval'''
    t, s = 0, min_shift + (depth << 1) + depth
    for l in range(depth + 1):
        b, e = t + (begin >> s), t + (end >> s)
        n = e - b + 1
        for k in range(b, e + 1):
            yield k
            n += 1
        t += 1 << ((l << 1) + l)
        s -= 3
