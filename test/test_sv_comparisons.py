from .utils import *
from vase.vcf_reader import VcfReader
import unittest
input = os.path.join(dir_path, 'test_data', 'ex6.bcf')

class SvComparisons(unittest.TestCase):

    @classmethod
    def setup_class(self):
        vcf = VcfReader(input)
        self.first_del = next(vcf).DECOMPOSED_ALLELES[0]
        self.same_del = next(vcf).DECOMPOSED_ALLELES[0]
        self.similar_del = next(vcf).DECOMPOSED_ALLELES[0]
        self.diff_del_end = next(vcf).DECOMPOSED_ALLELES[0]
        self.diff_del_start = next(vcf).DECOMPOSED_ALLELES[0]
        self.diff_dup = next(vcf).DECOMPOSED_ALLELES[0]
        self.first_ins = next(vcf).DECOMPOSED_ALLELES[0]
        self.same_ins = next(vcf).DECOMPOSED_ALLELES[0]
        self.diff_ins =  next(vcf).DECOMPOSED_ALLELES[0]
        self.first_bnd =  next(vcf).DECOMPOSED_ALLELES[0]
        self.same_bnd =  next(vcf).DECOMPOSED_ALLELES[0]
        self.diff_bnd =  next(vcf).DECOMPOSED_ALLELES[0]

    def test_identical_dels(self):
        self.first_del.breakpoint_precision = 0.1
        assert_equal(self.first_del, self.same_del)
        self.first_del.breakpoint_precision = 0.0
        assert_equal(self.first_del, self.same_del)

    def test_similar_dels(self):
        self.first_del.breakpoint_precision = 0.0
        assert_not_equal(self.first_del, self.similar_del)
        self.first_del.breakpoint_precision = 0.1
        assert_equal(self.first_del, self.similar_del)

    def test_different_dels(self):
        self.first_del.breakpoint_precision = 0.1
        assert_not_equal(self.first_del, self.diff_del_end)
        assert_not_equal(self.first_del, self.diff_del_start)

    def test_different_sv_types(self):
        self.first_del.breakpoint_precision = 0.2
        assert_not_equal(self.first_del, self.diff_dup)

    def test_same_ins(self):
        assert_equal(self.first_ins, self.same_ins)

    def test_diff_ins(self):
        assert_not_equal(self.first_ins, self.diff_ins)

    def test_same_bnd(self):
        assert_equal(self.first_bnd, self.same_bnd)

    def test_diff_bnd(self):
        assert_not_equal(self.first_bnd, self.diff_bnd)
