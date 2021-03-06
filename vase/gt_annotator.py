from .vcf_reader import VcfReader

_comp_fields = ['chrom', 'pos', 'ref', 'alts']


class GtAnnotator(object):
    '''
        An object for annotating genotype fields of a sample with those
        from another VCF.
    '''

    def __init__(self, vcf, format_fields, samples=None, prefix=None):
        '''
            Initialize object with VCF containing the fields for
            annotating records with and the FORMAT fields to annotate.
            Only matching records will be annotated.

            Args:
                vcf:    VCF file containing variants with annotation
                        fields to add to records. If a record passed to
                        the 'annotate' method has a matching record in
                        this file the specified FORMAT fields will be
                        added to overlapping samples' genotype fields.

                format_fields:
                        One or more FORMAT fields from 'vcf' to add to
                        matching samples and variants.

                samples:
                        One or more samples to annotate. Default is to
                        annotate all overlapping samples.

                prefix:
                        Prefix for added annotations.

        '''
        self.vcf = VcfReader(vcf)
        self.prefix = prefix
        self.format_fields = format_fields
        self.samples = samples
        self._check_args()
        self.header_fields = self._create_header_fields()

    def _create_header_fields(self):
        header_fields = dict()
        for f in self.format_fields:
            desc = ('{} INFO field parsed by {} object from file {}. '.format(
                    f, type(self).__name__, self.vcf.filename) +
                   'Original description was as follows: {}' .format(
                    self.vcf.header.formats[f].description.replace(
                        '"', '')))
            header_fields[f] = {'Number': self.vcf.header.formats[f].number,
                                'Type': self.vcf.header.formats[f].type,
                                'Description': desc}
        return header_fields

    def _check_args(self):
        not_found = [f for f in self.format_fields if f not in
                     self.vcf.header.formats]
        if not_found:
            raise RuntimeError("Could not find '{}' FORMAT field(s) in VCF "
                               .format(",".join(not_found)) +
                               "header for file {}".format(self.vcf))
        if self.samples:
            not_found = [s for s in self.samples if s not in
                         self.vcf.header.samples]
            if not_found:
                raise RuntimeError("Could not find '{}' sample(s) in VCF "
                                   .format(",".join(not_found)) +
                                   "header for file {}".format(
                                       self.vcf.filename))
        else:
            self.samples = list(self.vcf.header.samples)

    def annotate(self, record):
        '''
            Add FORMAT fields to sample calls if a matching record
            exists in self.vcf.
        '''
        match = self.find_matching_record(record)
        if match is None:
            return
        for f in self.format_fields:
            for s in (x for x in match.header.samples if x in self.samples):
                record.samples[s][f] = match.samples[s][f]

    def find_matching_record(self, record):
        overlapping = self.get_overlapping_records(record)
        # require CHROM, POS, REF and ALT fields to be identical
        for x in overlapping:
            if ([getattr(record, f) for f in _comp_fields] ==
                [getattr(x, f) for f in _comp_fields]):
                return x
        return None


    def get_overlapping_records(self, record):
        '''
            For a given record, returns a list of overlapping records
            in the class's VCF.
        '''
        self.vcf.set_region(record.chrom, record.start, record.stop)
        return (s for s in self.vcf)
