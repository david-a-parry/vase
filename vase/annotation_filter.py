import operator

ops = {
    ">": operator.gt,
    "<": operator.lt,
    ">=": operator.ge,
    "<=": operator.le,
    "==": operator.eq,
    "!=": operator.ne,
    "=": operator.eq,
    "!": operator.ne,
}


class AnnotationFilter(object):
    '''
        A class for filtering on given INFO fields in a VCF
    '''

    def __init__(self, vcf, field, filters):
        '''
            Args:
                vcf:    VcfReader object from vase

                field:  VCF field containing annotations. Either "INFO" or
                        "FORMAT".

                filters:
                        iterable of tuples of field names, operands and
                        values for filtering.

        '''
        self.vcf = vcf
        if field.lower() == 'format':
            self.record_field = 'samples'
            self.header_field = 'formats'
            self.field_name = 'FORMAT'
        elif field.lower() == 'info':
            self.record_field = 'info'
            self.header_field = 'info'
            self.field_name = 'INFO'
        else:
            raise ValueError("Unrecognised field: '{}'. ".format(field) +
                             "Only 'INFO' or 'FORMAT' fields are recognised.")
        self.metadata = getattr(self.vcf.header, self.header_field)
        self.fields = set()
        self.filters = self._parse_filters(filters)

    def _parse_filters(self, filters):
        '''
            Args:
                filters:
                    iterable of tuples of field names, operands and
                    values for filtering.
        '''
        f = []
        for field, operand, value in filters:
            try:
                op = ops[operand]
            except KeyError:
                raise ValueError("Unrecognised operand '{}' ".format(operand) +
                                 "in filter expression '{} {} {}'".format(
                                     field, operand, value))
            if field not in self.metadata:
                raise ValueError("{} field ".format(self.field_name) +
                                 "'{}' not in VCF ".format(field) +
                                 "header - can not be used for {} ".format(
                                     self.field_name) + "field filtering.")
            else:
                ftype = self.metadata[field].type
                coerc = None
                if ftype == 'Integer':
                    coerc = int
                elif ftype == 'Float':
                    coerc = float
                if coerc is not None:
                    try:
                        value = coerc(value)
                    except ValueError:
                        raise ValueError("Filter value for {} field ".format(
                                             self.field_name) +
                                         "'{}' could not be ".format(field) +
                                         "converted to {}, but".format(ftype) +
                                         "but field Type is {}".format(ftype) +
                                         "in VCF header.")
                elif ftype == 'Flag':
                    if value.lower() == 'true' or value == '1':
                        value = True
                    elif value.lower() == 'false' or value == '0':
                        value = False
                    else:
                        raise ValueError("INFO field '{}' is ".format(field) +
                                         "a Flag but value passed in filter " +
                                         "expression ('{}') ".format(value) +
                                         "can not be interpreted as a " +
                                         "boolean. Supported values are " +
                                         "'True', 'False', '1' or '0'.")
                    if op != operator.eq and op != operator.ne:
                        raise ValueError("INFO field '{}' is ".format(field) +
                                         "a Flag but '{}' ".format(operand) +
                                         "operand is neither '==' or '!='. " +
                                         "Flag based filters must be tests " +
                                         "of equality versus boolean values.")
                num = self.metadata[field].number
                f.append((field, op, value, num))
                self.fields.add(field)
        return f

    def filter(self, record, key=None):
        '''
            Read VcfRecord and return a list of booleans indicating
            whether each ALT allele should be filtered according to
            parameters from self.filters.

            record:    VcfRecord to assess according to self.filters.

            key:       If filtering on FORMAT fields, this key should
                       indicate which sample to assess.
        '''
        filter_alleles = [False] * len(record.alts)
        if key is not None:
            annots = getattr(record, self.record_field)[key]
        else:
            annots = getattr(record, self.record_field)
        for field, op, value, number in self.filters:
            if number == 0:  # flag
                flag_value = False
                if field in annots:
                    flag_value = True
                if not op(value, flag_value):
                    filter_alleles = [True] * len(record.alts)
                    break
                continue
            if field not in annots:
                continue
            if number == 'A' or number == 'R':
                if number == 'R':
                    alt_vals = annots[field][1:]
                else:
                    alt_vals = annots[field]
                for i in range(len(filter_alleles)):
                    if alt_vals[i] is None or not op(alt_vals[i], value):
                        filter_alleles[i] = True
            else:
                if number == 1:
                    if annots[field] is not None and not op(
                            annots[field], value):
                        filter_alleles = [True] * len(record.alts)
                else:
                    for x in annots[field]:
                        if x is not None and not op(x, value):
                            filter_alleles = [True] * len(record.alts)
                            break
            if all(filter_alleles):
                # all alleles filtered
                break
        return filter_alleles
