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


class InfoFilter(object):
    '''
        A class for filtering on given INFO fields in a VCF
    '''

    def __init__(self, vcf, filters, ):
        '''
            Args:
                vcf:    VcfReader object from parse_vcf.py

                filters:
                        iterable of tuples of field names, operands and
                        values for filtering.

        '''
        self.vcf = vcf
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
                raise ValueError("Unrecognised operand '{}'".format(operand))
            if field not in self.vcf.header.info:
                raise ValueError("INFO field '{}' not in VCF ".format(field) +
                                 "header - can not be used for INFO field " +
                                 "filtering.")
            else:
                ftype = self.vcf.header.info[field].type
                coerc = None
                if ftype == 'Integer':
                    coerc = int
                elif ftype == 'Float':
                    coerc = float
                if coerc is not None:
                    try:
                        value = coerc(value)
                    except ValueError:
                        raise ValueError("Filter value for INFO field " +
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
                num = self.vcf.header.info[field].number
                f.append((field, op, value, num))
                self.fields.add(field)
        return f

    def filter(self, record):
        '''
            Read and VcfRecord and return a list of booleans indicating
            whether each ALT allele should be filtered according to
            parameters from self.filters.
        '''
        filter_alleles = [False] * len(record.alts)
        for field, op, value, number in self.filters:
            if number == 0:  # flag
                flag_value = False
                if field in record.info:
                    flag_value = True
                if not op(value, flag_value):
                    filter_alleles = [True] * len(record.alts)
                    break
                continue
            if field not in record.info:
                continue
            if number == 'A' or number == 'R':
                if number == 'R':
                    alt_vals = record.info[field][1:]
                else:
                    alt_vals = record.info[field]
                for i in range(len(filter_alleles)):
                    if alt_vals[i] is None or not op(alt_vals[i], value):
                        filter_alleles[i] = True
            else:
                if number == 1:
                    if record.info[field] is not None and not op(
                            record.info[field], value):
                        filter_alleles = [True] * len(record.alts)
                else:
                    for x in record.info[field]:
                        if x is not None and not op(x, value):
                            filter_alleles = [True] * len(record.alts)
                            break
            if all(filter_alleles):
                # all alleles filtered
                break
        return filter_alleles
