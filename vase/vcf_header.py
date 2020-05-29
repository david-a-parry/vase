import re

_csq_format_re = re.compile(r'''.*Format:\s*((\S+\|)*\S+)''')
# for capturing CSQ format in Description field of metaheader
_common_csq_fields = ['CSQ', 'ANN', 'BCSQ', 'CQ']


class VcfHeader(object):
    ''' Header class storing metadata and sample information for a vcf '''

    def __init__(self, vcfreader):
        self.vcfreader = vcfreader
        self.header = self.vcfreader.variant_file.header
        self.formats = self.header.formats
        self.info = self.header.info
        self.filters = self.header.filters
        self.csq_fields = None

    @property
    def csq_label(self):
        '''
            String labelling the INFO field label of VEP consequence
            annotations. Will raise a KeyError if access is attempted
            but no VEP CSQ or ANN field is present in the header.
        '''
        if self.__csq_label is None:
            self.csq_fields
        return self.__csq_label

    @csq_label.setter
    def csq_label(self, c):
        self.__csq_label = c

    @property
    def csq_fields(self):
        '''
            A list of CSQ field names in the order they are represented
            in CSQ INFO field entries. Set to None on initialization.
            Will raise a KeyError if access is attempted but no VEP
            CSQ, ANN, BCSQ or CQ field is present in the header.
        '''

        if self.__csq_fields is None:
            if self.__csq_label is None:
                csq = None
                for x in _common_csq_fields:
                    if x in self.info:
                        csq = x
                        break
                if csq is None:
                    raise KeyError("No common CSQ fields found in INFO " +
                                   "header - unable to retrieve consequence " +
                                   "fields.")
                self.csq_label = csq
            else:
                csq = self.__csq_label
            csq_header = self.info[csq]
            match = self._csq_format_re.match(csq_header.description)
            if match:
                self.__csq_fields = match.group(1).split('|')
            else:
                raise KeyError("Could not parse {} Format in ".format(csq)
                               + "header. Unable to retrieve consequence "
                               + "annotations.")
        return self.__csq_fields

    @csq_fields.setter
    def csq_fields(self, csq):
        self.__csq_fields = csq
