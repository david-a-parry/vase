import re

_csq_format_re = re.compile(r'''.*Format:\s*((\S+\|)*\S+)''')
# for capturing CSQ format in Description field of metaheader
_common_csq_fields = ['CSQ', 'ANN', 'BCSQ', 'CQ']
_required_keys = {'info': ['number', 'type', 'description'],
                  'format': ['number', 'type', 'description'],
                  'filter': ['description'],
                  'alt': ['description']}
_field2pysam = {'info': 'info',
                'format': 'formats',
                'filter': 'filters',
                'alt': 'alts'}


class VcfHeader(object):
    ''' Header class storing metadata and sample information for a vcf '''

    __slots__ = ['vcfreader', 'header', '__csq_label', '__csq_fields']

    def __init__(self, vcfreader):
        self.vcfreader = vcfreader
        self.header = self.vcfreader.variant_file.header
        self.__csq_fields = None
        self.__csq_label = None

    @property
    def formats(self):
        return self.header.formats

    @property
    def info(self):
        return self.header.info

    @property
    def filters(self):
        return self.header.filters

    @property
    def samples(self):
        return self.header.samples

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
            match = _csq_format_re.match(csq_header.description)
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

    def add_header_field(self, name, string=None, field_type=None,
                         dictionary=None):
        '''
            Add a header field with given name and optional field type,
            and dictionary of properties.

            Args:
                name:   name of field to add

                string: string to add to field. Ignored if 'dictionary'
                        is provided.

                field_type:
                        type of field - e.g. if INFO/FILTER/FORMAT
                        field. Required if providing a dictionary.

                dictionary:
                        a dict of keys to values for the given field.
                        If 'field_type' is specified, this arg must be
                        provided and must contain all the essential keys
                        for that field type. For example, an 'INFO'
                        field must have 'Number', 'Type', and
                        'Description' keys.

        '''
        add_order = ['number', 'type', 'description']
        h_vals = []
        if dictionary is None and string is None:
            raise ValueError("Either dict or string argument is required")
        if field_type is not None and field_type in _required_keys:
            if dictionary is None:
                raise ValueError("Header type {} requires a dict.".format(
                    field_type))
        if dictionary:
            if not field_type:
                raise ValueError("field_type is required for use with " +
                                 "dictionary")
            dictionary.update([(k.lower(), v) for k, v in dictionary.items()])
            field_type = field_type.lower()
            field_header = getattr(self.header, _field2pysam[field_type])
            if name in field_header:
                field_header.remove_header(name)
                self.header = self.header.copy()
            if field_type in _required_keys:
                for k in add_order:
                    if k in _required_keys[field_type]:
                        try:
                            h_vals.append(dictionary[k])
                        except KeyError:
                            raise ValueError("Header type '" + field_type +
                                             "' requires '" + k + "' field")
                    else:
                        h_vals.append(None)
                h_vals.insert(0, name)
                getattr(self.header, _field2pysam[field_type]).add(*h_vals)
            else:
                raise ValueError("Field type {} not recognised".format(
                    field_type))
        else:
            self.header.add_meta(key=name, value=string)
