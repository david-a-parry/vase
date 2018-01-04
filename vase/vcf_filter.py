import sys
from parse_vcf import * 


class VcfFilter(object):
    ''' 
        An object that filters VCF records based on variant data in a 
        another VCF file.
    '''
    
    def __init__(self, vcf, prefix, freq=None, min_freq=None, 
                 freq_fields=("AF",), ac_fields=("AC",), an_fields=("AN",) ):
        ''' 
            Initialize object with a VCF file and optional filtering 
            arguments.
            
            Args:
                vcf:      VCF containing variants to use to filter or 
                          annotate records.
                    
                prefix:   Prefix to prepend to added INFO field 
                          annotations. Required.

                freq:     Filter alleles if allele frequency is greater
                          than this value. Optional.

                min_freq: Filter alleles if allele frequency is less 
                          than this value. Optional.

        '''

        self.vcf = VcfReader(vcf)
        self.prefix = prefix
        self.freq = freq
        self.min_freq = min_freq
        self.freq_info = freq_fields
        self.ac_info = ac_fields
        self.an_info = an_fields
        if self.freq is not None and self.min_freq is not None:
            if self.freq <= self.min_freq:
                raise Exception("freq argument must be greater than " +
                                "min_freq argument")
        self.freq_fields = dict()
        self.ac_fields = dict()
        self.an_fields = dict()
        self.annot_fields = dict()
        self.get_annot_fields()
        self.added_info = {}
        self.create_header_fields()

    def get_overlapping_records(self, record):
        ''' 
            For a given record, returns a list of overlapping records 
            in the class's VCF.
        '''

        start = record.POS
        end = record.SPAN
        self.vcf.set_region(record.CHROM, start - 1, end)
        return list(s for s in self.vcf.parser)

    def annotate_and_filter_record(self, record):
        ''' 
            For a given record add AF annotations per allele and returns 
            two lists of boolean values. The first list has a value for  
            each ALT allele, in order, indicating whether to filter the 
            allele or not, according to the settings for self.freq and 
            self.min_freq. The second list is the same as the first but 
            is meant for indicating whether alleles should be kept, 
            overriding values in the first list, which is unused in the 
            'VcfFilter' class but is provided for use by inheriting 
            classes.

        '''

        filter_alleles = []
        keep_alleles = []
        matched_alleles = []
        annotations = []
        hits = self.get_overlapping_records(record)
        all_annots = set() #all fields added - may not be present for every ALT
        for i in range(len(record.DECOMPOSED_ALLELES)):
            filt,keep,matched,annot = self._compare_var_values(
                                            record.DECOMPOSED_ALLELES[i], hits)
            filter_alleles.append(filt)
            keep_alleles.append(keep)
            matched_alleles.append(matched)
            annotations.append(annot)
            all_annots.update(annot.keys())
        info_to_add = {}
        for f in all_annots:
            f_name = self.prefix + "_" + f
            info_to_add[f_name] = []
            for i in range(len(record.DECOMPOSED_ALLELES)):
                if f in annotations[i]:
                    a_val = annotations[i][f]
                else:
                    a_val = '.'
                info_to_add[f_name].append(a_val)
        if info_to_add:
            record.add_info_fields(info_to_add)
        return filter_alleles, keep_alleles, matched_alleles

    def _compare_var_values(self, alt_allele, var_list):
        do_filter = False #only flag indicating should be filtered
        do_keep = False #flag to indicate that should be kept, for overriding 
                        #do_filter in downstream applications 
        annot = {}
        matched = False
        for var in var_list:
            for i in range(len(var.DECOMPOSED_ALLELES)):
                if alt_allele == var.DECOMPOSED_ALLELES[i]:
                    matched = True
                    for f,d in self.freq_fields.items():
                        val = self._get_value(f, d, var, i)
                        if val is not None:
                            annot[f] = val
                            if self.freq is not None:
                                try:
                                    if float(val) >= self.freq:
                                        do_filter = True
                                except ValueError: 
                                    pass 
                            if self.min_freq is not None:
                                try:
                                    if float(val) < self.min_freq:
                                        do_filter = True
                                except ValueError: 
                                    pass
                    for an_k,an_v in self.an_fields.items():
                        # will only have self.an_fields if there were no 
                        # self.freq_fields 
                        f = an_k.replace("AN", "AF", 1)
                        ac_k = an_k.replace("AN", "AC", 1)
                        ac_v = self.ac_fields[ac_k]
                        an = self._get_value(an_k, an_v, var, i)
                        ac = self._get_value(ac_k, ac_v, var, i)
                        if an is None or ac is None:
                            continue
                        try:
                            an = float(an)
                            ac = float(ac)
                            af = ac/an
                            annot[f] = "{:g}".format(af)
                            if self.freq is not None:
                                if af >= self.freq:
                                    do_filter = True
                            if self.min_freq is not None:
                                if af < self.min_freq:
                                    do_filter = True
                        except (ValueError, ZeroDivisionError): 
                            pass
                    for f,d in self.annot_fields.items():
                        val = self._get_value(f, d, var, i)
                        if val is not None:
                            annot[f] = val
                if matched: break
            if matched: break #bail out on first matching variant
        return (do_filter, do_keep, matched, annot)

    def _get_value(self, field_name, field_properties, variant, index):
        if field_name not in variant.INFO_FIELDS:
            return None
        a_offset = None
        val = None
        if field_properties['Number'] == 'A':
            a_offset = 0
        elif field_properties['Number'] == 'R':
            a_offset = 1
        if a_offset is not None:
            val = variant.INFO_FIELDS[field_name].split(',')[index + a_offset]
        elif len(variant.DECOMPOSED_ALLELES) == 1:
            # if not a value per allele, only process if 
            # var only has one ALT allele because we don't 
            # know which is the relevant allele
            val = variant.INFO_FIELDS[field_name]
        return val

    def get_annot_fields(self):
        '''
            Creates dicts of INFO field names to dicts of 'Type', 
            'Number'and 'Description' as found in the VCF metadata for 
            INFO field names for frequency (freq_fields = ['AF']),
            allele count fields (ac_fields = ['AC']) and allele 
            number fields (an_fields = ['AN'] or allele annotations 
            fields (annot_fields = ['AN', 'AC']).
            
            Override this method to use custom allele frequency and 
            annotation fields.
        '''

        annot_info = self.ac_info + self.an_info
        for f in self.freq_info:
            if f in self.vcf.metadata['INFO']:
                self.freq_fields[f] = self.vcf.metadata['INFO'][f][-1]
        for f in annot_info:
            if f in self.vcf.metadata['INFO']:
                self.annot_fields[f] = self.vcf.metadata['INFO'][f][-1]

        if not self.freq_fields and (self.freq is not None or 
                                     self.min_freq is not None):
            self._get_an_and_ac(self.an_info, self.ac_info)
            if not self.an_fields:
                raise Exception("ERROR: no frequency fields identified in " + 
                                "VCF header for file '{}'."
                                .format(self.vcf.filename) + " Unable to use" +
                                "freq/min_freq arguments for variant " + 
                                "filtering.")

    def create_header_fields(self):
        '''
            Create dict entries for all INFO fields added by this 
            instance, suitable for adding to a VcfHeader object.
        '''
        for f,v in self.annot_fields.items():
            self._make_metadata(f, v)
        if self.freq_fields:
            for f,v in self.freq_fields.items():
                self._make_metadata(f, v)
        elif self.an_fields:
            self._make_an_freq_metadata()

    def _make_an_freq_metadata(self):
        for an in self.an_fields:
            name = an.replace("AN", "AF", 1)
            ac = an.replace("AN", "AC", 1)
            desc = ('"{} INFO field parsed by {} object. '.format(
                    name, type(self).__name__) + "Calculated from {} and {} "
                    .format(ac, an) + "from file {}.".format(self.vcf.filename)
                   ) 
            self.added_info[self.prefix + "_" + name] = {'Number' : 'A', 
                                                         'Type' : 'Float',
                                                         'Description' : desc }

    def _make_metadata(self, name, properties):
       
        desc = ('"{} INFO field parsed by {} object from file {}. '.format(
                name, type(self).__name__, self.vcf.filename) + 
               'Original description was as follows: {}"' .format(
                properties['Description'].replace('"', '')))
        if properties['Type'] == 'Flag':
            f_type = 'Integer'
        else:
            f_type = properties['Type']
        self.added_info[self.prefix + "_" + name] = {'Number' : 'A', 
                                                     'Type' : f_type,
                                                     'Description' : desc }
         
    def _get_an_and_ac(self, an_info, ac_info):
        for i in range(len(ac_info)):
            an = an_info[i]
            ac = ac_info[i]
            if (an in self.vcf.metadata['INFO'] and 
                ac in self.vcf.metadata['INFO']):
                #ANs and ACs should be in same order
                self.an_fields[an] = self.vcf.metadata['INFO'][an][-1]
                self.ac_fields[ac] = self.vcf.metadata['INFO'][ac][-1]



