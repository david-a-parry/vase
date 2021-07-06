from collections import OrderedDict
import re

sv_fields = ['SVTYPE', 'CIPOS', 'CIEND', 'SVLEN', 'IMPRECISE',
             'LEFT_SVINSSEQ', 'RIGHT_SVINSSEQ']
_svalt_re = re.compile(r'<(\w+)(:\w+)*>')  # group 1 gives SV type
_bnd_re = re.compile(r'^(([ACTGN]*)[\[\]]\w+):\d+[\]\[]([ACGTN]*)$')
# group 1 gives VEP CSQ allele


class VaseRecord(object):
    """
        A class providing convenience methods for parsing features from
        VCF records.
    """

    __slots__ = ['record', 'caller', 'header', '__CSQ', '__ANN', '__is_sv',
                 '__DECOMPOSED_ALLELES', '_vep_allele']

    def __init__(self, record, vcfreader):
        """
            Create a new VaseRecord object.

            Args:
                record:    A libcbcf.VariantRecord object from pysam

                vcfreader: VcfReader object from which record was generated

        """
        self.record = record
        self.caller = vcfreader
        self.header = self.caller.header
        self.__CSQ = None
        self.__ANN = None
        self.__is_sv = None
        self.__DECOMPOSED_ALLELES = None
        self._vep_allele = {}

    def __str__(self):
        return str(self.record)

    @property
    def alleles(self):
        return self.record.alleles

    @alleles.setter
    def alleles(self, x):
        self.record.alleles = x

    @property
    def alts(self):
        return self.record.alts

    @alts.setter
    def alts(self, x):
        self.record.alts = x

    @alts.setter
    def alts(self, x):
        self.record.alts = x

    @property
    def alt(self):
        return ",".join(self.record.alts)

    @property
    def info(self):
        return self.record.info

    @info.setter
    def info(self, x):
        self.record.info = x

    @property
    def filter(self):
        return self.record.filter

    @filter.setter
    def filter(self, x):
        self.record.filter = x

    @property
    def format(self):
        return self.record.format

    @format.setter
    def format(self, x):
        self.record.format = x

    @property
    def chrom(self):
        return self.record.chrom

    @chrom.setter
    def chrom(self, x):
        self.record.chrom = x

    @property
    def pos(self):
        return self.record.pos

    @pos.setter
    def pos(self, x):
        self.record.pos = x

    @property
    def id(self):
        return self.record.id

    @id.setter
    def id(self, x):
        self.record.id = x

    @property
    def qual(self):
        return self.record.qual

    @property
    def ref(self):
        return self.record.ref

    @property
    def rlen(self):
        return self.record.rlen

    @property
    def samples(self):
        return self.record.samples

    @property
    def start(self):
        return self.record.start

    @property
    def stop(self):
        return self.record.stop

    @property
    def IS_SV(self):
        '''True if record represents a structural variant'''
        if self.__is_sv is None:
            self.__is_sv = 'SVTYPE' in self.record.info
        return self.__is_sv

    @IS_SV.setter
    def IS_SV(self, is_sv):
        self.__is_sv = is_sv

    @property
    def DECOMPOSED_ALLELES(self):
        '''
            list of AltAllele objects, one for each ALT allele in
            order, after reducing them to their minimal representations
            (i.e. by trimming redundant nucleotides).
        '''

        if self.__DECOMPOSED_ALLELES is None:
            self._minimize_alleles()
        return self.__DECOMPOSED_ALLELES

    @DECOMPOSED_ALLELES.setter
    def DECOMPOSED_ALLELES(self, alleles):
        self.__DECOMPOSED_ALLELES = alleles

    def _minimize_alleles(self):
        self.DECOMPOSED_ALLELES = []
        for alt in self.record.alleles[1:]:
            if self.IS_SV:
                self.DECOMPOSED_ALLELES.append(
                    AltAllele(chrom=self.record.chrom,
                              pos=self.record.pos,
                              ref=self.record.ref,
                              alt=alt,
                              is_sv=True,
                              record=self))
            else:
                ref = self.record.ref
                pos = self.record.pos
                while len(ref) > 1 and len(alt) > 1:
                    if ref[-1] == alt[-1]:  # remove identical suffixes
                        ref = ref[:-1]
                        alt = alt[:-1]
                    else:
                        break
                while len(ref) > 1 and len(alt) > 1:
                    if ref[0] == alt[0]:  # remove identical prefixes
                        ref = ref[1:]
                        alt = alt[1:]
                        pos += 1
                    else:
                        break
                self.DECOMPOSED_ALLELES.append(
                    AltAllele(chrom=self.record.chrom,
                              pos=pos,
                              ref=ref,
                              alt=alt))

    @property
    def CSQ(self):
        '''
            A list of dicts of CSQ annotations from VEP to values.
            Empty values are represented by empty Strings. Will raise
            a HeaderError if the associated VCF header does not contain
            CSQ information and a ParseError if the record being
            parsed does not contain a CSQ annotation in the INFO
            field.
        '''
        if self.__CSQ is None:
            lbl = self.header.csq_label
            try:
                csqs = self.record.info[lbl]
            except KeyError:
                raise ValueError("Could not find '{}' label in ".format(lbl) +
                                 "INFO field of record at {}:{}"
                                 .format(self.chrom, self.pos))
            self.__CSQ = []
            for c in csqs:
                d = OrderedDict([(k, v) for (k, v) in zip(
                    self.header.csq_fields, c.split('|'))])
                if len(self.record.alleles) == 2:  # only one ALT allele
                    d['alt_index'] = 1
                elif 'ALLELE_NUM' in d:
                    d['alt_index'] = int(d['ALLELE_NUM'])
                else:
                    d['alt_index'] = self._vep_to_alt(d)
                self.__CSQ.append(d)
        return self.__CSQ

    @CSQ.setter
    def CSQ(self, c):
        self.__CSQ = c

    @property
    def ANN(self):
        '''
            A list of dicts of ANN/EFF annotations from SnpEff to values.
            Empty values are represented by empty Strings. Will raise
            a HeaderError if the associated VCF header does not contain
            ANN/EFF information and a ParseError if the record being
            parsed does not contain a ANN/EFF annotation in the INFO
            field.
        '''
        if self.__ANN is None:
            lbl = self.header.ann_label
            try:
                anns = self.record.info[lbl]
            except KeyError:
                raise ValueError("Could not find '{}' label in ".format(lbl) +
                                 "INFO field of record at {}:{}"
                                 .format(self.chrom, self.pos))
            self.__ANN = [dict(zip(self.header.ann_fields, x.split('|'))) for x
                          in anns]
            for d in self.__ANN:
                d['alt_index'] = self.alleles.index(d['Allele'])
        return self.__ANN

    @ANN.setter
    def ANN(self, ann):
        self.__ANN = ann

    def _vep_to_alt(self, csq):
        # figure out how alleles will be handled by looking at the REF vs ALTs
        allele = csq['Allele']
        if allele in self._vep_allele:
            return self._vep_allele[allele]
        is_sv = False
        is_snv = False
        is_indel = False
        is_mnv = False
        ref = self.record.ref
        asterisk = False
        for i in range(1, len(self.record.alleles)):
            alt = self.record.alleles[i]
            if alt == '*':
                self._vep_allele[alt] = i
                asterisk = True
            else:
                matches_sv = _svalt_re.match(alt)
                matches_bnd = _bnd_re.match(alt)
                if matches_sv or matches_bnd:
                    is_sv = True
                    # sometimes VEP unhelpfully just uses '-'
                    if allele == '-':
                        sv_type = '-'
                    elif matches_sv:
                        sv_type = matches_sv.group(1)
                    else:
                        sv_type = matches_bnd.group(1)
                    if sv_type == 'DUP':
                        self._vep_allele['duplication'] = i
                    elif sv_type == 'INS':
                        self._vep_allele['insertion'] = i
                    elif sv_type == 'DEL':
                        self._vep_allele['deletion'] = i
                    else:
                        # should catch CNVs, INVs, BNDs
                        self._vep_allele[sv_type] = i
                else:
                    if len(alt) == 1 and len(ref) == 1:
                        if alt != ref:
                            is_snv = True
                    elif len(alt) == len(ref):
                        is_mnv = True
                    else:
                        is_indel = True
                    if is_indel:
                        # special case for longer non SV type 'deletion'
                        # 'insertion' or 'duplication' alleles which VEP
                        # sometimes annotates as deletion/insertion/duplication
                        # despite presence of REF/ALT sequences
                        if allele == 'deletion' and len(alt) < len(ref):
                            self._vep_allele[allele] = i
                            return i
                        elif allele == 'insertion' and len(alt) > len(ref):
                            self._vep_allele[allele] = i
                            return i
                        elif allele == 'duplication' and len(alt) > len(ref):
                            self._vep_allele[allele] = i
                            return i
                    self._vep_allele[alt] = i

        if is_sv:
            # no more editing required as long as
            # not at the same site as short variant
            if is_snv or is_mnv or is_indel:
                raise ValueError("Unable to parse structural variants at the "
                                 + "same site as a non-structural variant")
        else:
            if not is_snv and (is_indel or
                               (is_mnv and asterisk)):
                # VEP trims first base unless REF and ALT differ at first base
                first_base_differs = False
                ref_start = ref[:1]
                for alt in self.record.alleles[1:]:
                    if alt != '*':
                        alt_start = alt[:1]
                        if alt_start != ref_start:
                            first_base_differs = True
                            break
                if not first_base_differs:
                    # no trimming if first base differs for any ALT,
                    # otherwise first base is trimmed
                    trimmed = {}
                    pop = []
                    for alt in self._vep_allele:
                        if alt != '*':
                            i = self._vep_allele[alt]
                            pop.append(alt)
                            if len(alt) > 1:
                                alt = alt[1:]
                            else:
                                alt = '-'
                            trimmed[alt] = i
                    for p in pop:
                        self._vep_allele.pop(p, None)
                    self._vep_allele.update(trimmed)
        return self._vep_allele[allele]

    def in_cis_with(self, sample, allele, other, other_allele):
        '''
            Returns True if the two alleles are physically phased
            according to the GT or PID and PGT fields of both records.

            Args:
                sample: Sample ID to check phasing data for.

                allele: Allele number of this record.

                other:  Other record to compare with this record.

                other_allele:
                        Allele number for other record.

        '''
        # If GT field is phased in both records check this and only this
        if (self.record.samples[sample].phased and
                other.record.samples[sample].phased):
            i = self.record.samples[sample].allele_indices.index(allele)
            j = other.record.samples[sample].allele_indices.index(other_allele)
            return i == j
        # Otherwise check phase group annotations
        if 'PID' not in self.record.format or 'PID' not in other.record.format:
            return False
        if 'PGT' not in self.record.format or 'PGT' not in other.record.format:
            return False
        try:
            pid1 = self.record.samples[sample]['PID']
            pid2 = other.record.samples[sample]['PID']
            pgt1 = self.record.samples[sample]['PGT']
            pgt2 = other.record.samples[sample]['PGT']
        except KeyError:
            # when joining VCFs together only some samples may have PID/PGT
            return False
        if pid1 != pid2:
            return False
        if pgt1 == '.' or pgt2 == '.':
            return False
        try:
            phase1 = pgt1.split('|').index(str(allele))
            phase2 = pgt2.split('|').index(str(other_allele))
            return phase1 == phase2
        except ValueError:  # allele might not be in phase group
            return False

    def add_info_fields(self, info, append_existing=False):
        '''
            Requires a dict of INFO field names to a list of values.
            Adds or replaces existing INFO fields in the record with
            the items in given dict.

            Args:
                info: A dict of INFO field names to add with values
                      being list of values for the given field.

                append_existing:
                      Add values to existing INFO fields in a record.
                      If the field being added already exists and this
                      argument is True, the values provided will be
                      added to the existing values. If the Number
                      property is a fixed value, multiple values at the
                      same index will be separated by '|' characters.
                      Otherwise they will be separated by commas. Note,
                      that the latter method is only supported for
                      'String' types.
                      Default = False.

        '''
        for k, v in sorted(info.items()):
            if append_existing and k in self.info:
                self._append_to_existing_info(k, v)
            else:
                self.info[k] = v

    def _append_to_existing_info(self, field, values):
        if field not in self.header.info:
            raise KeyError("{} INFO field does not exist in ".format(field) +
                           "header - can not add to record")
        if self.header.info[field].type == 'Flag':
            self.info[field] = True
            return
        elif self.header.info[field].number == '.':
            self.info[field] = ",".join(str(x) for x in self.info[field] +
                                        tuple(values))
            return
        if self.header.info[field].type != 'String':
            raise ValueError("Appending to fixed length INFO fields is only " +
                             "supported for String types, not {} types".format(
                                 self.header.info[field].type))
        if self.header.info[field].number == '1':
            self.info[field] += "|" + values
            return
        if (len(self.info[field]) != len(new)):
            raise ValueError("New {} INFO field '{}'" .format(field, values) +
                             "has differing number of values to existing " +
                             "field '{}'" .format(self.info[field]))
        self.info[field] = str.join(",", (str.join("|", x) for x in zip(
            self.info[field], values)))

    def add_ids(self, ids, replace=False):
        '''
            Adds given IDs to the ID field of the VCF record. If the
            record already has an ID (i.e. is not '.') these IDs are
            added to the existing value(s) unless the replace
            argument is True.

            Args:
                ids:     A list of IDs to add.

                replace: If True, existing ID values are replaced,
                         otherwise the given IDs are added to.
                         Default = False.

        '''
        if replace or self.id is None:
            self.id = str.join(';', ids)
        else:
            uids = set(ids + self.id.split(';'))
            self.id = str.join(';', uids)


class AltAllele(object):
    '''
        Represents basic genomic features of a single alternative
        allele call. Features are 'CHROM', 'POS', 'REF' and 'ALT'.
    '''

    __slots__ = ['CHROM', 'POS', 'REF', 'ALT', 'is_sv', '__var_type',
                 'sv_info', 'breakpoint_precision']

    def __init__(self, chrom, pos, ref, alt, is_sv=False, record=None,
                 breakpoint_precision=0.1):
        '''
            Either created from a given VcfRecord and the index of the
            allele to be represented or from chrom, pos, ref and alt
            arguments.

            Args:
                chrom:     chromosome/contig.

                pos:       position of ref allele.

                ref:       reference allele.

                alt:       alternative allele.

                is_sv:     True if record represents a structural variant.
                           Default=False.

                record     VcfRecord object containing the allele of
                           interest. Required if allele is a structural
                           variant.

                breakpoint_precision:
                           Structural variants are considered equal if
                           breakpoints differ by this fraction or less
                           of the length of the smaller variant. For
                           example, when comparing a 1000 bp deletion to
                           another equally sized or larger deletion, if
                           breakpoint_precision is set to 0.1 then
                           breakpoints may differ by 100bp and the two
                           variants will be considered equal. If 'CIPOS'
                           and 'CIEND' fields are available, these will
                           be used to set the lower and upper bounds.
                           Set this to 0.0 to endpoints to be exactly
                           matching or within the 'CIPOS' and 'CIEND
                           intervals (if available). Default=0.1.


        '''
        self.CHROM = chrom
        self.POS = pos
        self.REF = ref
        self.ALT = alt
        self.__var_type = None
        self.is_sv = is_sv
        self.breakpoint_precision = breakpoint_precision
        self.sv_info = dict()
        if self.is_sv:
            if record is None:
                raise ValueError("record argument is required if AltAllele " +
                                 "is a structural variant (is_sv == True)")
            self.sv_info['END'] = record.stop
            for x in sv_fields:
                self.sv_info[x] = record.info[x] if x in record.info else None
            if isinstance(self.sv_info['SVLEN'], tuple): # Manta gives Number as '.'
                self.sv_info['SVLEN'] = self.sv_info['SVLEN'][0]

    @property
    def var_type(self):
        '''
            String indicating whether the ALT is a DELETION, INSERTION,
            SNV, MNV or SV.
        '''
        if self.__var_type is None:
            if self.is_sv:
                self.__var_type = 'SV'
            elif len(self.REF) < len(self.ALT):
                self.__var_type = 'INSERTION'
            elif len(self.REF) > len(self.ALT):
                self.__var_type = 'DELETION'
            else:
                if len(self.REF) == 1:
                    self.__var_type = 'SNV'
                else:
                    self.__var_type = 'MNV'
        return self.__var_type

    @var_type.setter
    def var_type(self, t):
        self.__var_type = t

    def __str__(self):
        if self.is_sv:
            return "{}:{}-{}_{}/{}".format(self.CHROM,
                                           self.POS,
                                           self.sv_info['END'],
                                           self.REF,
                                           self.ALT)
        return "{}:{}-{}_{}/{}".format(self.CHROM,
                                       self.POS,
                                       len(self.REF) + self.POS - 1,
                                       self.REF,
                                       self.ALT)

    def __eq__(self, other):
        if self.is_sv:
            return self._compare_svs(other)
        elif other.is_sv:
            return False
        return (self.CHROM == other.CHROM and self.POS == other.POS and
                self.REF == other.REF and self.ALT == other.ALT)

    def _compare_svs(self, other):
        '''
            Return True if both are SVs and have same type and
            breakpoints. Assumes self is a SV.

            Args:
                self:   This AltAllele

                other:  Another AltAllele

        '''
        if self.CHROM != other.CHROM:
            return False
        if self.sv_info['SVTYPE'] != other.sv_info['SVTYPE']:
            return False
        if self.sv_info['SVTYPE'] == 'BND':
            return self.compare_bnd(other)
        elif self.sv_info['SVTYPE'] == 'INS':
            return self.compare_svins(other)
        else:
            return self.compare_sv_pos_end(other)

    def compare_sv_pos_end(self, other):
        s_len = abs(self.sv_info['SVLEN'])
        o_len = abs(other.sv_info['SVLEN'])
        prec_margin = self.breakpoint_precision * min(s_len, o_len)
        # check start POS are close enough to each other
        if self.POS != other.POS:
            s_pos_int = (self.POS - prec_margin, self.POS + prec_margin)
            o_pos_int = (other.POS - prec_margin, other.POS + prec_margin)
            if self.sv_info['CIPOS'] is not None:
                s_pos_int = (self.POS - self.sv_info['CIPOS'][0] - prec_margin,
                             self.POS + self.sv_info['CIPOS'][1] + prec_margin)
            if other.sv_info['CIPOS'] is not None:
                o_pos_int = (other.POS - other.sv_info['CIPOS'][0] -
                             prec_margin,
                             other.POS + other.sv_info['CIPOS'][1] +
                             prec_margin)
            if s_pos_int[0] > o_pos_int[1] or s_pos_int[1] < o_pos_int[0]:
                return False
        # check end POS are close enough to each other
        if self.sv_info['END'] != other.sv_info['END']:
            s_end_int = (self.sv_info['END'] - prec_margin,
                         self.sv_info['END'] + prec_margin)
            o_end_int = (other.sv_info['END'] - prec_margin,
                         other.sv_info['END'] + prec_margin)
            if self.sv_info['CIEND'] is not None:
                s_end_int = (self.sv_info['END'] - self.sv_info['CIEND'][0] -
                             prec_margin,
                             self.sv_info['END'] + self.sv_info['CIEND'][1] +
                             prec_margin)
            if other.sv_info['CIEND'] is not None:
                o_end_int = (other.sv_info['END'] - other.sv_info['CIEND'][0],
                             other.sv_info['END'] + other.sv_info['CIEND'][1])
            if s_end_int[0] > o_end_int[1] or s_end_int[1] < o_end_int[0]:
                return False
        return True

    def compare_bnd(self, other):
        '''
            Crude check on two BNDs - only checks standard POS/
            REF/ALT fields for equality.
        '''
        return (self.POS == other.POS and self.REF == other.REF and
                self.ALT == other.ALT)

    def compare_svins(self, other):
        # for these purposes we don't care if exact insertion is the same,
        # merely whether it is of roughly equal length
        # if _svalt_re.match(self.ALT) and _svalt_re.match(other.ALT):
        if self.sv_info['LEFT_SVINSSEQ'] is not None:
            return self.compare_svinvseq(other)
        if other.sv_info['LEFT_SVINSSEQ'] is not None:
            return other.compare_svinvseq(self)
        if not self.compare_sv_pos_end(other):
            return False
        if (self.sv_info['SVLEN'] is not None and other.sv_info['SVLEN'] is not
                None):
            ln = sorted([self.sv_info['SVLEN'], other.sv_info['SVLEN']])
            if ln[1] <= ln[0] + (self.breakpoint_precision * ln[0]):
                return True
            else:
                return False
        return True

    def compare_svinvseq(self, other):
        '''
            For Manta large insertions that are not fully assembled,
            return True if LEFT_SVINVSEQ and RIGHT_SVINVSEQ are
            identical for self and other.
        '''
        if self.POS != other.POS:
            return False
        if 'LEFT_SVINSSEQ' not in other.sv_info:
            return False
        return (self.sv_info['LEFT_SVINSSEQ'] == other.sv_info['LEFT_SVINSSEQ']
                and self.sv_info['RIGHT_SVINSSEQ'] ==
                other.sv_info['RIGHT_SVINSSEQ'])
