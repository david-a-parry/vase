from .csq_filter import CsqFilter
import logging

class SnpEffFilter(CsqFilter):
    '''An object that filters VCF records based on annotated SnpEff data.'''

    def __init__(self, vcf, csq=[], impact=[], canonical=False, biotypes=[],
                 retain_labels=[], filter_flagged_features=False,
                 gene_filter=None, blacklist=None, g2p=None,
                 check_g2p_consequence=False, logging_level=logging.WARNING):
        '''
            Args:
                vcf:    input VcfReader object

                csq:    list of consequence types to keep. If 'default'
                        appears anywhere in this list then the default
                        consequence set (as indicated in
                        data/vep_classes.tsv) will be used. Similarly if
                        'all' appears anywhere in this list no filtering
                        on consequence type will occur.

                impact: list of variant impacts to retain.

                biotypes:
                        Filter consequences for features not of the
                        given biotypes. If not provided the default set
                        of biotypes (as indicated in data/biotypes.tsv)
                        will be used for biotype filtering.

                retain_labels:
                        Do not filter on consequence type if the
                        following values are present for a label. Labels
                        and values must be separated by '=' sign. For
                        example, to retain any consequence which has
                        a SnpEff annotation named 'FOO' with  value 'BAR'
                        use 'FOO=BAR'.

                filter_flagged_features:
                        Filter consequences on features which have warnings
                        from SnpEff's 'ERRORS / WARNINGS / INFO' field.

                gene_filter:
                        VarByRegion object from vase.var_by_region. If
                        provided, consequences will be filtered if they
                        do not alter the features specified in the
                        VarByRegion object for the current region.

                blacklist:
                        File containing a list of Feature IDs to ignore.

                g2p:
                        G2P object from vase.g2p for filtering on
                        presence and/or requirements from a G2P file.

                check_g2p_consequence:
                        If a G2P object is provided above, require that
                        that the observed consequence matches the
                        'mutation consequence' in the G2P file.

                logging_level:
                        Logging level to use. Default=logging.WARNING.

        '''
        self.logger = self._get_logger(logging_level)
        super().__init__(vcf=vcf, csq_attribute='ANN', csq=csq, impact=impact,
                         biotypes=biotypes, retain_labels=retain_labels,
                         filter_flagged_features=filter_flagged_features,
                         gene_filter=gene_filter, blacklist=blacklist, g2p=g2p,
                         check_g2p_consequence=check_g2p_consequence)

    def filter_csq(self, csq):
        '''
        Returns two boolean values. The first indicates whether the consequence
        annotation should be filtered. The second indicates whether the ALT
        allele should be filtered irrespective of the given or any other
        consequence annotation.
        '''
        if self.filter_flagged:
            try:
                eff_warnings = csq['ERRORS / WARNINGS / INFO']
                if eff_warnings is not None:
                    if eff_warnings.startswith(("ERROR", "WARNING")):
                        return True, False
            except KeyError:
                pass
        if (self.biotypes is not None and csq['Transcript_BioType'].lower()
                not in self.biotypes):
            return True, False
        if self.gene_filter:
            if not self.gene_filter.target_in_csq(csq):
                return True, False
        if self.g2p:
            if not any(x in self.g2p.g2p for x in csq['Gene_Name'].split('&')):
                return True, False
        if self.blacklist and csq['Feature_ID'] in self.blacklist:
            return True, False
        if (self.csq is None and self.impact is None and
                not self.check_g2p_consequence):
            # if only using biotypes/presence in G2P for filtering
            return False, False
        if self._retain_label_matched(csq):
            return False, False
        if self.check_g2p_consequence and self.g2p:
            filt_csq = self.g2p.consequences_from_gene(csq['Gene_Name'])
        else:
            filt_csq = self.csq
        for s_csq in (x.lower() for x in csq['Annotation'].split('&')):
            if filt_csq is not None and s_csq in filt_csq:
                return False, False
            if self.impact is not None and \
               csq['Annotation_Impact'] in self.impact:
                return False, False
        return True, False

    def get_required_header_fields(self):
        '''
        Check which CSQ/ANN annotation fields are required given arguments
        passed to __init__
        '''
        required = ['Allele', 'Annotation', 'Transcript_BioType']
        if self.impact:
            required.append('Annotation_Impact')
        if self.gene_filter or self.g2p:
            required.append('Gene_Name')
        if self.filter_flagged:
            required.append('ERRORS / WARNINGS / INFO')
        return required
