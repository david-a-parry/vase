import sys
from .parse_vcf.parse_vcf import * 

default_csq = set(
               ["TFBS_ablation",
                "TFBS_amplification",
                "inframe_deletion",
                "inframe_insertion",
                "frameshift_variant",
                "initiator_codon_variant",
                "missense_variant",
                "protein_altering_variant",
                "regulatory_region_ablation",
                "regulatory_region_amplification",
                "splice_acceptor_variant",
                "splice_donor_variant",
                "stop_gained",
                "stop_lost",
                "transcript_ablation",
                "transcript_amplification",]
)

valid_csq = set(
               [
                "3_prime_UTR_variant",
                "5_prime_UTR_variant",
                "NMD_transcript_variant",
                "TFBS_ablation",
                "TFBS_amplification",
                "TF_binding_site_variant",
                "coding_sequence_variant",
                "downstream_gene_variant",
                "feature_elongation",
                "feature_truncation",
                "frameshift_variant",
                "incomplete_terminal_codon_variant",
                "inframe_deletion",
                "inframe_insertion",
                "initiator_codon_variant",
                "intergenic_variant",
                "intron_variant",
                "mature_miRNA_variant",
                "missense_variant",
                "nc_transcript_variant",
                "non_coding_exon_variant",
                "protein_altering_variant",
                "regulatory_region_ablation",
                "regulatory_region_amplification",
                "regulatory_region_variant",
                "splice_acceptor_variant",
                "splice_donor_variant",
                "splice_region_variant",
                "stop_gained",
                "stop_lost",
                "stop_retained_variant",
                "synonymous_variant",
                "transcript_ablation",
                "transcript_amplification",
                "upstream_gene_variant",]
)

default_biotypes = set(
               [
                "3prime_overlapping_ncrna",
                "antisense",
                "CTCF_binding_site",
                "enhancer",
                "IG_C_gene",
                "IG_D_gene",
                "IG_J_gene",
                "IG_V_gene",
                "lincRNA",
                "miRNA",
                "misc_RNA",
                "Mt_rRNA",
                "Mt_tRNA",
                "open_chromatin_region",
                "polymorphic_pseudogene",
                "processed_transcript",
                "promoter",
                "promoter_flanking_region",
                "protein_coding",
                "rRNA",
                "sense_intronic",
                "sense_overlapping",
                "snoRNA",
                "snRNA",
                "TF_binding_site",
                "translated_processed_pseudogene",
                "TR_C_gene",
                "TR_D_gene",
                "TR_J_gene",
                "TR_V_gene",]
)

valid_biotypes = set(
               [
                "3prime_overlapping_ncrna",
                "antisense",
                "CTCF_binding_site",
                "enhancer",
                "IG_C_gene",
                "IG_C_pseudogene",
                "IG_D_gene",
                "IG_J_gene",
                "IG_J_pseudogene",
                "IG_V_gene",
                "IG_V_pseudogene",
                "lincRNA",
                "miRNA",
                "misc_RNA",
                "Mt_rRNA",
                "Mt_tRNA",
                "nonsense_mediated_decay",
                "non_stop_decay",
                "open_chromatin_region",
                "polymorphic_pseudogene",
                "processed_pseudogene",
                "processed_transcript",
                "promoter",
                "promoter_flanking_region",
                "protein_coding",
                "pseudogene",
                "retained_intron",
                "rRNA",
                "sense_intronic",
                "sense_overlapping",
                "snoRNA",
                "snRNA",
                "TF_binding_site",
                "transcribed_processed_pseudogene",
                "transcribed_unprocessed_pseudogene",
                "translated_processed_pseudogene",
                "TR_C_gene",
                "TR_D_gene",
                "TR_J_gene",
                "TR_J_pseudogene",
                "TR_V_gene",
                "TR_V_pseudogene",
                "unitary_pseudogene",
                "unprocessed_pseudogene",]
)


class VepFilter(object):
    '''An object that filters VCF records based on annotated VEP data.'''

    def __init__(self, csq=[], canonical=False, biotypes=[]):
        self.csq = set()
        self.biotypes = set()
        if len(csq) == 0:
            csq = ['default']
        for c in csq:
            if c.lower() == 'default':
                self.csq.update(default_csq)
            else:
                if c.lower() in valid_csq:
                    self.csq.add(c.lower())
                else:
                    raise Exception("ERROR: Unrecognised VEP consequence " +  
                                    "class '{}'".format(c))
        if len(biotypes) == 0:
            biotypes = ['default']
        for b in biotypes:
            if b.lower() == 'all':
                self.biotypes.update(valid_biotypes)
            elif b.lower() == 'default':
                self.biotypes.update(default_biotypes)
            else:
                if b.lower() in valid_biotypes:
                    self.biotypes.add(b.lower())
                else:
                    raise Exception("ERROR: Unrecognised VEP biotype " +  
                                    "'{}'".format(b))

        self.canonical = canonical

    def filter(self, record):
        filter_alleles = [True] * (len(record.ALLELES) -1)
        for c in record.CSQ:
            if self.canonical:
                try: 
                    if c['CANONICAL'] != 'YES':
                        continue
                except KeyError:
                    pass
            if c['BIOTYPE'] not in self.biotypes:
                continue
            consequence = c['Consequence'].split('&')
            for s_csq in consequence:
                if s_csq in self.csq:
                    filter_alleles[c['alt_index'] -1] = False
        return filter_alleles
