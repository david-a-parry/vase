import sys
sys.path.insert(0, '')
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


class VepFilter(object):
    '''An object that filters VCF records based on annotated VEP data.'''

    def __init__(self, csq, canonical=False, keep_nmd_transcripts=False):
        self.csq = set()
        for c in csq:
            if c.lower() == 'default':
                self.csq.update(default_csq)
            else:
                if c.lower in valid_csq:
                    self.csq.add(c.lower())
                else:
                    raise Exception("ERROR: Unrecognised VEP consequence " +  
                                    "class '{}'".format(c))
        self.canonical = canonical
        self.keep_nmd_transcripts = keep_nmd_transcripts

    def filter(self, record):
        filter_alleles = [True] * (len(record.ALLELES) -1)
        for c in record.CSQ:
            if self.canonical:
                try: 
                    if c['CANONICAL'] != 'YES':
                        continue
                except KeyError:
                        pass
            consequence = c['Consequence'].split('&')
            if not self.keep_nmd_transcripts:
                if 'NMD_transcript_variant' in consequence:
                    continue
            for s_csq in consequence:
                if s_csq in self.csq:
                    filter_alleles[c['alt_index'] -1] = False
        return filter_alleles
            
                
