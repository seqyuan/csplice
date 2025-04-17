"""csplice core modules for genomic data processing."""
from .alignment_utils import get_alignment_intervals
from .gtftobed import (parse_attribute, gene_region, 
                      transcript_region, exon_region)

__all__ = [
    'get_alignment_intervals',
    'parse_attribute',
    'gene_region',
    'transcript_region', 
    'exon_region'
]
