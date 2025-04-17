import pysam
import pandas as pd
import bioframe as bf

def get_alignment_intervals(read):
    """
    Generate alignment intervals from a pysam read object based on its position and CIGAR string.
    
    Args:
        read (pysam.AlignedSegment): A pysam read object
        
    Returns:
        pd.DataFrame: DataFrame with columns ['chrom', 'start', 'end'] representing alignment intervals
    """
    if not read.cigartuples:
        return []
    
    intervals = []
    current_pos = read.reference_start
    
    for operation, length in read.cigartuples:
        # M: match/mismatch (consumes reference)
        # D: deletion (consumes reference)
        # N: skipped region (consumes reference)
        # S: soft clipping (does not consume reference)
        # I: insertion (does not consume reference)
        # H: hard clipping (does not consume reference)
        # P: padding (does not consume reference)
        # =: sequence match (consumes reference)
        # X: sequence mismatch (consumes reference)
        
        if operation in (0, 7, 8):  # M/=/X: alignment match
            intervals.append((current_pos, current_pos + length))
            current_pos += length
        elif operation in (1, 4, 5):  # I/S/H: insertion/soft clip/hard clip
            continue  # doesn't consume reference
        elif operation in (2, 3):  # D/N: deletion/skipped region
            current_pos += length
        else:
            continue
    
    # Merge adjacent intervals
    if not intervals:
        return []
    
    merged = [intervals[0]]
    for current in intervals[1:]:
        last = merged[-1]
        if current[0] == last[1]:  # adjacent intervals
            merged[-1] = (last[0], current[1])
        else:
            merged.append(current)
    
    # Convert to DataFrame
    if not merged:
        return pd.DataFrame(columns=['chrom', 'start', 'end'])
    
    chrom = read.reference_name
    df = pd.DataFrame(merged, columns=['start', 'end'])
    df.insert(0, 'chrom', chrom)
    return df


def overlap_gene(read: pysam.AlignedSegment, read_pos: pd.DataFrame, gene: pd.DataFrame):
    gene['chrom'] = gene['chrom'].astype(str)
    gene = gene[gene['chrom']==read.reference_name]
    
    gene_ol = bf.overlap(read_pos, gene, return_overlap=True)
    gene_ol = gene_ol[gene_ol['gene_id_'] != None]
    
    gene_id = None
    gene_name = None
    gene_id_names = gene_ol[['gene_id_', 'gene_name_']]
    gene_id_names.index = gene_id_names['gene_id_']
    
    if gene_ol.shape[0] ==1:
        gene_id = gene_ol.iloc[0, 6]
        gene_name = gene_ol.iloc[0, 7]
    elif gene_ol.shape[0] >1:
        gene_ol['overlap_len'] = gene_ol['overlap_end'] - gene_ol['overlap_start']
        tmp = gene_ol.sort_values(['overlap_len'], ascending=False).reset_index().loc[0,:]  
        gene_id = tmp['gene_id_']
        gene_name = tmp['gene_name_']
    else:
        pass
    
    return gene_id, gene_name, gene_id_names

def read_overlap(read: pysam.AlignedSegment, gene: pd.DataFrame, intron: pd.DataFrame):
    r_df = get_alignment_intervals(read)
    
    gene_id, gene_name, gene_id_names = overlap_gene(read, r_df, gene)
    splice = 'splice'
    
    if gene_id is None:
        return None, None, None
    intron = intron[intron['gene_id'].isin(gene_id_names['gene_id_'])]
    intron_ol = bf.overlap(r_df, intron, return_overlap=True)
    intron_ol = intron_ol[intron_ol['gene_id_'] != None]
    intron_ol['overlap_len'] = intron_ol['overlap_end'] - intron_ol['overlap_start'] 
    
    if intron_ol.shape[0] == 1:
        if intron_ol['overlap_len'].head(1).values[0] > 10:
            splice = 'unsplice'
            gene_id = intron_ol.iloc[0, 7]
            gene_name = gene_id_names.loc[gene_id, "gene_name_"]
        
    elif intron_ol.shape[0] > 1:
        splice = 'unsplice'
        intron_ol['overlap_len'] = intron_ol['overlap_end'] - intron_ol['overlap_start']
        tmp = intron_ol.groupby(['gene_id_']).sum().sort_values(['overlap_len'], ascending=False).reset_index().loc[0,:] 
        if tmp['overlap_len'] > 10:
            splice = 'unsplice'
            gene_id = tmp['gene_id_']
            gene_name = gene_id_names.loc[gene_id, "gene_name_"]
    else:
        pass
    
    return gene_id, gene_name, splice
