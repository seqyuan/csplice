import click
import os
import pandas as pd
import bioframe as bf

from csplice.core.gtftobed import (parse_attribute, gene_region, transcript_region, exon_region)
from csplice.core.readoverlap import get_alignment_intervals, read_overlap

@click.group()
def main() -> None:
    """Command line for ."""
    pass

@main.command(name="gtf2bed")
@click.option('--gtf', '-g', required=True,
              help="gtf file path")
@click.option('--outdir', '-o', required=True,
              help="output directory")

@click.option('--gene_id', '-i', required=False, default="gene_id",
              help="gene id key word in column 9")
@click.option('--gene_name', '-n', required=False, default="gene_name",
              help="gene name key word in column 9")
@click.option('--transcript_id', '-t', required=False, default="transcript_id",
              help="transcript id key word in column 9")

def gtf_to_bed(gtf: str, outdir: str, gene_id: str, gene_name: str, transcript_id: str) -> None:
    """Process gtf to gene.bed transcript.bed exon.bed intron.bed."""
    os.makedirs(f"{outdir}", exist_ok=True)
    
    gene = gene_region(gtf, gid=gene_id, gname=gene_name)
    gene.to_csv(f"{outdir}/gene.bed", sep="\t", header=False, index=False)
    
    transcript = transcript_region(gtf, gid=gene_id, tid=transcript_id)
    transcript.to_csv(f"{outdir}/transcript.bed", sep="\t", header=False, index=False)
    
    exon = exon_region(gtf, gid=gene_id, tid=transcript_id)
    exon.to_csv(f"{outdir}/exon.bed", sep="\t", header=False, index=False)
    
    intron = bf.subtract(transcript, exon)
    intron.to_csv(f"{outdir}/intron.bed", sep="\t", header=False, index=False)


@main.command(name="bam2gene")
@click.option('--bam', '-b', required=True,
              help="alignment bam file")
@click.option('--outdir', '-o', required=True,
              help="output directory")
@click.option('--genebed', '-g', required=True,
              help="gene bed file path")
@click.option('--introned', '-i', required=True,
              help="intron bed file path")
@click.option('--cb', '-c', required=False, default="CB",
              help="Cell barcode tag (10X mobidrop should be CB")
@click.option('--ub', '-u', required=False, default="UB",
              help="UMI tag ((10X mobidrop should be UB")
@click.option('--threads', '-t', required=False, default=8, type=int,
              help="Number of threads to use (default: 8)")

from concurrent.futures import ThreadPoolExecutor
from threading import Lock

def process_read(read, gene, intron, cb, ub, splice_lock, unsplice_lock):
    result = read_overlap(read, gene, intron)
    if result[0] is None:  # Check if gene_id is None
        return None
    
    gene_id, gene_name, splice = result
    if not (read.has_tag(cb) and read.has_tag(ub)):
        return None
        
    cell_barcode = read.get_tag(cb)
    umi = read.get_tag(ub)
    return (cell_barcode, umi, gene_id, gene_name, splice)

def write_result(result, outdir, splice_lock, unsplice_lock):
    cell_barcode, umi, gene_id, gene_name, splice = result
    content = f'{cell_barcode}\t{umi}\t{gene_id}\t{gene_name}\n'
    
    if splice == 'splice':
        with splice_lock:
            with open(f"{outdir}/splice.txt", "a") as f:
                f.write(content)
    else:
        with unsplice_lock:
            with open(f"{outdir}/unsplice.txt", "a") as f:
                f.write(content)

def bam_to_gene(bam: str, outdir: str, genebed: str, introned: str, cb: str = None, ub: str = None, threads: int = 8) -> None:
    """Process bam to gene expression using multi-threading."""
    os.makedirs(f"{outdir}", exist_ok=True)
    
    # Initialize empty output files
    with open(f"{outdir}/splice.txt", "w") as f1, open(f"{outdir}/unsplice.txt", "w") as f2:
        pass
    
    import pysam
    samfile = pysam.AlignmentFile(bam, "rb")
    gene = pd.read_table(genebed, sep="\t", header=None, names=['chrom', 'start', 'end', 'gene_id', 'gene_name', 'strand'])
    intron = pd.read_table(introned, sep="\t", header=None, names=['chrom', 'start', 'end', 'transcript_id', 'gene_id', 'strand'])
    gene['chrom'] = gene['chrom'].astype(str)
    intron['chrom'] = intron['chrom'].astype(str)
    
    splice_lock = Lock()
    unsplice_lock = Lock()
    read_count = 0
    
    # Validate threads parameter
    import multiprocessing
    available_cores = multiprocessing.cpu_count()
    if threads <= 0:
        threads = 1
    elif threads > available_cores * 2:
        print(f"Warning: Reducing threads from {threads} to {available_cores * 2} (2x CPU cores)")
        threads = available_cores * 2
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        for read in samfile.fetch():
            read_count += 1
            if read_count % 10000 == 0:
                print(f"Processed {read_count} reads")
            
            future = executor.submit(
                process_read, 
                read, gene, intron, cb, ub, 
                splice_lock, unsplice_lock
            )
            future.add_done_callback(
                lambda f: write_result(f.result(), outdir, splice_lock, unsplice_lock) if f.result() else None
            )
            futures.append(future)
