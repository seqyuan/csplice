# csplice

Genomic data processing toolkit for converting GTF to BED format and analyzing BAM gene overlaps.

## Installation

```bash
# Install using poetry
poetry install

# Or install directly
pip install .
```

## Usage

### gtf2bed command
Convert GTF to BED files:
```bash
csplice gtf2bed -g input.gtf -o output_dir
```

This will generate 4 BED files in the output directory:
- gene.bed: Standard 6-column BED format (chrom, start, end, gene_id, gene_name, strand)
- transcript.bed: Standard 6-column BED format (chrom, start, end, transcript_id, gene_id, strand) 
- exon.bed: Standard 6-column BED format (chrom, start, end, transcript_id, gene_id, strand)
- intron.bed: Standard 6-column BED format (chrom, start, end, transcript_id, gene_id, strand)

### bam2gene command
Analyze BAM file gene overlaps:
```bash
csplice bam2gene -b input.bam -g genes.bed -o output_dir
```

This will generate 2 files in the output directory:
- splice.txt: TSV format with columns (barcode, umi, gene_id, gene_name), reads with out intron
- unsplice.txt: TSV format with columns (barcode, umi, gene_id, gene_name), reads with intron

## Options

### gtf2bed
```
Options:
  -g, --gtf TEXT         Input GTF file path  [required]
  -o, --outdir TEXT      Output directory  [required]
  -i, --gene_id TEXT     Gene ID key in attributes  [default: gene_id]
  -n, --gene_name TEXT   Gene name key in attributes  [default: gene_name]
  -t, --transcript_id TEXT  Transcript ID key in attributes  [default: transcript_id]
  --help                 Show this message and exit.
```

### bam2gene
```
Options:
  -b, --bam TEXT         Input BAM file path  [required]
  -o, --outdir TEXT      Output directory  [required]
  -g, --gene_bed TEXT    Gene BED file path  [required]
  --help                 Show this message and exit.
```

## Development

```bash
# Install dev dependencies
poetry install --with dev

# Run tests
poetry run pytest
