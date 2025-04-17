# csplice

Genomic data processing toolkit for converting GTF to BED format.

## Installation

```bash
# Install using poetry
poetry install

# Or install directly
pip install .
```

## Usage

Convert GTF to BED files:

```bash
csplice gtf2bed -g input.gtf -o output_dir
```

This will generate 4 BED files in the output directory:
- gene.bed
- transcript.bed 
- exon.bed
- intron.bed

## Options

```
Options:
  -g, --gtf TEXT         Input GTF file path  [required]
  -o, --outdir TEXT      Output directory  [required]
  -i, --gene_id TEXT     Gene ID key in attributes  [default: gene_id]
  -n, --gene_name TEXT   Gene name key in attributes  [default: gene_name]
  -t, --transcript_id TEXT  Transcript ID key in attributes  [default: transcript_id]
  --help                 Show this message and exit.
```

## Development

```bash
# Install dev dependencies
poetry install --with dev

# Run tests
poetry run pytest
