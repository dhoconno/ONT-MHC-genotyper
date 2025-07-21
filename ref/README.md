# Reference Files

This directory contains reference files required for the MHC genotyping pipeline.

## Required Files

1. **Reference FASTA** (`reference.fasta`)
   - MHC allele sequences in FASTA format
   - User must provide their own reference file
   - Example reference files available:
     - Macaque: `Mafa_MiSeq-IPD_17.06.01_2.2.0.0_plus_SP.fasta`
     - Custom: User-generated reference sequences

2. **Primer sequences** (`mhc_specific_primers.fa`)
   - Forward and reverse primers used in amplicon generation
   - Default primers included for standard MHC amplicons

3. **Fluidigm barcodes** (`fluidigm.txt`)
   - Tab-separated file with barcode names and sequences
   - Default Fluidigm barcode set included

4. **Sample mapping template** (`sample_mapping_template.csv`)
   - Example CSV format for mapping barcodes to sample names
   - Copy and modify for your samples

## Adding Your Reference

1. Copy your MHC reference FASTA to this directory
2. Update the pipeline command to use your reference:
   ```bash
   nextflow run workflow/mhc_genotyping.nf \
     --reference ref/your_reference.fasta \
     --barcode_dir input/
   ```

## Reference Format

The reference FASTA should follow standard MHC nomenclature:
```
>Mamu-A1*001:01
ATGCGGGTCACGGCGCCCCGAACCCTCCTCCTGCTGCTCTCGGCGGCCCTGGCCCTGACCGAGACCTGGGCC...
>Mamu-B*001:01
ATGCGGGTCACGGCACCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCGAGACCTGGGCG...
```