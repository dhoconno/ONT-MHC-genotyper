# Test Data

This directory should contain small test datasets for pipeline validation.

## Required Test Files

To run the test pipeline, you need:

1. **Barcode FASTQ files** (2-3 small files)
   - Example: `barcode01.fastq.gz`, `barcode02.fastq.gz`
   - Can be subsampled from real data

2. **Sample mapping file** (`test_sample_mapping.csv`):
   ```csv
   tag,GS ID
   BC01,Test_Sample_001
   BC02,Test_Sample_002
   ```

3. **Reference FASTA** (small subset)
   - A few MHC allele sequences for testing

## Creating Test Data

From existing data:
```bash
# Subsample FASTQ files (first 1000 reads)
seqkit head -n 1000 original_barcode01.fastq.gz | gzip > test_data/barcode01.fastq.gz

# Create small reference with few alleles
seqkit grep -n -p "Mamu-A1*001:01,Mamu-B*001:01" original_reference.fasta > test_data/test_reference.fasta
```

## Running Tests

```bash
nextflow run workflow/mhc_genotyping.nf \
  --barcode_dir test_data/ \
  --reference test_data/test_reference.fasta \
  --sample_sheet test_data/test_sample_mapping.csv \
  --outdir test_results \
  -profile test
```