# ONT-MHC-genotyper

A Nextflow pipeline for MHC genotyping from Oxford Nanopore sequencing data using Fluidigm barcode-based demultiplexing.

## Overview

This pipeline processes MiSeq amplicon data sequenced on Oxford Nanopore Technologies (ONT) platforms to determine MHC alleles. It performs:

- Barcode-based demultiplexing using Fluidigm tags
- Read orientation and quality control
- Primer filtering and trimming
- Reference-based allele calling
- Generation of allele count matrices

## Features

- **Automated demultiplexing**: Uses Fluidigm barcode sequences for sample identification
- **Primer-aware processing**: Filters and trims reads based on MHC-specific primers
- **Full-span alignment**: Ensures only complete amplicon sequences are counted
- **Flexible output**: Generates both detailed per-sample results and aggregated pivot tables
- **Multiple environment support**: Run locally, on HPC clusters, or in containers

## Requirements

### System Requirements
- Unix-like operating system (Linux, macOS)
- Nextflow (>=21.04.0)
- Python (>=3.9)
- 8GB RAM minimum (more recommended for large datasets)
- 4 CPU cores recommended

### Software Dependencies
The pipeline requires the following tools, which can be installed via conda/mamba or pixi:
- vsearch (>=2.21.0)
- bbmap (>=39.00)
- seqkit (>=2.10.0)
- samtools (>=1.17)
- pigz
- Python packages: pandas, pysam, xlsxwriter, biopython, openpyxl, numpy

## Installation

### Option 1: Using Pixi (Recommended)

```bash
# Clone the repository
git clone https://github.com/dholab/ONT-MHC-genotyper.git
cd ONT-MHC-genotyper

# Install pixi if not already installed
curl -fsSL https://pixi.sh/install.sh | bash

# Install dependencies
pixi install
```

### Option 2: Using Conda/Mamba

```bash
# Clone the repository
git clone https://github.com/dholab/ONT-MHC-genotyper.git
cd ONT-MHC-genotyper

# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate mhc-genotyping
```

### Option 3: Using Docker

```bash
# Pull the Docker image (once available)
docker pull dholab/ont-mhc-genotyper:latest
```

## Quick Start

1. **Prepare your input files:**
   - Barcode FASTQ files in a single directory
   - Sample mapping CSV file
   - Reference FASTA file
   - (Optional) Custom primer sequences

2. **Create a sample mapping file** (`sample_mapping.csv`):
   ```csv
   tag,GS ID
   BC01,Sample_001
   BC02,Sample_002
   BC03,Sample_003
   ```

3. **Run the pipeline:**
   ```bash
   nextflow run workflow/mhc_genotyping.nf \
     --barcode_dir /path/to/barcode/files \
     --reference /path/to/reference.fasta \
     --sample_sheet sample_mapping.csv \
     --outdir results
   ```

## Pipeline Parameters

### Required Parameters
- `--barcode_dir`: Directory containing barcode FASTQ files (*.fastq.gz)
- `--reference`: Path to reference FASTA file containing MHC allele sequences
- `--sample_sheet`: CSV file mapping barcode tags to sample names

### Optional Parameters
- `--outdir`: Output directory (default: `results`)
- `--primers`: FASTA file with primer sequences (default: `ref/mhc_specific_primers.fa`)
- `--fluidigm_barcodes`: Fluidigm barcode file (default: `ref/fluidigm.txt`)
- `--min_reads`: Minimum reads for allele calling (default: 10)
- `--mismatch`: Maximum mismatches allowed in primer matching (default: 2)

### Resource Parameters
- `--max_memory`: Maximum memory (default: '8.GB')
- `--max_cpus`: Maximum CPUs (default: 4)
- `--max_time`: Maximum time (default: '12.h')

## Input File Formats

### Barcode FASTQ Files
- Standard FASTQ format (can be gzipped)
- Named by barcode (e.g., `barcode01.fastq.gz`)

### Sample Mapping CSV
```csv
tag,GS ID
BC01,Sample_001
BC02,Sample_002
```

### Reference FASTA
Standard FASTA format with MHC allele sequences:
```
>Mamu-A1*001:01
ATGCGGGTCACGGCGCCCCGAACCCTCCTCCTGCTGCTCTCGGCGGCCCTGGCCCTGACCGAGACCTGGGCCGGCTCCCACTCCATGAGGTATTTCTCCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCGGGAGACACGGAACGCCAAGGGCCACGCACAGACTGACCGAGAGAACCTGCGGATCGCGCTCCGCTACTACAACCAGAGCGAGGCCGGGTCTCACACCCTCCAGAGGATGTACGGCTGCGACGTGGGGCCGGACGGGCGCCTCCTCCGCGGGCATGACCAGTCCGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCCTGGACCGCCGCGGACACGGCGGCTCAGATCACCCAGCGCAAGTTGGAGGCGGCCCGTGCGGCGGAGCAGCTGAGAGCCTACCTGGAGGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGGAACACCCAAAGACACACGTGACCCACCACCCCCTCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGCGAGGACCAAACTCAGGACACCGAGCTTGTGGAGACCAGGCCAGCAGGAGATGGAACCTTCCAGAAGTGGGCAGCTGTGGTGGTGCCTTCTGGAGAAGAGCAGAGATACACGTGCCATGTGCAGCACGAGGGGCTGCCGGAGCCCCTCACCCTGAGATGGGAGCCGTCTTCCCAGTCCACCGTCCCCATCGTGGGCATTGTTGCTGGCCTGGCTGTCCTAGCAGTTGTGGTCATCGGAGCTGTGGTCGCTGCTGTGATGTGTAGGAGGAAGAGCTCAGG
```

### Fluidigm Barcode File
Tab-separated file with barcode names and sequences:
```
tag	tag_seq
BC01	AAGAAAGTTGTCGGTGTCTTTGTG
BC02	TCGATTCCGTTTGTAGTCGTCTGT
```

## Output Files

The pipeline generates organized output in the following structure:

```
results/
├── 11_aggregated/
│   ├── all_samples_counts.csv    # Aggregated allele counts
│   └── sample_summary.txt         # Summary statistics
├── 12_pivot_table/
│   ├── allele_counts_pivot.xlsx  # Excel pivot table
│   └── pivot_summary.txt          # Pivot table summary
└── pipeline_info/
    ├── execution_report.html      # Detailed execution report
    ├── execution_timeline.html    # Visual timeline
    └── execution_trace.txt        # Resource usage trace
```

### Key Output Files

1. **all_samples_counts.csv**: CSV file with allele counts per sample
   ```csv
   sample,allele,count
   Sample_001,Mamu-A1*001:01,1523
   Sample_001,Mamu-B*001:01,892
   ```

2. **allele_counts_pivot.xlsx**: Excel spreadsheet with samples as columns and alleles as rows

3. **execution_report.html**: Interactive HTML report with pipeline statistics

## Advanced Usage

### Running with Different Profiles

```bash
# Run with Docker
nextflow run workflow/mhc_genotyping.nf -profile docker --barcode_dir input/

# Run on SLURM cluster
nextflow run workflow/mhc_genotyping.nf -profile slurm --barcode_dir input/

# Debug mode with detailed logging
nextflow run workflow/mhc_genotyping.nf -profile debug --barcode_dir input/
```

### Custom Reference Files

To use your own reference sequences:

1. Create a FASTA file with your MHC allele sequences
2. Ensure sequence names follow standard nomenclature
3. Specify with `--reference` parameter

### Adjusting Filtering Parameters

```bash
# Stricter primer matching (1 mismatch)
nextflow run workflow/mhc_genotyping.nf \
  --barcode_dir input/ \
  --mismatch 1

# Higher minimum read threshold
nextflow run workflow/mhc_genotyping.nf \
  --barcode_dir input/ \
  --min_reads 50
```

## Troubleshooting

### Common Issues

1. **No reads passing filter**
   - Check primer sequences match your amplicons
   - Verify read orientation
   - Try increasing `--mismatch` parameter

2. **Memory errors**
   - Increase `--max_memory` parameter
   - Use `-profile slurm` for cluster execution

3. **Missing alleles**
   - Verify reference FASTA contains expected sequences
   - Check minimum read threshold

### Getting Help

- Check the [execution_report.html] for detailed error messages
- Review the [execution_trace.txt] for resource usage
- Open an issue on GitHub with the error message and trace file

## Citation

If you use this pipeline in your research, please cite:

```
ONT-MHC-genotyper: A Nextflow pipeline for MHC genotyping from Oxford Nanopore sequencing data
[Citation details to be added upon publication]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Contact

For questions or support, please open an issue on GitHub or contact the DHO Lab.