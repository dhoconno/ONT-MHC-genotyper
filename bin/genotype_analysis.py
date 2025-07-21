#!/usr/bin/env python3
"""
MHC genotyping analysis for Oxford Nanopore sequencing data.
Adapted from MiSeq genotyping pipeline for single-sample processing.
"""

import argparse
import os
import sys
import subprocess
import tempfile
import logging
import re
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import xlsxwriter

# Import haplotype definitions
from haplotype_definitions import get_species_haplotypes, call_haplotypes

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)

def check_external_tools():
    """Check if required external tools are available."""
    logger = logging.getLogger(__name__)
    tools = ['vsearch', 'bbduk.sh', 'bbmap.sh']
    
    for tool in tools:
        try:
            result = subprocess.run([tool, '--version'], 
                                  capture_output=True, text=True, check=True)
            logger.info(f"Found {tool}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error(f"{tool} not found in PATH. Please install BBTools and vsearch.")
            return False
    return True

def count_reads(fastq_file):
    """
    Count reads in FASTQ file.
    
    Args:
        fastq_file (str): Path to FASTQ file
        
    Returns:
        int: Number of reads
    """
    import gzip
    
    read_count = 0
    
    # Determine if file is gzipped
    if fastq_file.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    with open_func(fastq_file, mode) as f:
        for i, line in enumerate(f):
            if i % 4 == 0:  # Every 4th line is a header
                read_count += 1
    
    return read_count

def remove_primers_ont(reads, primers, out_dir):
    """
    Remove primer sequences from ONT reads using BBDuk.
    Adapted for single-end ONT data.
    
    Args:
        reads (str): Path to input FASTQ file
        primers (str): Path to primer FASTA file or None/NO_FILE to skip
        out_dir (str): Output directory
        
    Returns:
        str: Path to primer-trimmed FASTQ file (or original if no primers)
    """
    logger = logging.getLogger(__name__)
    
    # Check if primers file is provided
    if not primers or primers == 'NO_FILE' or not os.path.exists(primers):
        logger.info("No primers file provided - skipping primer trimming")
        return reads
    
    os.makedirs(out_dir, exist_ok=True)
    
    # Output file names
    base_name = Path(reads).stem.replace('.fastq', '').replace('.fq', '')
    if base_name.endswith('.gz'):
        base_name = base_name[:-3]
    
    # Two-step trimming for ONT reads
    left_trimmed = os.path.join(out_dir, f"{base_name}_l.fastq.gz")
    final_trimmed = os.path.join(out_dir, f"{base_name}_trimmed.fastq.gz")
    
    logger.info("Removing primers from ONT reads...")
    
    # Step 1: Left-end trimming
    cmd_left = [
        'bbduk.sh',
        f'in={reads}',
        f'ref={primers}',
        f'out={left_trimmed}',
        'ktrim=l',
        'k=15',
        'restrictleft=50',  # Larger window for ONT
        'hdist=1',  # Allow 1 mismatch for ONT errors
        'overwrite=true'
    ]
    
    logger.info(f"Left trimming: {' '.join(cmd_left)}")
    result = subprocess.run(cmd_left, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Error in left trimming: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd_left)
    
    # Step 2: Right-end trimming
    cmd_right = [
        'bbduk.sh',
        f'in={left_trimmed}',
        f'ref={primers}',
        f'out={final_trimmed}',
        'ktrim=r',
        'k=15',
        'restrictright=50',  # Larger window for ONT
        'hdist=1',  # Allow 1 mismatch for ONT errors
        'overwrite=true'
    ]
    
    logger.info(f"Right trimming: {' '.join(cmd_right)}")
    result = subprocess.run(cmd_right, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Error in right trimming: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd_right)
    
    # Clean up intermediate file
    os.remove(left_trimmed)
    
    logger.info("Primer trimming completed")
    return final_trimmed

def quality_filter_ont(reads, out_dir, min_length=100, min_quality=7):
    """
    Quality filter ONT reads.
    
    Args:
        reads (str): Path to input FASTQ file
        out_dir (str): Output directory
        min_length (int): Minimum read length
        min_quality (int): Minimum average quality
        
    Returns:
        str: Path to quality-filtered FASTQ file
    """
    logger = logging.getLogger(__name__)
    
    base_name = Path(reads).stem.replace('.fastq', '').replace('.fq', '')
    if base_name.endswith('.gz'):
        base_name = base_name[:-3]
    
    filtered_reads = os.path.join(out_dir, f"{base_name}_filtered.fastq.gz")
    
    logger.info("Quality filtering ONT reads...")
    
    cmd = [
        'bbduk.sh',
        f'in={reads}',
        f'out={filtered_reads}',
        f'minlen={min_length}',
        f'maq={min_quality}',
        'overwrite=true'
    ]
    
    logger.info(f"Quality filtering: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Error in quality filtering: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd)
    
    logger.info("Quality filtering completed")
    return filtered_reads

def vsearch_unique(reads, out_dir):
    """
    Find unique sequences using vsearch.
    
    Args:
        reads (str): Path to input FASTQ file
        out_dir (str): Output directory
        
    Returns:
        str: Path to unique sequences FASTA file
    """
    logger = logging.getLogger(__name__)
    
    base_name = Path(reads).stem.replace('.fastq', '').replace('.fq', '')
    if base_name.endswith('.gz'):
        base_name = base_name[:-3]
    
    unique_fasta = os.path.join(out_dir, f"{base_name}_unique.fasta")
    
    logger.info("Finding unique sequences...")
    
    cmd = [
        'vsearch',
        '--fastx_uniques', reads,
        '--sizeout',
        '--relabel', 'Uniq',
        '--fastaout', unique_fasta
    ]
    
    logger.info(f"Unique sequences: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Error finding unique sequences: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd)
    
    logger.info("Unique sequence identification completed")
    return unique_fasta

def calculate_threshold(min_freq, total_reads):
    """
    Calculate minimum abundance threshold for denoising.
    
    Args:
        min_freq (float): Minimum frequency (default 0.0002 = 0.02%)
        total_reads (int): Total number of reads
        
    Returns:
        int: Minimum read count threshold (minimum 2)
    """
    threshold = max(2, int(total_reads * min_freq))
    return threshold

def vsearch_denoise(reads, total_reads, out_dir):
    """
    Denoise sequences using vsearch UNOISE algorithm.
    
    Args:
        reads (str): Path to unique sequences FASTA file
        total_reads (int): Total number of reads for threshold calculation
        out_dir (str): Output directory
        
    Returns:
        str: Path to denoised sequences FASTA file
    """
    logger = logging.getLogger(__name__)
    
    base_name = Path(reads).stem.replace('.fasta', '').replace('.fa', '')
    
    # Calculate threshold
    min_size = calculate_threshold(0.0002, total_reads)
    logger.info(f"Using minimum abundance threshold: {min_size} reads")
    
    # Output files
    unoise_fasta = os.path.join(out_dir, f"{base_name}_unoise.fasta")
    denoised_fasta = os.path.join(out_dir, f"{base_name}_denoised.fasta")
    
    logger.info("Denoising sequences...")
    
    # Step 1: Cluster denoising
    cmd_unoise = [
        'vsearch',
        '--cluster_unoise', reads,
        '--minsize', str(min_size),
        '--unoise_alpha', '2',
        '--centroids', unoise_fasta
    ]
    
    logger.info(f"UNOISE clustering: {' '.join(cmd_unoise)}")
    result = subprocess.run(cmd_unoise, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Error in UNOISE clustering: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd_unoise)
    
    # Step 2: Chimera removal
    cmd_chimera = [
        'vsearch',
        '--uchime_denovo', unoise_fasta,
        '--abskew', '16',
        '--nonchimeras', denoised_fasta
    ]
    
    logger.info(f"Chimera removal: {' '.join(cmd_chimera)}")
    result = subprocess.run(cmd_chimera, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Error in chimera removal: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd_chimera)
    
    logger.info("Denoising completed")
    return denoised_fasta

def map_to_reference(reads, reference, out_dir):
    """
    Map sequences to reference using BBMap.
    
    Args:
        reads (str): Path to sequences FASTA file
        reference (str): Path to reference FASTA file
        out_dir (str): Output directory
        
    Returns:
        str: Path to SAM output file
    """
    logger = logging.getLogger(__name__)
    
    base_name = Path(reads).stem.replace('.fasta', '').replace('.fa', '')
    sam_output = os.path.join(out_dir, f"{base_name}_mapped.sam")
    
    logger.info("Mapping sequences to reference...")
    
    cmd = [
        'bbmap.sh',
        f'in={reads}',
        f'outm={sam_output}',
        f'ref={reference}',
        'nodisk=t',
        'ambiguous=toss',
        'ordered=t',
        'semiperfectmode=t',
        'overwrite=true'
    ]
    
    logger.info(f"Mapping: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Error in mapping: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd)
    
    logger.info("Mapping completed")
    return sam_output

def parse_sam_file(sam_file, sample_name, total_reads, species):
    """
    Parse SAM file to extract genotype information.
    
    Args:
        sam_file (str): Path to SAM file
        sample_name (str): Sample name
        total_reads (int): Total number of reads
        species (str): Species identifier
        
    Returns:
        pd.DataFrame: Genotyping results
    """
    logger = logging.getLogger(__name__)
    
    logger.info("Parsing SAM file for genotype information...")
    
    # Read SAM file, skipping header lines
    sam_data = []
    with open(sam_file, 'r') as f:
        for line in f:
            if not line.startswith('@'):  # Skip header lines
                parts = line.strip().split('\t')
                # Take only the first 11 columns (standard SAM fields)
                if len(parts) >= 11:
                    sam_data.append(parts[:11])
    
    if not sam_data:
        logger.warning("No alignment data found in SAM file")
        return pd.DataFrame()
    
    # Create DataFrame with explicit string dtypes to avoid accessor issues
    df = pd.DataFrame(sam_data, columns=['query', 'flag', 'reference', 'pos', 'mapq', 
                                        'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'])
    
    # Ensure reference column is string type
    df['reference'] = df['reference'].astype(str)
    
    if df.empty:
        logger.warning("No mapped reads found in SAM file")
        return pd.DataFrame()
    
    # Extract read counts from query names (vsearch size annotations)
    size_pattern = r'size=(\d+)'
    extracted_sizes = df['query'].str.extract(size_pattern)[0]
    
    # Handle cases where size annotation is not found (use 1 as default)
    df['read_count'] = pd.to_numeric(extracted_sizes, errors='coerce').fillna(1).astype(int)
    
    # Group by reference allele and sum read counts
    allele_counts = df.groupby('reference')['read_count'].sum().reset_index()
    allele_counts.columns = ['allele', 'read_count']
    
    # Add metadata
    allele_counts['sample_name'] = sample_name
    allele_counts['species'] = species
    allele_counts['total_reads'] = total_reads
    allele_counts['percentage'] = (allele_counts['read_count'] / total_reads) * 100 if total_reads > 0 else 0
    
    # Parse allele information
    allele_counts['locus'] = allele_counts['allele'].str.extract(r'_([A-Z]+[0-9]*[A-Z]*)_')
    allele_counts['allele_group'] = allele_counts['allele'].str.extract(r'_([A-Z]+[0-9]*)')
    
    logger.info(f"Found {len(allele_counts)} alleles with reads")
    
    return allele_counts

def load_haplotype_definitions():
    """
    Load haplotype definitions for different species.
    This is a simplified version - in practice, these would be loaded from files.
    
    Returns:
        dict: Haplotype definitions by species and locus
    """
    # Simplified haplotype definitions - in practice these would be comprehensive
    haplotype_defs = {
        'Rhesus': {
            'A': {
                'Mamu-A1*001': ['A1*001:01', 'A1*001:02'],
                'Mamu-A1*002': ['A1*002:01', 'A1*002:02'],
                # Add more haplotype definitions...
            },
            'B': {
                'Mamu-B*001': ['B*001:01'],
                'Mamu-B*003': ['B*003:01'],
                # Add more haplotype definitions...
            }
        },
        'Cynomolgus': {
            'A1': {
                'Mafa-A1*001': ['A1*001:01', 'A1*001:02'],
                # Add more haplotype definitions...
            }
        },
        'Pig-tailed': {
            'A': {
                'Mane-A*001': ['A*001:01'],
                # Add more haplotype definitions...
            }
        }
    }
    
    return haplotype_defs

def call_haplotypes_for_sample(genotype_df, species):
    """
    Call haplotypes based on observed alleles using the original MiSeq algorithm.
    
    Args:
        genotype_df (pd.DataFrame): Genotyping results with columns: allele, read_count, sample_name, species, etc.
        species (str): Species identifier ('Rhesus', 'Cynomolgus', 'Pig-tailed')
        
    Returns:
        pd.DataFrame: Genotyping results with haplotype calls added
    """
    logger = logging.getLogger(__name__)
    
    logger.info("Calling haplotypes using original MiSeq algorithm...")
    
    # Get species haplotype definitions
    species_haplotypes = get_species_haplotypes(species)
    prefix = species_haplotypes['PREFIX']
    
    # Initialize haplotype columns
    genotype_df['MHC-A Haplotype 1'] = '-'
    genotype_df['MHC-A Haplotype 2'] = '-'
    genotype_df['MHC-B Haplotype 1'] = '-'
    genotype_df['MHC-B Haplotype 2'] = '-'
    genotype_df['MHC-DRB Haplotype 1'] = '-'
    genotype_df['MHC-DRB Haplotype 2'] = '-'
    genotype_df['MHC-DQA Haplotype 1'] = '-'
    genotype_df['MHC-DQA Haplotype 2'] = '-'
    genotype_df['MHC-DQB Haplotype 1'] = '-'
    genotype_df['MHC-DQB Haplotype 2'] = '-'
    genotype_df['MHC-DPA Haplotype 1'] = '-'
    genotype_df['MHC-DPA Haplotype 2'] = '-'
    genotype_df['MHC-DPB Haplotype 1'] = '-'
    genotype_df['MHC-DPB Haplotype 2'] = '-'
    genotype_df['Comments'] = ''
    
    # Call haplotypes for each locus
    try:
        # MHC-A
        MHC_A = call_haplotypes(locus=str(prefix + '-A'), 
                               locus_haplotype_definitions=species_haplotypes['MHC_A_HAPLOTYPES'], 
                               df=genotype_df)
        genotype_df['MHC-A Haplotype 1'] = MHC_A[0]
        genotype_df['MHC-A Haplotype 2'] = MHC_A[1]
        
        # MHC-B
        MHC_B = call_haplotypes(locus=str(prefix + '-B'), 
                               locus_haplotype_definitions=species_haplotypes['MHC_B_HAPLOTYPES'], 
                               df=genotype_df)
        genotype_df['MHC-B Haplotype 1'] = MHC_B[0]
        genotype_df['MHC-B Haplotype 2'] = MHC_B[1]
        
        # MHC-DRB
        MHC_DRB = call_haplotypes(locus=str(prefix + '-DRB'), 
                                 locus_haplotype_definitions=species_haplotypes['MHC_DRB_HAPLOTYPES'], 
                                 df=genotype_df)
        genotype_df['MHC-DRB Haplotype 1'] = MHC_DRB[0]
        genotype_df['MHC-DRB Haplotype 2'] = MHC_DRB[1]
        
        # MHC-DQA
        MHC_DQA = call_haplotypes(locus=str(prefix + '-DQA'), 
                                 locus_haplotype_definitions=species_haplotypes['MHC_DQA_HAPLOTYPES'], 
                                 df=genotype_df)
        genotype_df['MHC-DQA Haplotype 1'] = MHC_DQA[0]
        genotype_df['MHC-DQA Haplotype 2'] = MHC_DQA[1]
        
        # MHC-DQB
        MHC_DQB = call_haplotypes(locus=str(prefix + '-DQB'), 
                                 locus_haplotype_definitions=species_haplotypes['MHC_DQB_HAPLOTYPES'], 
                                 df=genotype_df)
        genotype_df['MHC-DQB Haplotype 1'] = MHC_DQB[0]
        genotype_df['MHC-DQB Haplotype 2'] = MHC_DQB[1]
        
        # MHC-DPA
        MHC_DPA = call_haplotypes(locus=str(prefix + '-DPA'), 
                                 locus_haplotype_definitions=species_haplotypes['MHC_DPA_HAPLOTYPES'], 
                                 df=genotype_df)
        genotype_df['MHC-DPA Haplotype 1'] = MHC_DPA[0]
        genotype_df['MHC-DPA Haplotype 2'] = MHC_DPA[1]
        
        # MHC-DPB
        MHC_DPB = call_haplotypes(locus=str(prefix + '-DPB'), 
                                 locus_haplotype_definitions=species_haplotypes['MHC_DPB_HAPLOTYPES'], 
                                 df=genotype_df)
        genotype_df['MHC-DPB Haplotype 1'] = MHC_DPB[0]
        genotype_df['MHC-DPB Haplotype 2'] = MHC_DPB[1]
        
        logger.info(f"Haplotype calling completed for {species} ({prefix})")
        
    except Exception as e:
        logger.warning(f"Error during haplotype calling: {e}")
        # Fill with default values if haplotyping fails
        genotype_df['Comments'] = f'Haplotype calling failed: {str(e)}'
    
    return genotype_df

def generate_excel_report(genotype_df, output_file):
    """
    Generate Excel report with genotyping results.
    
    Args:
        genotype_df (pd.DataFrame): Genotyping results
        output_file (str): Path to output Excel file
    """
    logger = logging.getLogger(__name__)
    
    logger.info("Generating Excel report...")
    
    # Create workbook and worksheet with options to handle NaN/Inf values
    workbook = xlsxwriter.Workbook(output_file, {'nan_inf_to_errors': True})
    worksheet = workbook.add_worksheet('Genotyping Results')
    
    # Define formats
    header_format = workbook.add_format({
        'bold': True,
        'bg_color': '#4F81BD',
        'font_color': 'white',
        'border': 1
    })
    
    data_format = workbook.add_format({
        'border': 1,
        'num_format': '0'
    })
    
    percentage_format = workbook.add_format({
        'border': 1,
        'num_format': '0.00%'
    })
    
    # Write headers
    headers = ['Sample', 'Species', 'Locus', 'Allele', 'Read Count', 
               'Percentage', 'Haplotype 1', 'Haplotype 2']
    
    for col, header in enumerate(headers):
        worksheet.write(0, col, header, header_format)
    
    # Write data
    for row, (_, data) in enumerate(genotype_df.iterrows(), 1):
        worksheet.write(row, 0, data.get('sample_name', ''), data_format)
        worksheet.write(row, 1, data.get('species', ''), data_format)
        worksheet.write(row, 2, data.get('locus', ''), data_format)
        worksheet.write(row, 3, data.get('allele', ''), data_format)
        worksheet.write(row, 4, data.get('read_count', 0), data_format)
        worksheet.write(row, 5, data.get('percentage', 0) / 100, percentage_format)
        worksheet.write(row, 6, data.get('haplotype_1', ''), data_format)
        worksheet.write(row, 7, data.get('haplotype_2', ''), data_format)
    
    # Auto-adjust column widths
    for col in range(len(headers)):
        worksheet.set_column(col, col, 15)
    
    # Freeze header row
    worksheet.freeze_panes(1, 0)
    
    workbook.close()
    
    logger.info(f"Excel report generated: {output_file}")

def write_stats_report(stats, output_dir, sample_name=None):
    """Write processing statistics to report file."""
    if sample_name:
        report_file = os.path.join(output_dir, f'{sample_name}_genotyping_stats.txt')
    else:
        report_file = os.path.join(output_dir, 'genotyping_stats.txt')
    
    with open(report_file, 'w') as f:
        f.write("MHC Genotyping Statistics\n")
        f.write("=" * 40 + "\n\n")
        
        for key, value in stats.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:.2f}\n")
            else:
                f.write(f"{key}: {value}\n")

def process_sample(input_fastq, reference_fasta, primers_fasta, sample_name, 
                  species, output_dir):
    """
    Process a single sample through the complete genotyping pipeline.
    
    Args:
        input_fastq (str): Path to input FASTQ file
        reference_fasta (str): Path to reference sequences
        primers_fasta (str): Path to primer sequences
        sample_name (str): Sample identifier
        species (str): Species identifier
        output_dir (str): Output directory
        
    Returns:
        pd.DataFrame: Genotyping results
    """
    logger = logging.getLogger(__name__)
    
    logger.info(f"Processing sample: {sample_name}")
    
    # Create intermediate directories
    work_dir = os.path.join(output_dir, 'intermediate')
    os.makedirs(work_dir, exist_ok=True)
    
    stats = {'sample_name': sample_name, 'species': species}
    
    try:
        # Count initial reads
        initial_reads = count_reads(input_fastq)
        stats['initial_reads'] = initial_reads
        logger.info(f"Initial reads: {initial_reads}")
        
        # Step 1: Remove primers
        trimmed_fastq = remove_primers_ont(input_fastq, primers_fasta, work_dir)
        trimmed_reads = count_reads(trimmed_fastq)
        stats['trimmed_reads'] = trimmed_reads
        logger.info(f"Reads after primer trimming: {trimmed_reads}")
        
        # Step 2: Quality filtering
        filtered_fastq = quality_filter_ont(trimmed_fastq, work_dir)
        filtered_reads = count_reads(filtered_fastq)
        stats['filtered_reads'] = filtered_reads
        logger.info(f"Reads after quality filtering: {filtered_reads}")
        
        # Step 3: Find unique sequences
        unique_fasta = vsearch_unique(filtered_fastq, work_dir)
        
        # Step 4: Denoise sequences
        denoised_fasta = vsearch_denoise(unique_fasta, filtered_reads, work_dir)
        
        # Step 5: Map to reference
        sam_file = map_to_reference(denoised_fasta, reference_fasta, work_dir)
        
        # Step 6: Parse results
        genotype_df = parse_sam_file(sam_file, sample_name, filtered_reads, species)
        
        if not genotype_df.empty:
            # Step 7: Call haplotypes
            genotype_df = call_haplotypes_for_sample(genotype_df, species)
            
            # Calculate mapping statistics
            mapped_reads = genotype_df['read_count'].sum()
            stats['mapped_reads'] = mapped_reads
            stats['mapping_rate'] = mapped_reads / filtered_reads if filtered_reads > 0 else 0
            stats['unique_alleles'] = len(genotype_df)
            
            logger.info(f"Mapped reads: {mapped_reads} ({stats['mapping_rate']:.2%})")
            logger.info(f"Unique alleles detected: {stats['unique_alleles']}")
        else:
            logger.warning("No alleles detected for this sample")
            stats['mapped_reads'] = 0
            stats['mapping_rate'] = 0
            stats['unique_alleles'] = 0
        
        # Write statistics
        write_stats_report(stats, output_dir, sample_name)
        
        return genotype_df
        
    except Exception as e:
        logger.error(f"Error processing sample {sample_name}: {e}")
        raise

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="MHC genotyping analysis for Oxford Nanopore data"
    )
    parser.add_argument(
        'input_fastq',
        help='Path to input FASTQ.gz file'
    )
    parser.add_argument(
        'reference_fasta',
        help='Path to MHC reference sequences'
    )
    parser.add_argument(
        'primers_fasta',
        help='Path to primer sequences'
    )
    parser.add_argument(
        'output_dir',
        help='Output directory for results'
    )
    parser.add_argument(
        '--sample-name',
        required=True,
        help='Sample name'
    )
    parser.add_argument(
        '--species',
        required=True,
        choices=['Rhesus', 'Cynomolgus', 'Pig-tailed'],
        help='Species identifier'
    )
    
    args = parser.parse_args()
    
    # Set up logging
    logger = setup_logging()
    
    try:
        # Check external tools
        if not check_external_tools():
            sys.exit(1)
        
        # Ensure output directory exists
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Process sample
        genotype_df = process_sample(
            args.input_fastq,
            args.reference_fasta,
            args.primers_fasta,
            args.sample_name,
            args.species,
            args.output_dir
        )
        
        # Generate Excel report with sample name prefix
        if not genotype_df.empty:
            excel_file = os.path.join(args.output_dir, f'{args.sample_name}_genotyping_results.xlsx')
            generate_excel_report(genotype_df, excel_file)
            
            # Save CSV version with sample name prefix
            csv_file = os.path.join(args.output_dir, f'{args.sample_name}_genotyping_results.csv')
            genotype_df.to_csv(csv_file, index=False)
        
        logger.info("Genotyping analysis completed successfully")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()