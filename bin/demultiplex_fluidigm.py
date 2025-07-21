#!/usr/bin/env python3
"""
Demultiplex reads by Fluidigm tags with exact sequence matching.
Searches for tags near 5' and 3' ends with zero error tolerance.
"""

import argparse
import gzip
import os
import sys
import logging
from pathlib import Path
from collections import defaultdict
import pandas as pd

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)

def load_fluidigm_tags(tag_file):
    """
    Load Fluidigm tags from file.
    
    Args:
        tag_file (str): Path to Fluidigm tags file
        
    Returns:
        dict: Mapping of tag sequence to tag name
    """
    logger = logging.getLogger(__name__)
    tags = {}
    
    logger.info(f"Loading Fluidigm tags from {tag_file}")
    
    with open(tag_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Handle different file formats
            if '\t' in line:
                parts = line.split('\t')
            elif ',' in line:
                parts = line.split(',')
            else:
                parts = line.split()
            
            if len(parts) >= 2:
                tag_name = parts[0].strip()
                tag_sequence = parts[1].strip().upper()
                
                # Remove any arrow notation from tag names (like "1→FLD0001")
                if '→' in tag_name:
                    tag_name = tag_name.split('→')[1]
                
                tags[tag_sequence] = tag_name
                logger.debug(f"Loaded tag {tag_name}: {tag_sequence}")
            else:
                logger.warning(f"Skipping malformed line {line_num}: {line}")
    
    logger.info(f"Loaded {len(tags)} Fluidigm tags")
    return tags

def reverse_complement(sequence):
    """Generate reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def search_for_tags(sequence, tags, search_window=50):
    """
    Search for Fluidigm tags in a sequence.
    
    Args:
        sequence (str): DNA sequence to search
        tags (dict): Dictionary of tag sequences to tag names
        search_window (int): Window size to search from ends
        
    Returns:
        tuple: (found_tag_name, found_position, found_sequence) or (None, None, None)
    """
    sequence = sequence.upper().strip()
    
    # Define search regions
    start_region = sequence[:search_window] if len(sequence) > search_window else sequence
    end_region = sequence[-search_window:] if len(sequence) > search_window else sequence
    
    # Search in start region (5' end)
    for tag_seq, tag_name in tags.items():
        # Forward orientation
        if tag_seq in start_region:
            pos = start_region.find(tag_seq)
            return tag_name, pos, tag_seq
        
        # Reverse complement
        rc_tag = reverse_complement(tag_seq)
        if rc_tag in start_region:
            pos = start_region.find(rc_tag)
            return tag_name, pos, rc_tag
    
    # Search in end region (3' end)
    for tag_seq, tag_name in tags.items():
        # Forward orientation
        if tag_seq in end_region:
            pos = len(sequence) - search_window + end_region.find(tag_seq)
            return tag_name, pos, tag_seq
        
        # Reverse complement
        rc_tag = reverse_complement(tag_seq)
        if rc_tag in end_region:
            pos = len(sequence) - search_window + end_region.find(rc_tag)
            return tag_name, pos, rc_tag
    
    return None, None, None

def demultiplex_fastq(input_fastq, output_dir, tags, sample_mapping=None):
    """
    Demultiplex FASTQ file by Fluidigm tags.
    
    Args:
        input_fastq (str): Path to input FASTQ file
        output_dir (str): Output directory for demultiplexed files
        tags (dict): Dictionary of tag sequences to tag names
        sample_mapping (dict): Optional mapping of tag names to sample names
        
    Returns:
        dict: Statistics about demultiplexing
    """
    logger = logging.getLogger(__name__)
    
    logger.info(f"Demultiplexing {input_fastq}")
    logger.info(f"Output directory: {output_dir}")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Statistics tracking
    stats = {
        'total_reads': 0,
        'tagged_reads': 0,
        'untagged_reads': 0,
        'tag_counts': defaultdict(int),
        'tag_rate': 0.0
    }
    
    # Open output files
    output_files = {}
    untagged_file = gzip.open(os.path.join(output_dir, 'untagged.fastq.gz'), 'wt')
    
    # Determine if input is gzipped
    if input_fastq.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    try:
        with open_func(input_fastq, mode) as inf:
            while True:
                # Read FASTQ record (4 lines)
                header = inf.readline()
                if not header:
                    break
                
                sequence = inf.readline()
                plus = inf.readline()
                quality = inf.readline()
                
                if not (header and sequence and plus and quality):
                    break
                
                stats['total_reads'] += 1
                
                # Search for tags
                tag_name, tag_pos, found_seq = search_for_tags(sequence.strip(), tags)
                
                if tag_name:
                    # Tagged read
                    stats['tagged_reads'] += 1
                    stats['tag_counts'][tag_name] += 1
                    
                    # Get sample name from mapping if available
                    sample_name = sample_mapping.get(tag_name, tag_name) if sample_mapping else tag_name
                    
                    # Open output file for this tag if not already open
                    if tag_name not in output_files:
                        output_filename = f"{sample_name}.fastq.gz"
                        output_path = os.path.join(output_dir, output_filename)
                        output_files[tag_name] = gzip.open(output_path, 'wt')
                    
                    # Write to tagged output
                    output_files[tag_name].write(header)
                    output_files[tag_name].write(sequence)
                    output_files[tag_name].write(plus)
                    output_files[tag_name].write(quality)
                    
                else:
                    # Untagged read
                    stats['untagged_reads'] += 1
                    untagged_file.write(header)
                    untagged_file.write(sequence)
                    untagged_file.write(plus)
                    untagged_file.write(quality)
                
                # Progress reporting
                if stats['total_reads'] % 10000 == 0:
                    logger.info(f"Processed {stats['total_reads']:,} reads...")
    
    finally:
        # Close all output files
        untagged_file.close()
        for f in output_files.values():
            f.close()
    
    # Calculate final statistics
    if stats['total_reads'] > 0:
        stats['tag_rate'] = stats['tagged_reads'] / stats['total_reads']
    
    logger.info(f"Demultiplexing complete:")
    logger.info(f"  Total reads: {stats['total_reads']:,}")
    logger.info(f"  Tagged reads: {stats['tagged_reads']:,} ({stats['tag_rate']:.2%})")
    logger.info(f"  Untagged reads: {stats['untagged_reads']:,}")
    logger.info(f"  Unique tags found: {len(stats['tag_counts'])}")
    
    return stats

def load_sample_mapping(mapping_file):
    """
    Load sample mapping from Excel file.
    
    Args:
        mapping_file (str): Path to Excel file with tag-to-sample mapping
        
    Returns:
        dict: Mapping of tag names to sample names
    """
    logger = logging.getLogger(__name__)
    
    if not mapping_file or not os.path.exists(mapping_file):
        logger.info("No sample mapping file provided or file not found")
        return {}
    
    logger.info(f"Loading sample mapping from {mapping_file}")
    
    try:
        # Try to read Excel file
        df = pd.read_excel(mapping_file)
        
        # Look for columns that might contain tag and sample information
        tag_col = None
        sample_col = None
        
        for col in df.columns:
            col_lower = col.lower()
            if 'tag' in col_lower or 'barcode' in col_lower or 'fld' in col_lower:
                tag_col = col
            elif 'sample' in col_lower or 'name' in col_lower:
                sample_col = col
        
        if tag_col is None or sample_col is None:
            logger.warning(f"Could not identify tag and sample columns in {mapping_file}")
            logger.info(f"Available columns: {list(df.columns)}")
            return {}
        
        # Create mapping
        mapping = {}
        for _, row in df.iterrows():
            tag = str(row[tag_col]).strip()
            sample = str(row[sample_col]).strip()
            if tag and sample and tag != 'nan' and sample != 'nan':
                mapping[tag] = sample
        
        logger.info(f"Loaded {len(mapping)} sample mappings")
        return mapping
        
    except Exception as e:
        logger.error(f"Error loading sample mapping: {e}")
        return {}

def write_stats_report(stats, output_dir, tags):
    """Write demultiplexing statistics to report files."""
    
    # Write summary statistics
    summary_file = os.path.join(output_dir, 'demultiplexing_stats.txt')
    with open(summary_file, 'w') as f:
        f.write("Fluidigm Demultiplexing Statistics\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Total reads processed: {stats['total_reads']:,}\n")
        f.write(f"Successfully tagged: {stats['tagged_reads']:,}\n")
        f.write(f"Untagged reads: {stats['untagged_reads']:,}\n")
        f.write(f"Tagging rate: {stats['tag_rate']:.2%}\n")
        f.write(f"Unique tags detected: {len(stats['tag_counts'])}\n\n")
        
        f.write("Per-tag statistics:\n")
        f.write("-" * 20 + "\n")
        for tag_name in sorted(stats['tag_counts'].keys()):
            count = stats['tag_counts'][tag_name]
            percentage = count / stats['total_reads'] if stats['total_reads'] > 0 else 0
            f.write(f"{tag_name}: {count:,} reads ({percentage:.2%})\n")
    
    # Write detailed CSV report
    csv_file = os.path.join(output_dir, 'demultiplexing_details.csv')
    tag_data = []
    
    for tag_seq, tag_name in tags.items():
        count = stats['tag_counts'].get(tag_name, 0)
        percentage = count / stats['total_reads'] if stats['total_reads'] > 0 else 0
        
        tag_data.append({
            'tag_name': tag_name,
            'tag_sequence': tag_seq,
            'read_count': count,
            'percentage': percentage,
            'detected': count > 0
        })
    
    df = pd.DataFrame(tag_data)
    df.to_csv(csv_file, index=False)

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Demultiplex reads by Fluidigm tags with exact matching"
    )
    parser.add_argument(
        'input_fastq',
        help='Path to input oriented FASTQ.gz file'
    )
    parser.add_argument(
        'tag_file',
        help='Path to file containing Fluidigm tag sequences'
    )
    parser.add_argument(
        'output_dir',
        help='Output directory for demultiplexed FASTQ files'
    )
    parser.add_argument(
        '--sample-mapping',
        help='Excel file with tag-to-sample name mapping',
        default=None
    )
    parser.add_argument(
        '--search-window',
        type=int,
        default=50,
        help='Window size to search for tags from read ends (default: 50)'
    )
    
    args = parser.parse_args()
    
    # Set up logging
    logger = setup_logging()
    
    try:
        # Load Fluidigm tags
        tags = load_fluidigm_tags(args.tag_file)
        if not tags:
            logger.error("No tags loaded from tag file")
            sys.exit(1)
        
        # Load sample mapping if provided
        sample_mapping = load_sample_mapping(args.sample_mapping)
        
        # Demultiplex reads
        stats = demultiplex_fastq(
            args.input_fastq,
            args.output_dir,
            tags,
            sample_mapping
        )
        
        # Write statistics reports
        write_stats_report(stats, args.output_dir, tags)
        
        logger.info("Demultiplexing completed successfully")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()