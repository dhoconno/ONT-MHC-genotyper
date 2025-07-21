#!/usr/bin/env python3
"""
Concatenate FASTQ files for MHC genotyping pipeline.
Handles both folder input (multiple files) and single file input.
"""

import argparse
import gzip
import os
import sys
from pathlib import Path
import logging

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)

def find_fastq_files(input_path):
    """
    Find all FASTQ.gz files in the input path.
    
    Args:
        input_path (str): Path to folder or file
        
    Returns:
        list: List of FASTQ file paths
    """
    input_path = Path(input_path)
    
    if input_path.is_file():
        if input_path.suffix.lower() in ['.fastq', '.fq'] or \
           (input_path.suffix.lower() == '.gz' and 
            input_path.stem.lower().endswith(('.fastq', '.fq'))):
            return [input_path]
        else:
            raise ValueError(f"Input file {input_path} is not a FASTQ file")
    
    elif input_path.is_dir():
        # Find all FASTQ.gz files in directory
        fastq_files = []
        for pattern in ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']:
            fastq_files.extend(input_path.glob(pattern))
        
        if not fastq_files:
            raise ValueError(f"No FASTQ files found in directory {input_path}")
        
        # Sort files for consistent ordering
        return sorted(fastq_files)
    
    else:
        raise ValueError(f"Input path {input_path} does not exist")

def count_reads_in_file(file_path):
    """
    Count number of reads in a FASTQ file.
    
    Args:
        file_path (Path): Path to FASTQ file
        
    Returns:
        int: Number of reads
    """
    read_count = 0
    
    # Determine if file is gzipped
    if file_path.suffix.lower() == '.gz':
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    try:
        with open_func(file_path, mode) as f:
            for i, line in enumerate(f):
                if i % 4 == 0:  # Every 4th line is a header
                    read_count += 1
    except Exception as e:
        logging.warning(f"Could not count reads in {file_path}: {e}")
        return 0
    
    return read_count

def concatenate_fastq_files(input_files, output_file):
    """
    Concatenate multiple FASTQ files into a single gzipped file.
    
    Args:
        input_files (list): List of input FASTQ file paths
        output_file (str): Path to output concatenated file
        
    Returns:
        dict: Statistics about the concatenation
    """
    logger = logging.getLogger(__name__)
    
    stats = {
        'input_files': len(input_files),
        'total_reads': 0,
        'file_stats': []
    }
    
    logger.info(f"Concatenating {len(input_files)} FASTQ files to {output_file}")
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir:  # Only create directory if there's a directory component
        os.makedirs(output_dir, exist_ok=True)
    
    with gzip.open(output_file, 'wt') as outf:
        for i, input_file in enumerate(input_files):
            logger.info(f"Processing file {i+1}/{len(input_files)}: {input_file}")
            
            # Count reads in this file
            read_count = count_reads_in_file(input_file)
            
            # Determine if input file is gzipped
            if input_file.suffix.lower() == '.gz':
                open_func = gzip.open
                mode = 'rt'
            else:
                open_func = open
                mode = 'r'
            
            # Copy file contents
            try:
                with open_func(input_file, mode) as inf:
                    for line in inf:
                        outf.write(line)
                
                file_stat = {
                    'filename': str(input_file),
                    'reads': read_count,
                    'size_mb': input_file.stat().st_size / (1024*1024)
                }
                stats['file_stats'].append(file_stat)
                stats['total_reads'] += read_count
                
                logger.info(f"  Added {read_count:,} reads from {input_file.name}")
                
            except Exception as e:
                logger.error(f"Error processing {input_file}: {e}")
                raise
    
    logger.info(f"Concatenation complete: {stats['total_reads']:,} total reads")
    return stats

def write_stats_report(stats, output_dir):
    """
    Write concatenation statistics to a report file.
    
    Args:
        stats (dict): Statistics dictionary
        output_dir (str): Output directory for report
    """
    report_file = os.path.join(output_dir, 'concatenation_stats.txt')
    
    with open(report_file, 'w') as f:
        f.write("FASTQ Concatenation Statistics\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Total input files: {stats['input_files']}\n")
        f.write(f"Total reads: {stats['total_reads']:,}\n\n")
        f.write("Per-file statistics:\n")
        f.write("-" * 20 + "\n")
        
        for file_stat in stats['file_stats']:
            f.write(f"File: {os.path.basename(file_stat['filename'])}\n")
            f.write(f"  Reads: {file_stat['reads']:,}\n")
            f.write(f"  Size: {file_stat['size_mb']:.2f} MB\n\n")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Concatenate FASTQ files for MHC genotyping pipeline"
    )
    parser.add_argument(
        'input_path',
        help='Path to folder containing FASTQ.gz files or single FASTQ file'
    )
    parser.add_argument(
        'output_file',
        help='Path to output concatenated FASTQ.gz file'
    )
    parser.add_argument(
        '--stats-dir',
        help='Directory to write statistics report (default: same as output file)',
        default=None
    )
    
    args = parser.parse_args()
    
    # Set up logging
    logger = setup_logging()
    
    try:
        # Find input files
        input_files = find_fastq_files(args.input_path)
        logger.info(f"Found {len(input_files)} FASTQ files to concatenate")
        
        # If only one file and it's already the target format, we could just copy/link
        if len(input_files) == 1:
            input_file = input_files[0]
            if input_file.suffix.lower() == '.gz':
                logger.info("Single gzipped file detected - copying to output location")
                import shutil
                shutil.copy2(input_file, args.output_file)
                
                # Still generate stats
                read_count = count_reads_in_file(input_file)
                stats = {
                    'input_files': 1,
                    'total_reads': read_count,
                    'file_stats': [{
                        'filename': str(input_file),
                        'reads': read_count,
                        'size_mb': input_file.stat().st_size / (1024*1024)
                    }]
                }
            else:
                # Need to gzip the single file
                stats = concatenate_fastq_files(input_files, args.output_file)
        else:
            # Concatenate multiple files
            stats = concatenate_fastq_files(input_files, args.output_file)
        
        # Write statistics report
        stats_dir = args.stats_dir or os.path.dirname(args.output_file)
        write_stats_report(stats, stats_dir)
        
        logger.info("Concatenation completed successfully")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()