#!/usr/bin/env python3
"""
Orient reads using vsearch against MHC reference sequences.
This step ensures all reads are in the same orientation before demultiplexing.
"""

import argparse
import os
import sys
import subprocess
import tempfile
import logging
from pathlib import Path

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)

def check_vsearch():
    """Check if vsearch is available in PATH."""
    try:
        result = subprocess.run(['vsearch', '--version'], 
                              capture_output=True, text=True, check=True)
        logging.info(f"Found vsearch: {result.stderr.split()[0]}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.error("vsearch not found in PATH. Please install vsearch.")
        return False

def orient_reads_vsearch(input_fastq, reference_fasta, output_fastq, temp_dir=None):
    """
    Orient reads using vsearch against reference sequences.
    
    Args:
        input_fastq (str): Path to input FASTQ file
        reference_fasta (str): Path to reference FASTA file
        output_fastq (str): Path to output oriented FASTQ file
        temp_dir (str): Temporary directory for intermediate files
        
    Returns:
        dict: Statistics about the orientation process
    """
    logger = logging.getLogger(__name__)
    
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp(prefix='orient_reads_')
    else:
        os.makedirs(temp_dir, exist_ok=True)
    
    logger.info("Starting read orientation with vsearch")
    logger.info(f"Input FASTQ: {input_fastq}")
    logger.info(f"Reference: {reference_fasta}")
    logger.info(f"Output FASTQ: {output_fastq}")
    
    # Create intermediate file names
    oriented_fasta = os.path.join(temp_dir, "oriented.fasta")
    oriented_fastq_temp = os.path.join(temp_dir, "oriented_temp.fastq")
    
    try:
        # Step 1: Convert FASTQ to FASTA for vsearch orientation
        logger.info("Converting FASTQ to FASTA for orientation...")
        
        # Use vsearch to convert FASTQ to FASTA
        cmd_fastq_to_fasta = [
            'vsearch',
            '--fastq_filter', input_fastq,
            '--fastaout', os.path.join(temp_dir, 'input.fasta'),
            '--fastq_qmax', '50'  # Allow high quality scores
        ]
        
        logger.info(f"Running: {' '.join(cmd_fastq_to_fasta)}")
        result = subprocess.run(cmd_fastq_to_fasta, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"Error converting FASTQ to FASTA: {result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, cmd_fastq_to_fasta)
        
        input_fasta = os.path.join(temp_dir, 'input.fasta')
        
        # Step 2: Orient sequences using vsearch
        logger.info("Orienting sequences against reference...")
        
        cmd_orient = [
            'vsearch',
            '--orient', input_fasta,
            '--db', reference_fasta,
            '--fastaout', oriented_fasta,
            '--notmatched', os.path.join(temp_dir, 'notmatched.fasta'),
            '--threads', '4'
        ]
        
        logger.info(f"Running: {' '.join(cmd_orient)}")
        result = subprocess.run(cmd_orient, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"Error orienting sequences: {result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, cmd_orient)
        
        # Parse vsearch output for statistics
        stats = parse_vsearch_orient_output(result.stderr)
        
        # Step 3: Convert oriented FASTA back to FASTQ with original quality scores
        logger.info("Converting oriented sequences back to FASTQ...")
        
        # Read original FASTQ to create quality mapping
        oriented_reads = read_fasta_headers(oriented_fasta)
        quality_map = create_quality_mapping(input_fastq, oriented_reads)
        
        # Create oriented FASTQ file
        create_oriented_fastq(oriented_fasta, quality_map, output_fastq)
        
        logger.info(f"Orientation complete. Statistics: {stats}")
        
        return stats
        
    except Exception as e:
        logger.error(f"Error during orientation: {e}")
        raise
    finally:
        # Clean up temporary files if we created the temp directory
        if temp_dir.startswith('/tmp/') or temp_dir.startswith(tempfile.gettempdir()):
            import shutil
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                logger.warning(f"Could not clean up temp directory {temp_dir}: {e}")

def parse_vsearch_orient_output(stderr_output):
    """Parse vsearch orient output to extract statistics."""
    stats = {
        'total_sequences': 0,
        'oriented_sequences': 0,
        'not_matched': 0,
        'orientation_rate': 0.0
    }
    
    for line in stderr_output.split('\n'):
        if 'sequences' in line.lower():
            # Parse lines like "12345 sequences processed"
            parts = line.strip().split()
            if len(parts) >= 2 and parts[1] == 'sequences':
                try:
                    count = int(parts[0])
                    if 'processed' in line:
                        stats['total_sequences'] = count
                    elif 'matched' in line and 'not' not in line:
                        stats['oriented_sequences'] = count
                except ValueError:
                    continue
    
    # Calculate derived statistics
    if stats['total_sequences'] > 0:
        stats['not_matched'] = stats['total_sequences'] - stats['oriented_sequences']
        stats['orientation_rate'] = stats['oriented_sequences'] / stats['total_sequences']
    
    return stats

def read_fasta_headers(fasta_file):
    """Read FASTA file and return set of sequence headers."""
    headers = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Remove '>' and any additional annotations
                header = line[1:].split()[0]
                headers.add(header)
    return headers

def create_quality_mapping(input_fastq, oriented_headers):
    """Create mapping of sequence ID to quality scores from original FASTQ."""
    import gzip
    
    quality_map = {}
    
    # Determine if file is gzipped
    if input_fastq.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    with open_func(input_fastq, mode) as f:
        while True:
            header_line = f.readline()
            if not header_line:
                break
            
            seq_line = f.readline()
            plus_line = f.readline() 
            qual_line = f.readline()
            
            if not (header_line and seq_line and plus_line and qual_line):
                break
            
            # Extract sequence ID (remove @ and any additional info)
            seq_id = header_line[1:].split()[0]
            
            # Only keep quality for sequences that were oriented
            if seq_id in oriented_headers:
                quality_map[seq_id] = {
                    'sequence': seq_line.strip(),
                    'quality': qual_line.strip()
                }
    
    return quality_map

def create_oriented_fastq(oriented_fasta, quality_map, output_fastq):
    """Create FASTQ file from oriented FASTA using original quality scores."""
    import gzip
    
    # Output will be gzipped
    with gzip.open(output_fastq, 'wt') as outf:
        with open(oriented_fasta, 'r') as inf:
            seq_id = None
            sequence = ""
            
            for line in inf:
                line = line.strip()
                if line.startswith('>'):
                    # Process previous sequence if exists
                    if seq_id and sequence:
                        write_fastq_record(outf, seq_id, sequence, quality_map)
                    
                    # Start new sequence
                    seq_id = line[1:].split()[0]
                    sequence = ""
                else:
                    sequence += line
            
            # Process last sequence
            if seq_id and sequence:
                write_fastq_record(outf, seq_id, sequence, quality_map)

def write_fastq_record(outf, seq_id, sequence, quality_map):
    """Write a single FASTQ record."""
    if seq_id in quality_map:
        quality = quality_map[seq_id]['quality']
        # Ensure quality string matches sequence length
        if len(quality) != len(sequence):
            # If lengths don't match, create dummy quality
            quality = 'I' * len(sequence)  # High quality placeholder
    else:
        # Create dummy quality if original not found
        quality = 'I' * len(sequence)
    
    outf.write(f"@{seq_id}\n")
    outf.write(f"{sequence}\n")
    outf.write(f"+\n")
    outf.write(f"{quality}\n")

def write_stats_report(stats, output_dir):
    """Write orientation statistics to a report file."""
    report_file = os.path.join(output_dir, 'orientation_stats.txt')
    
    with open(report_file, 'w') as f:
        f.write("Read Orientation Statistics\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Total sequences processed: {stats['total_sequences']:,}\n")
        f.write(f"Successfully oriented: {stats['oriented_sequences']:,}\n")
        f.write(f"Not matched to reference: {stats['not_matched']:,}\n")
        f.write(f"Orientation rate: {stats['orientation_rate']:.2%}\n")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Orient reads using vsearch against MHC reference"
    )
    parser.add_argument(
        'input_fastq',
        help='Path to input FASTQ.gz file'
    )
    parser.add_argument(
        'reference_fasta',
        help='Path to reference FASTA file for orientation'
    )
    parser.add_argument(
        'output_fastq',
        help='Path to output oriented FASTQ.gz file'
    )
    parser.add_argument(
        '--temp-dir',
        help='Temporary directory for intermediate files',
        default=None
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
        # Check if vsearch is available
        if not check_vsearch():
            sys.exit(1)
        
        # Ensure output directory exists
        output_dir = os.path.dirname(args.output_fastq)
        if output_dir:  # Only create directory if there's a directory component
            os.makedirs(output_dir, exist_ok=True)
        
        # Orient reads
        stats = orient_reads_vsearch(
            args.input_fastq,
            args.reference_fasta, 
            args.output_fastq,
            args.temp_dir
        )
        
        # Write statistics report
        stats_dir = args.stats_dir or os.path.dirname(args.output_fastq)
        write_stats_report(stats, stats_dir)
        
        logger.info("Read orientation completed successfully")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()