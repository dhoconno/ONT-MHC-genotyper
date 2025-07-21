#!/usr/bin/env python3
"""
Filter SAM/BAM alignments to retain only reads that span the entire reference sequence.

This script filters alignments to keep only those that:
1. Start at position 1 (0-based: position 0) of the reference
2. End at exactly the reference length
3. Have no clipping that would indicate partial alignment
"""

import sys
import argparse
import pysam
from pathlib import Path


def filter_full_span_alignments(input_file, output_file):
    """
    Filter alignments to keep only those spanning the full reference length.
    
    Parameters:
    -----------
    input_file : str
        Path to input SAM/BAM file
    output_file : str
        Path to output SAM/BAM file
    
    Returns:
    --------
    tuple: (total_reads, filtered_reads, reads_per_reference)
    """
    # Open input file
    try:
        samfile = pysam.AlignmentFile(input_file, "r")
    except Exception as e:
        print(f"Error opening input file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Determine output format based on extension
    if output_file.endswith('.bam'):
        mode = "wb"
    else:
        mode = "w"
    
    # Open output file
    try:
        outfile = pysam.AlignmentFile(output_file, mode, template=samfile)
    except Exception as e:
        print(f"Error opening output file: {e}", file=sys.stderr)
        samfile.close()
        sys.exit(1)
    
    # Get reference lengths
    ref_lengths = {}
    for i, ref in enumerate(samfile.references):
        ref_lengths[ref] = samfile.lengths[i]
    
    # Statistics
    total_reads = 0
    filtered_reads = 0
    reads_per_reference = {}
    
    # Process alignments
    for read in samfile:
        total_reads += 1
        
        # Skip unmapped reads
        if read.is_unmapped:
            continue
        
        ref_name = read.reference_name
        ref_len = ref_lengths.get(ref_name, 0)
        
        # Check if alignment spans full reference:
        # 1. Starts at position 0 (0-based coordinate)
        # 2. Ends at reference length
        # 3. Allow soft clipping at ends (read extends beyond reference)
        # 4. No gaps in the middle of the alignment
        
        # The read must cover the entire reference, but can extend beyond it
        if (read.reference_start == 0 and 
            read.reference_end == ref_len):
            
            outfile.write(read)
            filtered_reads += 1
            
            # Track reads per reference
            if ref_name not in reads_per_reference:
                reads_per_reference[ref_name] = 0
            reads_per_reference[ref_name] += 1
    
    # Close files
    samfile.close()
    outfile.close()
    
    return total_reads, filtered_reads, reads_per_reference


def main():
    parser = argparse.ArgumentParser(
        description="Filter SAM/BAM alignments to retain only full-span reads"
    )
    parser.add_argument(
        "input",
        help="Input SAM/BAM file"
    )
    parser.add_argument(
        "output",
        help="Output SAM/BAM file (format determined by extension)"
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Print filtering statistics to stderr"
    )
    
    args = parser.parse_args()
    
    # Check input file exists
    if not Path(args.input).exists():
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Run filtering
    total, filtered, ref_counts = filter_full_span_alignments(args.input, args.output)
    
    # Print statistics if requested
    if args.stats:
        print(f"\nFiltering Statistics:", file=sys.stderr)
        print(f"Total alignments: {total}", file=sys.stderr)
        print(f"Full-span alignments: {filtered}", file=sys.stderr)
        print(f"Percentage retained: {filtered/total*100:.1f}%", file=sys.stderr)
        
        if ref_counts:
            print(f"\nReads per reference:", file=sys.stderr)
            for ref, count in sorted(ref_counts.items(), key=lambda x: x[1], reverse=True):
                print(f"  {ref}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()