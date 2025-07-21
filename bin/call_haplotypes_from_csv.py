#!/usr/bin/env python3
"""
Call MHC haplotypes from allele count CSV files and generate Excel reports.
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import logging
from datetime import datetime

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


def read_allele_counts(csv_file):
    """
    Read allele counts from CSV file.
    
    Expected format: sample,allele,count
    """
    try:
        df = pd.read_csv(csv_file, names=['sample', 'allele', 'count'])
        # Clean allele names - remove quotes if present
        df['allele'] = df['allele'].str.strip('"')
        
        # Extract the simplified allele name from complex format
        # e.g., "05_M1M2M3_A1_063g" from "01_Mafa_A1_063g|A1_063_01,_A1_063_02..."
        df['allele_simple'] = df['allele'].str.extract(r'^\d+_(?:Mafa_|Mamu_)?(.+?)(?:\||$)', expand=False)
        
        # For haplotype matching, we'll use the simplified name
        df['allele_for_haplotype'] = df['allele_simple'].fillna(df['allele'])
        
        return df
    except Exception as e:
        logging.error(f"Error reading CSV file: {e}")
        return None


def process_genotypes(df, species='Rhesus', min_reads=10):
    """
    Process genotypes and call haplotypes.
    
    Args:
        df: DataFrame with allele counts
        species: Species name for haplotype definitions
        min_reads: Minimum read count threshold
        
    Returns:
        DataFrame with genotypes and haplotypes
    """
    logger = logging.getLogger(__name__)
    
    # Get species-specific haplotype definitions
    species_haps = get_species_haplotypes(species)
    prefix = species_haps['PREFIX']
    
    # Filter by minimum read count
    df_filtered = df[df['count'] >= min_reads].copy()
    
    if len(df_filtered) == 0:
        logger.warning(f"No alleles with >= {min_reads} reads")
        return None
    
    # Add total reads per sample
    total_reads = df.groupby('sample')['count'].sum().to_dict()
    df_filtered['total_reads'] = df_filtered['sample'].map(total_reads)
    
    # Calculate percentage
    df_filtered['percent_reads'] = (df_filtered['count'] / df_filtered['total_reads'] * 100).round(2)
    
    # Process each sample
    results = []
    
    for sample in df_filtered['sample'].unique():
        sample_df = df_filtered[df_filtered['sample'] == sample].copy()
        
        # For haplotype calling, use simplified allele names
        sample_df_for_hap = sample_df.copy()
        sample_df_for_hap['allele'] = sample_df_for_hap['allele_for_haplotype']
        
        # Add haplotype calls for each locus
        sample_dict = {
            'sample': sample,
            'total_reads': sample_df['total_reads'].iloc[0],
            'num_alleles': len(sample_df)
        }
        
        # Call haplotypes for each locus
        loci = ['A', 'B', 'DRB', 'DQA', 'DQB', 'DPA', 'DPB']
        
        for locus in loci:
            locus_name = f'{prefix}-{locus}'
            hap_defs = species_haps.get(f'MHC_{locus}_HAPLOTYPES', {})
            
            if hap_defs:
                h1, h2 = call_haplotypes(locus_name, hap_defs, sample_df_for_hap)
                sample_dict[f'MHC-{locus} Haplotype 1'] = h1
                sample_dict[f'MHC-{locus} Haplotype 2'] = h2
        
        results.append(sample_dict)
    
    # Create results dataframe
    results_df = pd.DataFrame(results)
    
    # Merge with allele data
    genotypes_df = df_filtered.merge(
        results_df[['sample'] + [col for col in results_df.columns if 'Haplotype' in col]], 
        on='sample'
    )
    
    return genotypes_df


def create_pivot_table(df):
    """
    Create pivot table from genotyping DataFrame.
    """
    # Add locus identifier column - extract from simplified allele name
    df['locus'] = df['allele_simple'].str.extract(r'^([A-Z]+\d*)', expand=False)
    
    # If locus extraction failed, try from original allele name
    df['locus'] = df['locus'].fillna(df['allele'].str.extract(r'^[A-Z]+\d*_([A-Z]+)', expand=False))
    
    # Select columns for pivot
    haplotype_cols = [col for col in df.columns if 'Haplotype' in col]
    
    # Check if we have data to pivot
    if len(df) == 0:
        logger = logging.getLogger(__name__)
        logger.warning("No data available for pivot table")
        return pd.DataFrame()
    
    # Create pivot table - simpler version
    try:
        pivot_df = pd.pivot_table(
            df,
            index=['allele'],
            columns=['sample'],
            values='count',
            fill_value=0
        )
    except Exception as e:
        logger = logging.getLogger(__name__)
        logger.warning(f"Could not create pivot table: {e}")
        pivot_df = pd.DataFrame()
    
    return pivot_df


def generate_excel_report(genotypes_df, pivot_df, output_file, sample_info=None):
    """
    Generate Excel report with genotypes and pivot table.
    """
    logger = logging.getLogger(__name__)
    
    try:
        with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
            # Get workbook
            workbook = writer.book
            
            # Add formats
            header_format = workbook.add_format({
                'bold': True,
                'bg_color': '#4472C4',
                'font_color': 'white',
                'border': 1
            })
            
            cell_format = workbook.add_format({
                'border': 1
            })
            
            # Sheet 1: Summary
            summary_data = []
            for sample in genotypes_df['sample'].unique():
                sample_df = genotypes_df[genotypes_df['sample'] == sample]
                
                summary_dict = {
                    'Sample': sample,
                    'Total Reads': sample_df['total_reads'].iloc[0],
                    'Alleles Detected': len(sample_df),
                }
                
                # Add haplotypes
                hap_cols = [col for col in sample_df.columns if 'Haplotype' in col]
                for col in hap_cols:
                    summary_dict[col] = sample_df[col].iloc[0]
                
                summary_data.append(summary_dict)
            
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_excel(writer, sheet_name='Summary', index=False)
            
            # Format summary sheet
            worksheet = writer.sheets['Summary']
            for col_num, value in enumerate(summary_df.columns.values):
                worksheet.write(0, col_num, value, header_format)
                
            # Sheet 2: Genotypes
            genotypes_df.to_excel(writer, sheet_name='Genotypes', index=False)
            
            # Sheet 3: Pivot Table
            if len(pivot_df) > 0:
                pivot_df.to_excel(writer, sheet_name='Pivot Table')
            else:
                # Write empty dataframe with message
                pd.DataFrame({'Message': ['No pivot table data available']}).to_excel(
                    writer, sheet_name='Pivot Table', index=False
                )
            
            # Sheet 4: Sample Info (if provided)
            if sample_info is not None:
                sample_info.to_excel(writer, sheet_name='Sample Info', index=False)
            
            logger.info(f"Excel report generated: {output_file}")
            
    except Exception as e:
        logger.error(f"Error generating Excel report: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Call MHC haplotypes from allele count CSV files"
    )
    parser.add_argument(
        "csv_file",
        help="CSV file with allele counts (format: sample,allele,count)"
    )
    parser.add_argument(
        "-o", "--output",
        default="haplotype_report.xlsx",
        help="Output Excel file name (default: haplotype_report.xlsx)"
    )
    parser.add_argument(
        "-s", "--species",
        choices=['Rhesus', 'Cynomolgus', 'Pig-tailed'],
        default='Rhesus',
        help="Species for haplotype definitions (default: Rhesus)"
    )
    parser.add_argument(
        "-m", "--min-reads",
        type=int,
        default=10,
        help="Minimum read count threshold (default: 10)"
    )
    parser.add_argument(
        "--sample-info",
        help="Optional CSV file with sample metadata"
    )
    
    args = parser.parse_args()
    
    # Set up logging
    logger = setup_logging()
    
    # Check input file
    if not Path(args.csv_file).exists():
        logger.error(f"Input file not found: {args.csv_file}")
        sys.exit(1)
    
    # Read allele counts
    logger.info(f"Reading allele counts from {args.csv_file}")
    df = read_allele_counts(args.csv_file)
    
    if df is None or len(df) == 0:
        logger.error("No data found in CSV file")
        sys.exit(1)
    
    # Process genotypes and call haplotypes
    logger.info(f"Processing genotypes for {args.species}")
    genotypes_df = process_genotypes(df, species=args.species, min_reads=args.min_reads)
    
    if genotypes_df is None or len(genotypes_df) == 0:
        logger.error("No genotypes after filtering")
        sys.exit(1)
    
    # Create pivot table
    logger.info("Creating pivot table")
    pivot_df = create_pivot_table(genotypes_df)
    
    # Load sample info if provided
    sample_info = None
    if args.sample_info:
        try:
            sample_info = pd.read_csv(args.sample_info)
        except Exception as e:
            logger.warning(f"Could not read sample info: {e}")
    
    # Generate Excel report
    logger.info(f"Generating Excel report: {args.output}")
    generate_excel_report(genotypes_df, pivot_df, args.output, sample_info)
    
    # Print summary
    print(f"\nHaplotype calling complete!")
    print(f"Samples processed: {len(genotypes_df['sample'].unique())}")
    print(f"Total alleles: {len(genotypes_df)}")
    print(f"Output file: {args.output}")


if __name__ == "__main__":
    main()