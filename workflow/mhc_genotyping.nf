#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * MHC Genotyping Pipeline
 * This pipeline performs MHC genotyping from MiSeq sequencing data using Fluidigm barcodes
 */

// Pipeline parameters
params.barcode_dir = null
params.reference = "${projectDir}/ref/31588_MCM_MHC_miseq_deduplicated_annotated.fasta"
params.primers = "${projectDir}/ref/mhc_specific_primers.fa"
params.fluidigm_barcodes = "${projectDir}/ref/fluidigm.txt"
params.sample_sheet = "${projectDir}/ref/sample_mapping.csv"
params.outdir = "${projectDir}/results"
params.min_reads = 10
params.mismatch = 2  // Allow up to 2 mismatches in primer matching

/*
 * Process: Concatenate barcode files
 */
process CONCATENATE_BARCODES {
    publishDir "${params.outdir}/01_concatenated", mode: 'copy', enabled: false
    
    input:
    path barcode_files
    
    output:
    path "concatenated.fq.gz", emit: concatenated
    path "concatenation_stats.txt", emit: stats
    
    script:
    """
    cat ${barcode_files} > concatenated.fq.gz
    
    # Generate stats
    echo "Concatenation Statistics" > concatenation_stats.txt
    echo "========================" >> concatenation_stats.txt
    echo "Input files: \$(ls -1 ${barcode_files} | wc -l)" >> concatenation_stats.txt
    echo "Total reads: \$(zcat concatenated.fq.gz | wc -l | awk '{print \$1/4}')" >> concatenation_stats.txt
    """
}

/*
 * Process: Orient reads to reference
 */
process ORIENT_READS {
    publishDir "${params.outdir}/02_oriented", mode: 'copy', enabled: false
    
    input:
    path concatenated
    path reference
    
    output:
    path "oriented.fq.gz", emit: oriented
    path "orientation_stats.txt", emit: stats
    
    script:
    """
    # Orient reads
    vsearch --orient ${concatenated} \
        --db ${reference} \
        --fastqout oriented.fq \
        --log orientation_stats.txt \
        --threads ${task.cpus}
    
    # Compress
    pigz -c oriented.fq > oriented.fq.gz
    rm oriented.fq
    """
}

/*
 * Process: Prepare barcode information
 */
process PREPARE_BARCODES {
    input:
    path barcode_file
    
    output:
    path "barcode_info.tsv", emit: barcode_info
    
    script:
    """
    # Create barcode info file without header
    tail -n +2 ${barcode_file} > barcode_info.tsv
    """
}

/*
 * Process: Demultiplex individual barcode
 */
process DEMULTIPLEX_SINGLE {
    tag "${barcode_name}"
    publishDir "${params.outdir}/03_demultiplexed", mode: 'copy', pattern: "${barcode_name}.fq.gz", enabled: false
    
    input:
    tuple val(barcode_name), val(barcode_seq)
    path oriented
    
    output:
    tuple val(barcode_name), path("${barcode_name}.fq.gz"), emit: demuxed
    path "${barcode_name}_stats.txt", emit: stats
    
    script:
    """
    # Extract reads for this specific barcode
    seqkit grep -j ${task.cpus} -s -p "${barcode_seq}" ${oriented} -o ${barcode_name}.fq.gz 2>/dev/null || touch ${barcode_name}.fq.gz
    
    # Count reads
    if [[ -s ${barcode_name}.fq.gz ]]; then
        read_count=\$(seqkit stats ${barcode_name}.fq.gz | tail -n 1 | awk '{print \$4}' | tr -d ',')
    else
        read_count=0
    fi
    
    # Write stats
    echo "${barcode_name},${barcode_seq},\${read_count}" > ${barcode_name}_stats.txt
    """
}

/*
 * Process: Filter for forward primers
 */
process FILTER_FORWARD_PRIMER {
    tag "${sample_id}"
    publishDir "${params.outdir}/04_filtered/${sample_id}", mode: 'copy', enabled: false
    
    input:
    tuple val(sample_id), path(reads)
    path primers
    
    output:
    tuple val(sample_id), path("${sample_id}_with_fwd.fq.gz"), emit: filtered
    path "${sample_id}_fwd_stats.txt", emit: stats
    
    script:
    """
    # Filter for reads WITH forward primers
    bbduk.sh in=${reads} \
        outm=${sample_id}_with_fwd.fq.gz \
        ref=${primers} \
        k=15 \
        hdist=${params.mismatch} \
        restrictleft=100 \
        copyundefined=t \
        stats=${sample_id}_fwd_stats.txt
    
    # Check if output has reads
    if [[ ! -s ${sample_id}_with_fwd.fq.gz ]]; then
        echo "No reads with forward primer found" > ${sample_id}_fwd_stats.txt
    fi
    """
}

/*
 * Process: Trim forward primers
 */
process TRIM_FORWARD_PRIMER {
    tag "${sample_id}"
    publishDir "${params.outdir}/05_fwd_trimmed/${sample_id}", mode: 'copy', enabled: false
    
    input:
    tuple val(sample_id), path(reads)
    path primers
    
    output:
    tuple val(sample_id), path("${sample_id}_fwd_trimmed.fq.gz"), emit: trimmed
    
    script:
    """
    # Trim forward primers
    bbduk.sh in=${reads} \
        out=${sample_id}_fwd_trimmed.fq.gz \
        ref=${primers} \
        k=15 \
        ktrim=l \
        hdist=${params.mismatch} \
        restrictleft=100 \
        copyundefined=t
    """
}

/*
 * Process: Filter for reverse primers
 */
process FILTER_REVERSE_PRIMER {
    tag "${sample_id}"
    publishDir "${params.outdir}/06_rev_filtered/${sample_id}", mode: 'copy', enabled: false
    
    input:
    tuple val(sample_id), path(reads)
    path primers
    
    output:
    tuple val(sample_id), path("${sample_id}_with_both.fq.gz"), emit: filtered
    path "${sample_id}_rev_stats.txt", emit: stats
    
    script:
    """
    # Filter for reads WITH reverse primers
    bbduk.sh in=${reads} \
        outm=${sample_id}_with_both.fq.gz \
        ref=${primers} \
        k=15 \
        hdist=${params.mismatch} \
        restrictright=100 \
        rcomp=t \
        copyundefined=t \
        stats=${sample_id}_rev_stats.txt
    
    # Check if output has reads
    if [[ ! -s ${sample_id}_with_both.fq.gz ]]; then
        echo "No reads with reverse primer found" > ${sample_id}_rev_stats.txt
    fi
    """
}

/*
 * Process: Trim reverse primers
 */
process TRIM_REVERSE_PRIMER {
    tag "${sample_id}"
    publishDir "${params.outdir}/07_trimmed/${sample_id}", mode: 'copy', enabled: false
    
    input:
    tuple val(sample_id), path(reads)
    path primers
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fq.gz"), emit: trimmed
    path "${sample_id}_trim_summary.txt", emit: summary
    
    script:
    """
    # Trim reverse primers
    bbduk.sh in=${reads} \
        out=${sample_id}_trimmed.fq.gz \
        ref=${primers} \
        k=15 \
        ktrim=r \
        hdist=${params.mismatch} \
        restrictright=100 \
        rcomp=t \
        copyundefined=t
    
    # Create summary
    echo "Sample: ${sample_id}" > ${sample_id}_trim_summary.txt
    echo "Final trimmed reads: \$(seqkit stats ${sample_id}_trimmed.fq.gz | tail -n 1 | awk '{print \$4}')" >> ${sample_id}_trim_summary.txt
    """
}

/*
 * Process: Map reads to reference
 */
process MAP_READS {
    tag "${sample_id}"
    publishDir "${params.outdir}/08_mapped/${sample_id}", mode: 'copy', enabled: false
    
    input:
    tuple val(sample_id), path(trimmed)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}_mapped.sam"), emit: mapped
    path "${sample_id}_mapping_stats.txt", emit: stats
    
    script:
    """
    # Map reads in semi-perfect mode
    mapPacBio.sh in=${trimmed} \
        outm=${sample_id}_mapped.sam \
        ref=${reference} \
        semiperfectmode=t \
        threads=${task.cpus} \
        statsfile=${sample_id}_mapping_stats.txt
    """
}

/*
 * Process: Filter full-span alignments
 */
process FILTER_FULLSPAN {
    tag "${sample_id}"
    publishDir "${params.outdir}/09_fullspan/${sample_id}", mode: 'copy', enabled: false
    
    input:
    tuple val(sample_id), path(mapped)
    
    output:
    tuple val(sample_id), path("${sample_id}_fullspan.sam"), emit: fullspan
    path "${sample_id}_fullspan_stats.txt", emit: stats
    
    script:
    """
    python ${projectDir}/../bin/filter_fullspan.py \
        ${mapped} \
        ${sample_id}_fullspan.sam \
        --stats > ${sample_id}_fullspan_stats.txt 2>&1
    """
}

/*
 * Process: Count alleles
 */
process COUNT_ALLELES {
    tag "${sample_id}"
    publishDir "${params.outdir}/10_counts/${sample_id}", mode: 'copy', enabled: false
    
    input:
    tuple val(sample_id), path(fullspan)
    
    output:
    path "${sample_id}_allele_counts.csv", emit: counts
    
    script:
    """
    # Count reads per allele
    samtools view -F 4 ${fullspan} | \
        cut -f3 | \
        sort | \
        uniq -c | \
        while read count allele; do 
            # Properly quote the allele name to handle embedded commas
            echo "${sample_id},\"\${allele}\",\${count}"
        done > ${sample_id}_allele_counts.csv
    """
}

/*
 * Process: Aggregate all counts and rename samples
 */
process AGGREGATE_COUNTS {
    publishDir "${params.outdir}/11_aggregated", mode: 'copy'
    
    input:
    path count_files
    path sample_sheet
    
    output:
    path "all_samples_counts.csv", emit: aggregated
    path "sample_summary.txt", emit: summary
    
    script:
    """
    # Aggregate all counts (samples already renamed)
    echo "sample,allele,count" > all_samples_counts.csv
    
    # Simply concatenate all count files
    for file in ${count_files}; do
        cat \$file >> all_samples_counts.csv
    done
    
    # Generate summary
    echo "Sample Summary" > sample_summary.txt
    echo "==============" >> sample_summary.txt
    echo "Total samples processed: \$(ls -1 ${count_files} | wc -l)" >> sample_summary.txt
    echo "" >> sample_summary.txt
    echo "Reads per sample:" >> sample_summary.txt
    tail -n +2 all_samples_counts.csv | \
        awk -F',' '{sum[\$1]+=\$3} END {for (s in sum) print s ": " sum[s] " reads"}' | \
        sort >> sample_summary.txt
    """
}

/*
 * Process: Create pivot table
 */
process CREATE_PIVOT_TABLE {
    publishDir "${params.outdir}/12_pivot_table", mode: 'copy'
    
    input:
    path aggregated_counts
    
    output:
    path "allele_counts_pivot.xlsx", emit: pivot_table
    path "pivot_summary.txt", emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import sys
    import re

    # Read the aggregated counts manually to handle complex allele names
    data = []
    with open('${aggregated_counts}', 'r') as f:
        header = f.readline().strip()  # Skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Parse each line carefully
            # Format is: sample,allele,count where allele may contain commas
            parts = line.split(',')
            if len(parts) >= 3:
                sample = parts[0]
                count = parts[-1]  # Last element is always the count
                # Everything in between is the allele name
                allele = ','.join(parts[1:-1])
                # Remove quotes if present
                allele = allele.strip('"')
                try:
                    count = int(count)
                    data.append({'sample': sample, 'allele': allele, 'count': count})
                except ValueError:
                    print(f"Skipping line with invalid count: {line}")
    
    # Create DataFrame
    df = pd.DataFrame(data)
    
    if df.empty:
        print("No data found in the aggregated counts file")
        sys.exit(1)
    
    # Create pivot table with samples as columns and alleles as rows
    pivot_df = df.pivot_table(
        index='allele', 
        columns='sample', 
        values='count', 
        fill_value=0,
        aggfunc='sum'
    )
    
    # Sort by total count across all samples
    pivot_df['Total'] = pivot_df.sum(axis=1)
    pivot_df = pivot_df.sort_values('Total', ascending=False)
    
    # Write to Excel
    with pd.ExcelWriter('allele_counts_pivot.xlsx', engine='xlsxwriter') as writer:
        pivot_df.to_excel(writer, sheet_name='Allele Counts')
        
        # Get workbook and worksheet
        workbook = writer.book
        worksheet = writer.sheets['Allele Counts']
        
        # Format the header row
        header_format = workbook.add_format({
            'bold': True,
            'text_wrap': True,
            'valign': 'top',
            'fg_color': '#D7E4BD',
            'border': 1
        })
        
        # Write headers with formatting
        for col_num, value in enumerate(pivot_df.columns):
            worksheet.write(0, col_num + 1, value, header_format)
        
        # Adjust column widths
        worksheet.set_column(0, 0, 80)  # Allele column (wider for long names)
        worksheet.set_column(1, len(pivot_df.columns), 12)  # Sample columns
    
    # Write summary
    with open('pivot_summary.txt', 'w') as f:
        f.write(f"Pivot Table Summary\\n")
        f.write(f"==================\\n")
        f.write(f"Total samples: {len(df['sample'].unique())}\\n")
        f.write(f"Total alleles: {len(df['allele'].unique())}\\n")
        f.write(f"Total reads: {df['count'].sum():,}\\n")
        f.write(f"\\nTop 10 alleles by total count:\\n")
        for allele, count in pivot_df['Total'].head(10).items():
            f.write(f"  {allele}: {count:,}\\n")
    
    print("Pivot table created successfully")
    """
}

/*
 * Main workflow
 */
workflow {
    // Validate required parameters
    if (!params.barcode_dir) {
        error "Please provide a barcode directory using --barcode_dir"
    }

    // Print parameters
    log.info """
    MHC Genotyping Pipeline
    =======================
    barcode_dir      : ${params.barcode_dir}
    reference        : ${params.reference}
    primers          : ${params.primers}
    fluidigm_barcodes: ${params.fluidigm_barcodes}
    sample_sheet     : ${params.sample_sheet}
    outdir           : ${params.outdir}
    min_reads        : ${params.min_reads}
    mismatch         : ${params.mismatch}
    """

    // Get barcode files
    barcode_files = Channel
        .fromPath("${params.barcode_dir}/*.fastq.gz")
        .collect()
    
    // Concatenate barcodes
    CONCATENATE_BARCODES(barcode_files)
    
    // Orient reads
    ORIENT_READS(
        CONCATENATE_BARCODES.out.concatenated,
        file(params.reference)
    )
    
    // Prepare barcodes
    PREPARE_BARCODES(file(params.fluidigm_barcodes))
    
    // Create channel of barcode info
    barcode_ch = PREPARE_BARCODES.out.barcode_info
        .splitCsv(sep: '\t', header: false)
        .map { row -> tuple(row[0], row[1]) }
    
    // Demultiplex individual barcodes
    DEMULTIPLEX_SINGLE(
        barcode_ch,
        ORIENT_READS.out.oriented
    )
    
    // Filter out empty files and rename samples
    sample_mapping = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row['tag'], row['GS ID']) }
    
    // Join demultiplexed files with sample names
    demux_renamed = DEMULTIPLEX_SINGLE.out.demuxed
        .filter { barcode_name, file -> file.size() > 0 }
        .join(sample_mapping, by: 0, remainder: true)
        .map { barcode_name, file, sample_name ->
            def final_name = sample_name ?: barcode_name
            tuple(final_name, file)
        }
        .filter { sample_name, file -> file != null }
    
    // Filter for forward primers
    FILTER_FORWARD_PRIMER(
        demux_renamed,
        file(params.primers)
    )
    
    // Trim forward primers
    TRIM_FORWARD_PRIMER(
        FILTER_FORWARD_PRIMER.out.filtered.filter { sample_id, file -> file.size() > 0 },
        file(params.primers)
    )
    
    // Filter for reverse primers
    FILTER_REVERSE_PRIMER(
        TRIM_FORWARD_PRIMER.out.trimmed,
        file(params.primers)
    )
    
    // Trim reverse primers
    TRIM_REVERSE_PRIMER(
        FILTER_REVERSE_PRIMER.out.filtered.filter { sample_id, file -> file.size() > 0 },
        file(params.primers)
    )
    
    // Map reads
    MAP_READS(
        TRIM_REVERSE_PRIMER.out.trimmed,
        file(params.reference)
    )
    
    // Filter full-span
    FILTER_FULLSPAN(MAP_READS.out.mapped)
    
    // Count alleles
    COUNT_ALLELES(FILTER_FULLSPAN.out.fullspan)
    
    // Aggregate counts
    AGGREGATE_COUNTS(
        COUNT_ALLELES.out.counts.collect(),
        file(params.sample_sheet)
    )
    
    // Create pivot table
    CREATE_PIVOT_TABLE(AGGREGATE_COUNTS.out.aggregated)
}

/*
 * Completion handler
 */
workflow.onComplete = {
    log.info """
    Pipeline completed!
    ===================
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    Results dir : ${params.outdir}
    """
}