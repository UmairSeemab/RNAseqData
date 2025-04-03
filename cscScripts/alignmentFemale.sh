#!/bin/bash
# Directory containing STAR genome index
GENOME_DIR="/scratch/project_2003924/bulkData_20250327/dataFiles/genome/human/STAR"
# Directory where your FASTQ files are stored
FASTQ_DIR="/scratch/project_2003924/bulkData_20250327/fastQC_files"
# Directory to store the output
OUTPUT_DIR="/scratch/project_2003924/bulkData_20250327/STARoutput20250331"
# STAR executable path (if not in PATH)
#STAR="/path/to/STAR"
# Number of threads to use
THREADS=16
# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR
# Loop through all paired-end read files in FASTQ directory
# Assuming file names are like sample1_1.fq.gz and sample1_2.fq.gz
# Loop through all forward read files in FASTQ directory
# Assuming file names are like sample_1_f.fastq.gz for forward reads
for R1 in $FASTQ_DIR/*_1_f.fastq.gz; do
    # Replace '_1_f' with '_2_f' to find the matching reverse file
    R2=${R1/_1_f.fastq.gz/_2_f.fastq.gz}
    
    # Extract the sample name based on the file naming convention
    sample_name=$(basename $R1 _1_f.fastq.gz)
    
    # Define output prefix
    output_prefix="$OUTPUT_DIR/${sample_name}_f_"
    
    echo "Processing sample: $sample_name"
    
    # Run STAR aligner
    STAR --runThreadN $THREADS \
          --genomeDir $GENOME_DIR \
          --readFilesIn $R1 $R2 \
          --readFilesCommand zcat \
          --outFileNamePrefix $output_prefix \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMstrandField intronMotif \
          --outSAMattributes All \
          --outFilterType BySJout
    echo "Finished processing sample: $sample_name"
done
echo "All samples processed."
