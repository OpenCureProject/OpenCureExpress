#!/bin/bash

# Set up environment variables
BASE_DIR=~/OpenCure-Fold
TRIMMED_DIR=${BASE_DIR}/trimmed_reads
TRIMMED_NORMAL_DIR=${TRIMMED_DIR}/normal
TRIMMED_TUMOR_DIR=${TRIMMED_DIR}/tumor
ALIGNMENTS_DIR=${BASE_DIR}/alignments
ALIGNMENTS_NORMAL_DIR=${ALIGNMENTS_DIR}/normal
ALIGNMENTS_TUMOR_DIR=${ALIGNMENTS_DIR}/tumor
VARIANTS_DIR=${BASE_DIR}/variants
QC_REPORTS_DIR=${BASE_DIR}/qc_reports
QC_NORMAL_DIR=${QC_REPORTS_DIR}/normal
QC_TUMOR_DIR=${QC_REPORTS_DIR}/tumor
REF_DIR=${BASE_DIR}/reference
REF_FA=${REF_DIR}/hg38.fa
INDEX_PREFIX=${REF_DIR}/hg38
ADAPTER_DIR=${BASE_DIR}/adapters
ADAPTER_PATH=${ADAPTER_DIR}/TruSeq3-PE.fa
KNOWN_SITES_DIR=${BASE_DIR}/known_sites
DBSNP_VCF=${KNOWN_SITES_DIR}/dbsnp.vcf.gz
SNPEFF_DIR=${BASE_DIR}/snpEff
SNPEFF_CONFIG=${SNPEFF_DIR}/snpEff.config
SNPEFF_JAR=${SNPEFF_DIR}/snpEff.jar

# Create directories
mkdir -p ${ALIGNMENTS_NORMAL_DIR} ${ALIGNMENTS_TUMOR_DIR} ${VARIANTS_DIR}
mkdir -p ${QC_NORMAL_DIR} ${QC_TUMOR_DIR}
mkdir -p ${REF_DIR} ${ADAPTER_DIR} ${KNOWN_SITES_DIR}
mkdir -p ${SNPEFF_DIR}

# Define lanes
NORMAL_LANES=(1 2 3)
TUMOR_LANES=(1 2 3 4 5)

# Download adapter file if it doesn't exist
if [ ! -f ${ADAPTER_PATH} ]; then
  wget -O ${ADAPTER_PATH} https://raw.githubusercontent.com/usadellab/Trimmomatic/master/adapters/TruSeq3-PE.fa
fi

# Run FastQC on trimmed reads
fastqc -t 28 ${TRIMMED_NORMAL_DIR}/WGS_Norm_Lane*_R*_paired.fastq.gz -o ${QC_NORMAL_DIR}
fastqc -t 28 ${TRIMMED_TUMOR_DIR}/WGS_Tumor_Lane*_R*_paired.fastq.gz -o ${QC_TUMOR_DIR}

# Download and prepare reference genome
# Check if the reference genome file exists
if [ ! -f ${REF_FA} ]; then
  echo "Downloading reference genome..."
  wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O ${REF_FA}.gz
  gunzip ${REF_FA}.gz
else
  echo "Reference genome file ${REF_FA} already exists. Skipping download."
fi

# Index the reference genome with BWA-MEM2 if index files don't exist
if [ ! -f ${INDEX_PREFIX}.0123 ]; then
  echo "Indexing reference genome with BWA-MEM2..."
  bwa-mem2 index -p ${INDEX_PREFIX} ${REF_FA}
else
  echo "BWA-MEM2 index files already exist. Skipping indexing."
fi

# Create SAMtools index if it doesn't exist
if [ ! -f ${REF_FA}.fai ]; then
  echo "Creating SAMtools index..."
  samtools faidx ${REF_FA}
else
  echo "SAMtools index ${REF_FA}.fai already exists. Skipping."
fi

# Create sequence dictionary for GATK if it doesn't exist
if [ ! -f ${REF_DIR}/hg38.dict ]; then
  echo "Creating sequence dictionary for GATK..."
  gatk CreateSequenceDictionary -R ${REF_FA} -O ${REF_DIR}/hg38.dict
else
  echo "Sequence dictionary ${REF_DIR}/hg38.dict already exists. Skipping."
fi

# Align normal samples
for lane in "${NORMAL_LANES[@]}"; do
  echo "Aligning Normal Sample Lane ${lane}..."
  bwa-mem2 mem -t 28 ${INDEX_PREFIX} \
    ${TRIMMED_NORMAL_DIR}/WGS_Norm_Lane${lane}_R1_paired.fastq.gz \
    ${TRIMMED_NORMAL_DIR}/WGS_Norm_Lane${lane}_R2_paired.fastq.gz | \
    samtools view -@ 4 -bS - > ${ALIGNMENTS_NORMAL_DIR}/normal_lane${lane}.bam

  echo "Sorting and Indexing Normal Sample Lane ${lane}..."
  samtools sort -@ 28 ${ALIGNMENTS_NORMAL_DIR}/normal_lane${lane}.bam -o ${ALIGNMENTS_NORMAL_DIR}/normal_lane${lane}_sorted.bam
  samtools index ${ALIGNMENTS_NORMAL_DIR}/normal_lane${lane}_sorted.bam
done

echo "Merging Normal Sample BAM files..."
samtools merge -@ 28 ${ALIGNMENTS_NORMAL_DIR}/normal_merged.bam ${ALIGNMENTS_NORMAL_DIR}/normal_lane*_sorted.bam
samtools sort -@ 28 ${ALIGNMENTS_NORMAL_DIR}/normal_merged.bam -o ${ALIGNMENTS_NORMAL_DIR}/normal_sorted.bam
samtools index ${ALIGNMENTS_NORMAL_DIR}/normal_sorted.bam

# Align tumor samples
for lane in "${TUMOR_LANES[@]}"; do
  echo "Aligning Tumor Sample Lane ${lane}..."
  bwa-mem2 mem -t 28 ${INDEX_PREFIX} \
    ${TRIMMED_TUMOR_DIR}/WGS_Tumor_Lane${lane}_R1_paired.fastq.gz \
    ${TRIMMED_TUMOR_DIR}/WGS_Tumor_Lane${lane}_R2_paired.fastq.gz | \
    samtools view -@ 4 -bS - > ${ALIGNMENTS_TUMOR_DIR}/tumor_lane${lane}.bam

  echo "Sorting and Indexing Tumor Sample Lane ${lane}..."
  samtools sort -@ 28 ${ALIGNMENTS_TUMOR_DIR}/tumor_lane${lane}.bam -o ${ALIGNMENTS_TUMOR_DIR}/tumor_lane${lane}_sorted.bam
  samtools index ${ALIGNMENTS_TUMOR_DIR}/tumor_lane${lane}_sorted.bam
done

echo "Merging Tumor Sample BAM files..."
samtools merge -@ 28 ${ALIGNMENTS_TUMOR_DIR}/tumor_merged.bam ${ALIGNMENTS_TUMOR_DIR}/tumor_lane*_sorted.bam
samtools sort -@ 28 ${ALIGNMENTS_TUMOR_DIR}/tumor_merged.bam -o ${ALIGNMENTS_TUMOR_DIR}/tumor_sorted.bam
samtools index ${ALIGNMENTS_TUMOR_DIR}/tumor_sorted.bam

# Mark duplicates
export _JAVA_OPTIONS="-Xmx16g"

# Mark duplicates for normal sample if output doesn't exist
if [ ! -f ${ALIGNMENTS_NORMAL_DIR}/normal_dedup.bam ]; then
  echo "Marking duplicates in normal sample..."
  picard MarkDuplicates \
    I=${ALIGNMENTS_NORMAL_DIR}/normal_sorted.bam \
    O=${ALIGNMENTS_NORMAL_DIR}/normal_dedup.bam \
    M=${ALIGNMENTS_NORMAL_DIR}/normal_markduplicates_metrics.txt \
    CREATE_INDEX=true
else
  echo "Deduplicated BAM for normal sample already exists. Skipping."
fi

# Mark duplicates for tumor sample if output doesn't exist
if [ ! -f ${ALIGNMENTS_TUMOR_DIR}/tumor_dedup.bam ]; then
  echo "Marking duplicates in tumor sample..."
  picard MarkDuplicates \
    I=${ALIGNMENTS_TUMOR_DIR}/tumor_sorted.bam \
    O=${ALIGNMENTS_TUMOR_DIR}/tumor_dedup.bam \
    M=${ALIGNMENTS_TUMOR_DIR}/tumor_markduplicates_metrics.txt \
    CREATE_INDEX=true
else
  echo "Deduplicated BAM for tumor sample already exists. Skipping."
fi

# Download known sites if they don't exist
if [ ! -f ${DBSNP_VCF} ] || [ ! -f ${DBSNP_VCF}.tbi ]; then
  echo "Downloading dbSNP VCF..."
  wget ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz -O ${DBSNP_VCF}
  wget ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz.tbi -O ${DBSNP_VCF}.tbi
else
  echo "dbSNP VCF and index already exist. Skipping download."
fi

# Base Quality Score Recalibration (BQSR)
# BQSR for normal sample
if [ ! -f ${ALIGNMENTS_NORMAL_DIR}/normal_bqsr.bam ]; then
  echo "Performing BQSR on normal sample..."
  gatk --java-options "-Xmx16g" BaseRecalibrator \
    -I ${ALIGNMENTS_NORMAL_DIR}/normal_dedup.bam \
    -R ${REF_FA} \
    --known-sites ${DBSNP_VCF} \
    -O ${ALIGNMENTS_NORMAL_DIR}/normal_recal_data.table

  gatk --java-options "-Xmx16g" ApplyBQSR \
    -I ${ALIGNMENTS_NORMAL_DIR}/normal_dedup.bam \
    -R ${REF_FA} \
    --bqsr-recal-file ${ALIGNMENTS_NORMAL_DIR}/normal_recal_data.table \
    -O ${ALIGNMENTS_NORMAL_DIR}/normal_bqsr.bam
else
  echo "BQSR BAM for normal sample already exists. Skipping."
fi

# BQSR for tumor sample
if [ ! -f ${ALIGNMENTS_TUMOR_DIR}/tumor_bqsr.bam ]; then
  echo "Performing BQSR on tumor sample..."
  gatk --java-options "-Xmx16g" BaseRecalibrator \
    -I ${ALIGNMENTS_TUMOR_DIR}/tumor_dedup.bam \
    -R ${REF_FA} \
    --known-sites ${DBSNP_VCF} \
    -O ${ALIGNMENTS_TUMOR_DIR}/tumor_recal_data.table

  gatk --java-options "-Xmx16g" ApplyBQSR \
    -I ${ALIGNMENTS_TUMOR_DIR}/tumor_dedup.bam \
    -R ${REF_FA} \
    --bqsr-recal-file ${ALIGNMENTS_TUMOR_DIR}/tumor_recal_data.table \
    -O ${ALIGNMENTS_TUMOR_DIR}/tumor_bqsr.bam
else
  echo "BQSR BAM for tumor sample already exists. Skipping."
fi

# Get sample IDs
TUMOR_SAMPLE_ID=$(samtools view -H ${ALIGNMENTS_TUMOR_DIR}/tumor_bqsr.bam | grep '^@RG' | awk -F'\t' '{for(i=1;i<=NF;i++){if($i ~ /^SM:/){print substr($i,4)}}}' | uniq)
NORMAL_SAMPLE_ID=$(samtools view -H ${ALIGNMENTS_NORMAL_DIR}/normal_bqsr.bam | grep '^@RG' | awk -F'\t' '{for(i=1;i<=NF;i++){if($i ~ /^SM:/){print substr($i,4)}}}' | uniq)

# Variant calling with Mutect2
if [ ! -f ${VARIANTS_DIR}/somatic_raw.vcf.gz ]; then
  echo "Running Mutect2 for somatic variant calling..."
  gatk --java-options "-Xmx16g" Mutect2 \
    -R ${REF_FA} \
    -I ${ALIGNMENTS_TUMOR_DIR}/tumor_bqsr.bam \
    -I ${ALIGNMENTS_NORMAL_DIR}/normal_bqsr.bam \
    --tumor-sample ${TUMOR_SAMPLE_ID} \
    --normal-sample ${NORMAL_SAMPLE_ID} \
    --output ${VARIANTS_DIR}/somatic_raw.vcf.gz \
    --native-pair-hmm-threads 28
else
  echo "Somatic raw VCF already exists. Skipping Mutect2."
fi

# Filter Mutect2 calls
if [ ! -f ${VARIANTS_DIR}/somatic_filtered.vcf.gz ]; then
  echo "Filtering Mutect2 variant calls..."
  gatk --java-options "-Xmx16g" FilterMutectCalls \
    -V ${VARIANTS_DIR}/somatic_raw.vcf.gz \
    -O ${VARIANTS_DIR}/somatic_filtered.vcf.gz
else
  echo "Filtered somatic VCF already exists. Skipping FilterMutectCalls."
fi

# Annotate variants with SnpEff
if [ ! -f ${VARIANTS_DIR}/somatic_annotated.vcf ]; then
  echo "Annotating variants with SnpEff..."

  # Install SnpEff if not already installed
  if [ ! -f ${SNPEFF_JAR} ]; then
    echo "Downloading SnpEff..."
    wget -O ${SNPEFF_DIR}/snpEff_latest_core.zip https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip ${SNPEFF_DIR}/snpEff_latest_core.zip -d ${SNPEFF_DIR}
    rm ${SNPEFF_DIR}/snpEff_latest_core.zip
  fi

  # Download SnpEff database for GRCh38 if not already present
  if [ ! -d ${SNPEFF_DIR}/data/GRCh38.99 ]; then
    echo "Downloading SnpEff database for GRCh38.99..."
    java -jar ${SNPEFF_JAR} download -v GRCh38.99
  fi

  # Run SnpEff annotation
  java -Xmx16g -jar ${SNPEFF_JAR} GRCh38.99 \
    -v \
    -stats ${VARIANTS_DIR}/snpEff_summary.html \
    ${VARIANTS_DIR}/somatic_filtered.vcf.gz > ${VARIANTS_DIR}/somatic_annotated.vcf
else
  echo "Annotated VCF already exists. Skipping SnpEff."
fi

# Generate MultiQC reports --- Revisit this
echo "Generating MultiQC reports..."
multiqc ${QC_NORMAL_DIR} -o ${QC_NORMAL_DIR}
multiqc ${QC_TUMOR_DIR} -o ${QC_TUMOR_DIR}