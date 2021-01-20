#!/bin/bash

REFBASE=/dfs3/bio/abchase/genomes/evolution/timeseries

cd $REFBASE

# reads were already filtered and mapped to the reference strain from de novo assembly step
# just need to reference those reads for this and apply a local alignment (i.e. BWA)

count=1

# index Reference - only need to do once!
REFGENOME=/dfs3/bio/abchase/genomes/ref_genome/MMLR14002_pacbio-corrected.fasta
# bwa index $REFGENOME
# samtools faidx $REFGENOME
# gatk CreateSequenceDictionary --REFERENCE=$REFGENOME --OUTPUT=${REFGENOME%.fasta}.dict

for f in *.R1.mapped.fq.gz 
do
	output=$(echo $f | cut -f1 -d'.')
	read2=${output}.R2.mapped.fq.gz 

	for hapcount in {1..8} 
	do

		OUTDIR=$REFBASE/gatkH${hapcount}
		mkdir -p $OUTDIR

		echo "#!/bin/bash
#$ -N ${output}.gat${hapcount}
#$ -m a
#$ -q bio,pub64,pub8i
#$ -j y
#$ -pe openmp 8

# Identifying genomic variants, such as single nucleotide polymorphisms (SNPs) and DNA insertions and deletions (indels)

module load gatk/4.0.0.0
module load bwa/0.7.8
module load samtools/1.8-11  
module load R/3.4.1
module load bedtools/2.25.0  

REFBASE=$REFBASE
OUTDIR=\$REFBASE/gatkH${hapcount}
REFGENOME=$REFGENOME

output=${output}.H${hapcount}

mkdir -p \$OUTDIR

cd \$OUTDIR

# Alignment – Map to Reference
# Need to provide the -M flag to BWA, this tells it to consider split reads as secondary, need this for GATK variant calling/Picard support. 
# Readgroup info is provided with the -R flag. This information is key for downstream GATK functionality. 
# GATK will not work without a read group tag.

bwa mem -M -t 8 -R '@RG\tID:sample_${count}\tLB:sample_${count}\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_${count}' \\
\$REFGENOME \$REFBASE/${f} \$REFBASE/${read2} > \${output}.aligned_reads.sam

# Sort SAM file by coordinate, convert to BAM
gatk SortSam --INPUT=\${output}.aligned_reads.sam --OUTPUT=\${output}.samsorted_reads.bam --SORT_ORDER=coordinate

rm \${output}.aligned_reads.sam

# Collect Alignment & Insert Size Metrics
gatk CollectAlignmentSummaryMetrics -R=\$REFGENOME -I=\${output}.samsorted_reads.bam -O=\${output}.alignment_metrics.txt

gatk CollectInsertSizeMetrics --INPUT=\${output}.samsorted_reads.bam \\
--OUTPUT=\${output}.insert_metrics.txt --Histogram_FILE=\${output}.insert_size_histogram.pdf

samtools depth -a \${output}.samsorted_reads.bam > \${output}.depth_out.txt

# Mark Duplicates
gatk MarkDuplicates --INPUT=\${output}.samsorted_reads.bam --OUTPUT=\${output}.dedup_reads.bam --METRICS_FILE=\${output}.metrics.txt

rm \${output}.samsorted_reads.bam

# Build BAM Index
gatk BuildBamIndex --INPUT=\${output}.dedup_reads.bam

# Call Variants
gatk HaplotypeCaller --reference \$REFGENOME --input \${output}.dedup_reads.bam --sample-ploidy ${hapcount} \\
--standard-min-confidence-threshold-for-calling 30 \\
--native-pair-hmm-threads 8 \\
--output \${output}.raw.vcf

# Extract SNPs & Indels
gatk SelectVariants -R \$REFGENOME --variant \${output}.raw.vcf \\
--select-type-to-include SNP --output \${output}.raw_snps.vcf

gatk SelectVariants -R \$REFGENOME --variant \${output}.raw.vcf \\
--select-type-to-include INDEL --output \${output}.raw_indels.vcf

# Filter SNPs
gatk VariantFiltration -R \$REFGENOME --variant \${output}.raw_snps.vcf \\
--filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \\
--filter-name \"basic_snp_filter\" --output \${output}.filtered_snps.vcf

# filter indels
gatk VariantFiltration -R \$REFGENOME --variant \${output}.raw_indels.vcf \\
--filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \\
--filter-name \"basic_indel_filter\" --output \${output}.filtered_indels.vcf

rm \${output}.raw_snps.vcf
rm \${output}.raw_indels.vcf
rm \${output}.raw.vcf

# Base Quality Score Recalibration (BQSR) #1
# BQSR is performed twice. The second pass is optional, but is required to produce a recalibration report.
gatk BaseRecalibrator -R \$REFGENOME --input \${output}.dedup_reads.bam --known-sites \${output}.filtered_snps.vcf \\
--known-sites \${output}.filtered_indels.vcf --output \${output}.recal_data.table

gatk ApplyBQSR -I \${output}.dedup_reads.bam --bqsr-recal-file \${output}.recal_data.table -O \${output}.dedup_reads.recalibrated.bam

rm -f \${output}.dedup_reads.bam
rm -f \${output}.dedup_reads.bai

# Base Quality Score Recalibration (BQSR) #2
gatk BaseRecalibrator -R \$REFGENOME --input \${output}.dedup_reads.recalibrated.bam --known-sites \${output}.filtered_snps.vcf \\
--known-sites \${output}.filtered_indels.vcf --output \${output}.post_recal_data.table

# Analyze Covariates
gatk AnalyzeCovariates -before \${output}.recal_data.table -after \${output}.post_recal_data.table \\
-plots \${output}.recalibration_plots.pdf

# Apply BQSR
gatk ApplyBQSR -I \${output}.dedup_reads.recalibrated.bam --bqsr-recal-file \${output}.post_recal_data.table -O \${output}.dedup_reads.recalibrated2.bam

rm -f \${output}.dedup_reads.recalibrated.bam
rm -f \${output}.dedup_reads.recalibrated.bai

rm -f \${output}.filtered_indels.vcf
rm -f \${output}.filtered_snps.vcf
rm -f \${output}.filtered_indels.vcf.idx
rm -f \${output}.filtered_snps.vcf.idx
rm -f \${output}.raw_snps.vcf.idx
rm -f \${output}.raw_indels.vcf.idx

# Call Variants - Second round of variant calling performed on recalibrated bam
gatk HaplotypeCaller --reference \$REFGENOME --input \${output}.dedup_reads.recalibrated2.bam \\
--native-pair-hmm-threads 8 -ploidy ${hapcount} \\
--output \${output}.raw_variants_recal.vcf

# Extract SNPs & Indels
gatk SelectVariants -R \$REFGENOME --variant \${output}.raw_variants_recal.vcf \\
--select-type-to-include SNP --output \${output}.raw_snps_recal.vcf

gatk SelectVariants -R \$REFGENOME --variant \${output}.raw_variants_recal.vcf \\
--select-type-to-include INDEL --output \${output}.raw_indels_recal.vcf

# filter
gatk VariantFiltration -R \$REFGENOME --variant \${output}.raw_snps_recal.vcf \\
--filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \\
--filter-name \"basic_snp_filter\" --output \${output}.filtered_snps_final.vcf

gatk VariantFiltration -R \$REFGENOME --variant \${output}.raw_indels_recal.vcf \\
--filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \\
--filter-name \"basic_indel_filter\" --output \${output}.filtered_indels_recal.vcf

# merge results for each
gatk MergeVcfs --INPUT \${output}.filtered_snps_final.vcf \\
--INPUT \${output}.filtered_indels_recal.vcf \\
--OUTPUT \${output}.filtered_recal_final.vcf

# Compute Coverage Statistics
bedtools genomecov -bga -ibam \${output}.dedup_reads.recalibrated2.bam > \${output}.bedgraph


rm -f \${output}.raw_variants_recal.vcf
rm -f \${output}.raw_snps_recal.vcf
rm -f \${output}.raw_indels_recal.vcf

rm -f \${output}.dedup_reads.recalibrated2.bam
rm -f \${output}.dedup_reads.recalibrated2.bai

rm -f \${output}.filtered_indels.vcf.idx
rm -f \${output}.filtered_snps.vcf.idx
rm -f \${output}.raw_snps.vcf.idx
rm -f \${output}.raw_indels.vcf.idx
rm -f \${output}.raw_snps_recal.vcf.idx
rm -f \${output}.raw_indels_recal.vcf.idx
rm -f \${output}.raw_variants_recal.vcf.idx
rm -f \${output}.filtered_indels_recal.vcf.idx
rm -f \${output}.filtered_snps_final.vcf.idx
rm -f \${output}.metrics.txt
rm -f \${output}.post_recal_data.table
rm -f \${output}.raw.vcf.idx
rm -f \${output}.recal_data.table
rm -f \${output}.insert_metrics.txt
rm -f \${output}.depth_out.txt
rm -f \${output}.alignment_metrics.txt

		" > $OUTDIR/${output}.gatkH${hapcount}.sh

		qsub $OUTDIR/${output}.gatkH${hapcount}.sh

	done

	count=`expr $count + 1`

done