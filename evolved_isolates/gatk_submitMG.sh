#!/bin/bash

REFBASE=/dfs5/bio/abchase/genomes/evolution_complete

cd $REFBASE

# reads were already filtered and mapped to the reference strain from de novo assembly step
# just need to reference those reads for this and apply a local alignment (i.e. BWA)

count=1

# index Reference - only need to do once!
REFGENOME=/dfs5/bio/abchase/genomes/ref_genome/MMLR14002_pacbio-corrected.fasta
# bwa index $REFGENOME
# samtools faidx $REFGENOME
# gatk CreateSequenceDictionary --REFERENCE=$REFGENOME --OUTPUT=${REFGENOME%.fasta}.dict

while read sampleID FFILE RFILE 
do
	
	echo "#!/bin/bash
#$ -N ${sampleID}
#$ -q mic,pub8i
#$ -pe openmp 8

module load BBMap/37.50
module load gatk/4.0.0.0
module load bwa/0.7.8
module load samtools/1.8-11  
module load R/3.4.1
module load bedtools/2.25.0  

REF=${sampleID}
READBASE=$REFBASE
BBMAPDIR=/data/users/abchase
FFILE=\$READBASE/genome_raw/$FFILE
RFILE=\$READBASE/genome_raw/$RFILE

OUTDIR=\$READBASE/gatk/$REF
THREAD=8

REFGENOME=$REFGENOME

cd \$READBASE/gatk

bbduk.sh qtrim=rl trimq=10 -Xmx20g threads=\$THREAD \\
minlen=25 ktrim=r k=25 ref=\$BBMAPDIR/nextera.fa.gz hdist=1 \\
in1=\$FFILE in2=\$RFILE \\
out1=\$REF.filter.R1.fq out2=\$REF.filter.R2.fq stats=\$REF.stats.txt

inputreads=\$(cat \$REF.stats.txt | grep \"Total\" | cut -f2)
rm -f \$REF.stats.txt

repair.sh in=\$REF.filter.R1.fq in2=\$REF.filter.R2.fq \\
out=\$REF.filter.clean.R1.fq.gz out2=\$REF.filter.clean.R2.fq.gz

rm -f \$REF.filter.R1.fq
rm -f \$REF.filter.R2.fq

######### Alignment – Map to Reference
## Need to provide the -M flag to BWA, this tells it to consider split reads as secondary, need this for GATK variant calling/Picard support. 
## Readgroup info is provided with the -R flag. This information is key for downstream GATK functionality. 
## GATK will not work without a read group tag.

bwa mem -M -t 8 -R '@RG\tID:sample_${count}\tLB:sample_${count}\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_${count}' \\
\$REFGENOME \$REF.filter.clean.R1.fq.gz \$REF.filter.clean.R2.fq.gz > \${REF}.aligned_reads.sam

# Sort SAM file by coordinate, convert to BAM
gatk SortSam --INPUT=\${REF}.aligned_reads.sam --OUTPUT=\${REF}.samsorted_reads.bam --SORT_ORDER=coordinate

rm \${REF}.aligned_reads.sam

# Collect Alignment & Insert Size Metrics
gatk CollectAlignmentSummaryMetrics -R=\$REFGENOME -I=\${REF}.samsorted_reads.bam -O=\${REF}.alignment_metrics.txt

gatk CollectInsertSizeMetrics --INPUT=\${REF}.samsorted_reads.bam \\
--OUTPUT=\${REF}.insert_metrics.txt --Histogram_FILE=\${REF}.insert_size_histogram.pdf

samtools depth -a \${REF}.samsorted_reads.bam > \${REF}.depth_out.txt

# Mark Duplicates
gatk MarkDuplicates --INPUT=\${REF}.samsorted_reads.bam --OUTPUT=\${REF}.dedup_reads.bam --METRICS_FILE=\${REF}.metrics.txt

rm \${REF}.samsorted_reads.bam

# Build BAM Index
gatk BuildBamIndex --INPUT=\${REF}.dedup_reads.bam

# Call Variants
gatk HaplotypeCaller --reference \$REFGENOME --input \${REF}.dedup_reads.bam --sample-ploidy 1 \\
--standard-min-confidence-threshold-for-calling 30 \\
--native-pair-hmm-threads 8 \\
--output \${REF}.raw.vcf

# Extract SNPs & Indels
gatk SelectVariants -R \$REFGENOME --variant \${REF}.raw.vcf \\
--select-type-to-include SNP --output \${REF}.raw_snps.vcf

gatk SelectVariants -R \$REFGENOME --variant \${REF}.raw.vcf \\
--select-type-to-include INDEL --output \${REF}.raw_indels.vcf

# Filter SNPs
gatk VariantFiltration -R \$REFGENOME --variant \${REF}.raw_snps.vcf \\
--filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \\
--filter-name \"basic_snp_filter\" --output \${REF}.filtered_snps.vcf

# filter indels
gatk VariantFiltration -R \$REFGENOME --variant \${REF}.raw_indels.vcf \\
--filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \\
--filter-name \"basic_indel_filter\" --output \${REF}.filtered_indels.vcf

rm \${REF}.raw_snps.vcf
rm \${REF}.raw_indels.vcf
rm \${REF}.raw.vcf

# Base Quality Score Recalibration (BQSR) #1
# BQSR is performed twice. The second pass is optional, but is required to produce a recalibration report.
gatk BaseRecalibrator -R \$REFGENOME --input \${REF}.dedup_reads.bam --known-sites \${REF}.filtered_snps.vcf \\
--known-sites \${REF}.filtered_indels.vcf --output \${REF}.recal_data.table

gatk ApplyBQSR -I \${REF}.dedup_reads.bam --bqsr-recal-file \${REF}.recal_data.table -O \${REF}.dedup_reads.recalibrated.bam

rm -f \${REF}.dedup_reads.bam
rm -f \${REF}.dedup_reads.bai

# Base Quality Score Recalibration (BQSR) #2
gatk BaseRecalibrator -R \$REFGENOME --input \${REF}.dedup_reads.recalibrated.bam --known-sites \${REF}.filtered_snps.vcf \\
--known-sites \${REF}.filtered_indels.vcf --output \${REF}.post_recal_data.table

# Analyze Covariates
gatk AnalyzeCovariates -before \${REF}.recal_data.table -after \${REF}.post_recal_data.table \\
-plots \${REF}.recalibration_plots.pdf

# Apply BQSR
gatk ApplyBQSR -I \${REF}.dedup_reads.recalibrated.bam --bqsr-recal-file \${REF}.post_recal_data.table -O \${REF}.dedup_reads.recalibrated2.bam

rm -f \${REF}.dedup_reads.recalibrated.bam
rm -f \${REF}.dedup_reads.recalibrated.bai

rm -f \${REF}.filtered_indels.vcf
rm -f \${REF}.filtered_snps.vcf
rm -f \${REF}.filtered_indels.vcf.idx
rm -f \${REF}.filtered_snps.vcf.idx
rm -f \${REF}.raw_snps.vcf.idx
rm -f \${REF}.raw_indels.vcf.idx

# Call Variants - Second round of variant calling performed on recalibrated bam
gatk HaplotypeCaller --reference \$REFGENOME --input \${REF}.dedup_reads.recalibrated2.bam \\
--native-pair-hmm-threads 8 -ploidy 1 \\
--output \${REF}.raw_variants_recal.vcf

# Extract SNPs & Indels
gatk SelectVariants -R \$REFGENOME --variant \${REF}.raw_variants_recal.vcf \\
--select-type-to-include SNP --output \${REF}.raw_snps_recal.vcf

gatk SelectVariants -R \$REFGENOME --variant \${REF}.raw_variants_recal.vcf \\
--select-type-to-include INDEL --output \${REF}.raw_indels_recal.vcf

# filter
gatk VariantFiltration -R \$REFGENOME --variant \${REF}.raw_snps_recal.vcf \\
--filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \\
--filter-name \"basic_snp_filter\" --output \${REF}.filtered_snps_final.vcf

gatk VariantFiltration -R \$REFGENOME --variant \${REF}.raw_indels_recal.vcf \\
--filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \\
--filter-name \"basic_indel_filter\" --output \${REF}.filtered_indels_recal.vcf

# merge results for each
gatk MergeVcfs --INPUT \${REF}.filtered_snps_final.vcf \\
--INPUT \${REF}.filtered_indels_recal.vcf \\
--OUTPUT \${REF}.filtered_recal_final.vcf

# Compute Coverage Statistics
bedtools genomecov -d -ibam \${REF}.dedup_reads.recalibrated2.bam > \${REF}.bedgraph


rm -f \${REF}.raw_variants_recal.vcf
rm -f \${REF}.raw_snps_recal.vcf
rm -f \${REF}.raw_indels_recal.vcf

rm -f \${REF}.dedup_reads.recalibrated2.bam
rm -f \${REF}.dedup_reads.recalibrated2.bai

rm -f \${REF}.filtered_indels.vcf.idx
rm -f \${REF}.filtered_snps.vcf.idx
rm -f \${REF}.raw_snps.vcf.idx
rm -f \${REF}.raw_indels.vcf.idx
rm -f \${REF}.raw_snps_recal.vcf.idx
rm -f \${REF}.raw_indels_recal.vcf.idx
rm -f \${REF}.raw_variants_recal.vcf.idx
rm -f \${REF}.filtered_indels_recal.vcf.idx
rm -f \${REF}.filtered_snps_final.vcf.idx
rm -f \${REF}.metrics.txt
rm -f \${REF}.post_recal_data.table
rm -f \${REF}.raw.vcf.idx
rm -f \${REF}.recal_data.table
rm -f \${REF}.insert_metrics.txt
rm -f \${REF}.depth_out.txt
rm -f \${REF}.alignment_metrics.txt

	" > $REFBASE/gatk/${sampleID}.gatk.sh

	count=`expr $count + 1`

done < good-reads.txt 

### when done 
### bedtools unionbedg -i curto*.bedgraph -header > total.coverage.bedgraph
# module load bcftools/1.9
# for f in *.filtered_recal_final.vcf
# do 
# 	bcftools view -Oz -o ${f%.filtered_recal_final.vcf}.compressed.vcf.gz $f
# 	htsfile ${f%.filtered_recal_final.vcf}.compressed.vcf.gz
# 	bcftools index ${f%.filtered_recal_final.vcf}.compressed.vcf.gz
# done
# bcftools merge *.compressed.vcf.gz -Oz -o merged.vcf.gz



