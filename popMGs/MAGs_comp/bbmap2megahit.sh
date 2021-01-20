#!/bin/bash

# REFDIR=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-evolution/MG-timeseries
# MGDIR=/Volumes/Bennett_BACKUP/Research/curtobacterium/raw_data/metagenomes/2017_evolutionMGs/novaseq_highres/libraries

REFDIR=/dfs5/bio/abchase/genomes/evolution/timeseries
MGDIR=$REFDIR/raw_reads

THREAD=16

cd $REFDIR

manifest=$REFDIR/MGmanifest.txt

echo -e "sampleID\tfilteredreads\tmappedreads" > $REFDIR/readdata.txt

# module load BBMap/37.50 
# bbmap.sh ref=/dfs5/bio/abchase/genomes/ref_genome/MMLR14002_pacbio-corrected.fasta path=$REFDIR

while read sampleID wellID i5 i7 primer1 primer2 runID
do

	read1=$(ls $MGDIR/*-${runID}-*${primer2}-${primer1}*.txt.gz | grep "READ1")
	read2=$(ls $MGDIR/*-${runID}-*${primer2}-${primer1}*.txt.gz | grep "READ2")

	echo "#!/bin/bash
#$ -N ${sampleID}.mega
#$ -m a
#$ -q bio
#$ -pe openmp $THREAD

module load samtools/1.3
module load BBMap/37.50 
module load megahit/1.1.1
module load blast/2.2.30
module load bowtie2/2.2.7

REF=${sampleID}
FFILE=$read1
RFILE=$read2

REFGENOME=/dfs5/bio/abchase/genomes/ref_genome/MMLR14002_pacbio-corrected.fasta
REFDIR=$REFDIR
OUTDIR=\$REFDIR/\$REF
BBMAPDIR=/data/users/abchase

THREAD=$THREAD

cd \$REFDIR

#############################################
#### QC filter the reads first
#### also will remove Illumina adapters 
#############################################
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


#############################################
#### map the reads to the reference genome 
#### will use a soft 97% mapping for now
#############################################
bbmap.sh in1=\$REF.filter.clean.R1.fq.gz in2=\$REF.filter.clean.R2.fq.gz \\
minid=0.97 \\
ref=\$REFGENOME path=\$REFDIR \\
out=\${REF}.mapped.sam \\
outm=\${REF}.mapped.fq scafstats=\$REF.stats.txt 

sortbyname.sh in=\${REF}.mapped.sam out=\${REF}.mapped.sorted.sam

\$BBMAPDIR/shrinksam-master/shrinksam \\
-i \${REF}.mapped.sorted.sam -k \${REF}.mapped.sortedSH.sam

rm -f \${REF}.mapped.sam
rm -f \${REF}.mapped.sorted.sam

mappedreads=\$(cat \$REF.stats.txt | grep \"ctg7180000000001_pilon\" | cut -f6)
rm -f \$REF.stats.txt

echo -e \"\${REF}\\t\${inputreads}\\t\${mappedreads}\" >> \$REFDIR/readdata.txt

reformat.sh in=\$REF.mapped.fq \\
out1=\$REF.R1.mapped.fq.gz out2=\$REF.R2.mapped.fq.gz ow=t

#### \$REF.R12.mapped.fq.gz files will be used for the SNP calling later so keep them

rm -f \$REF.filter.clean.R1.fq.gz
rm -f \$REF.filter.clean.R2.fq.gz

#############################################
#### assemble reads from the mapping
#############################################
rm -rf \$OUTDIR

megahit \\
-1 \$REF.R1.mapped.fq.gz \\
-2 \$REF.R2.mapped.fq.gz \\
-t \$THREAD \\
--min-count 3 \\
--k-list 31,41,51,61,71,81,91,95,101,105,111 \\
--min-contig-len 1000 \\
--memory 0.95 \\
--out-dir \$OUTDIR \\
--continue

#############################################
#### generate BLOB plot of assembly
#############################################

BIN=/data/users/abchase/bin
BLASTDB=/data/apps/commondata/blastdb
TAXDUMP=/dfs5/bio/abchase/taxdump_ncbi

cd \$OUTDIR

ASSEMBLY=final.contigs.fa

blastn -task megablast -query \$ASSEMBLY -db \$BLASTDB/nt -evalue 1e-5 -max_target_seqs 1 \\
-num_threads \$THREAD -outfmt '6 qseqid staxids' -out \$REF.nt.1e-5.megablast

bowtie2-build \$ASSEMBLY \$ASSEMBLY

bowtie2 -x \$ASSEMBLY --very-fast-local -k 1 -t -p \$THREAD --reorder --mm \\
-U \$REFDIR/\${REF}.mapped.fq | samtools view -S -b -T \$ASSEMBLY - > \$REF.bowtie2.bam

\$BIN/gc_cov_annotate.pl \\
--blasttaxid \$REF.nt.1e-5.megablast \\
--assembly \$ASSEMBLY --bam *.bam --out \$REF.bowtie.blobplot.txt \\
--taxdump \$TAXDUMP --taxlist genus family phylum

cp \$ASSEMBLY \$REFDIR/MAGs/\$REF.megahit.fasta
cp \$REF.bowtie.blobplot.txt \$REFDIR/MAGs/\$REF.bowtie.blobplot.txt

cd \$REFDIR
rm -rf \$OUTDIR/

rm -rf \$REFDIR/\${REF}.mapped.fq
rm -rf \$REFDIR/\${REF}.mapped.norm.fq


	" > $REFDIR/${sampleID}.sh

	qsub $REFDIR/${sampleID}.sh

done < $manifest
