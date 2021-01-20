#!/bin/bash

REFBASE=/dfs5/bio/abchase/genomes/evolution/timeseries
OUTDIR=$REFBASE/breseq

rm -rf $OUTDIR
mkdir -p $OUTDIR

THREAD=16

cd $REFBASE

# reads were already filtered and mapped to the reference strain from de novo assembly step
REFGENOME=/dfs5/bio/abchase/genomes/ref_genome/Curtobacterium_MMLR14002.gbk


for f in *.R1.mapped.fq.gz 
do
	output=$(echo $f | cut -f1 -d'.')
	read2=${output}.R2.mapped.fq.gz 

	gread1=${f%.gz}
	gread2=${read2%.gz}

	echo "#!/bin/bash
#$ -N ${output}.BRESEQ
#$ -q bio
#$ -j y
#$ -pe openmp $THREAD

module load samtools/1.3
module load bowtie2/2.2.7
module load R/3.4.1
module load breseq/0.26.1
module load BBMap/37.50

REF=${output}
READBASE=$REFBASE
OUTDIR=\$READBASE/breseq/\$REF

FFILE=\$READBASE/$f
RFILE=\$READBASE/$read2

REFGENOME=$REFGENOME

cd \$READBASE/breseq

bbnorm.sh \\
in=\$FFILE in2=\$RFILE \\
target=300 mindepth=5 threads=${THREAD} \\
out=\$READBASE/breseq/$gread1 out2=\$READBASE/breseq/$gread2

breseq -j $THREAD -p -r \$REFGENOME \\
--require-match-fraction 0.97 \\
--brief-html-output --no-junction-prediction \\
-o \$OUTDIR \\
\$READBASE/breseq/$gread1 \$READBASE/breseq/$gread2

rm -f \$READBASE/breseq/$gread1
rm -f \$READBASE/breseq/$gread2


cp \$OUTDIR/output/output.gd \$READBASE/breseq/\${REF}.gd
cp -R \$OUTDIR/output/ \$READBASE/breseq/\${REF}_output/
rm -rf \$OUTDIR/

cd \$READBASE/breseq
tar -zcvf \${REF}.breseq.tar.gz \${REF}_output/
rm -rf \${REF}_output/


	" > $OUTDIR/${output}.sh

	qsub $OUTDIR/${output}.sh

sleep 2s

done

# ##### once complete
# gdtools COMPARE -o compare.html -r $REFGENOME `ls *.gd`
# gdtools COUNT -o count.csv -r $REFGENOME `ls *.gd`
# gdtools COMPARE -o compare.gd -r $REFGENOME `ls *.gd`

# for sample in desert grass pine subalpine wood
# do

# 	gdtools COMPARE -o compare_${sample}.html -r $REFGENOME `ls ${sample}*.gd`
# 	gdtools COUNT -o count_${sample}.csv -r $REFGENOME `ls ${sample}*.gd`
# 	gdtools COMPARE -o compare_${sample}.gd -r $REFGENOME `ls ${sample}*.gd`

# done

