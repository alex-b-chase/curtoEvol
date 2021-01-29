#!/bin/bash

### THIS WILL APPLY A SECONDARY FILTER FOR THE MG READS FOR CURTO CORE GENES
### THE INPUT WILL BE THE FILTERED BLASTp READS *.blast.faa

# so each gene will most likely need an unique e-value
# we do not want to lose too much information so apply the loosest value to pass to PPLACER

# This will be OK since we applied a strict BLASTp filter (e-20)

REFDIR=/bio/abchase/MG/elevation_transplant/curto_coregenes
HMMPROF=/bio/abchase/refDB/coregenes/hmmprofiles
REFDB=/bio/abchase/refDB/coregenes/pplacer.refpkgs

OUTDIR=$REFDIR/pplacer

THREAD=4

cd $OUTDIR

# remove old files
rm -f *.e*
rm -f *.o*
rm -f *.jplace
rm -f *.csv
rm -f *.xml
rm -f *.log

cd $REFDIR/blast_filter

# remove old files from HMMer scan earlier (not needed right now)

echo -e "coregene\te-valueHMMer" > $OUTDIR/hmmerfiltered.txt

for blastout in *.blast.faa 
do 

	protein=${blastout%.blast.faa}

	echo "#!/bin/bash
#$ -N ${protein}_pplacer
#$ -q bio,pub8i,mic,abio
#$ -ckpt restart
#$ -j y
#$ -pe openmp ${THREAD}
# RESTART_EMAIL_SUMMARY

module load enthought_python/7.3.2
module load ea-utils/r811
module load phylosift/1.0.1
module load clcbio/8.5.1
module load BBMap/37.50
module load hmmer/3.1b2

REFDIR=${REFDIR}
HMMPROF=${HMMPROF}
REFDB=${REFDB}

OUTDIR=\$REFDIR/pplacer 

protein=${protein}

for minID in {40,45,50,55,60}
do

	cd \$REFDIR/blast_filter

	hmmsearch --tblout \${protein}.hmm.txt -E 1e-\${minID} \\
	--cpu ${THREAD} \$HMMPROF/\${protein}.hmm \${protein}.blast.faa > \${protein}.log.txt

	cat \${protein}.hmm.txt | cut -f1 -d' ' | sort | uniq > \${protein}.temp.txt

	# subset new HMMer filtered reads
	filterbyname.sh \\
	in=\${protein}.blast.faa \\
	out=total_\${protein}.\${minID}.hmm.faa \\
	names=\${protein}.temp.txt ow=t include=t 2>/dev/null

	rm -f \${protein}.temp.txt
	rm -f \${protein}.hmm.txt
	rm -f \${protein}.log.txt

	cat total_\${protein}.\${minID}.hmm.faa | tr -d '*' > \$OUTDIR/\${protein}.\${minID}.temp.hmm.faa

	cd \$OUTDIR

	rm -f \$OUTDIR/\${protein}.\${minID}.fa

	clustalo-1.2.0 --profile1 \$REFDB/\${protein}.refpkg/\${protein}.total.aln \\
	-i \${protein}.\${minID}.temp.hmm.faa -o \$OUTDIR/\${protein}.\${minID}.fa

	rm \${protein}.\${minID}.temp.hmm.faa

	### now can test input for pplacer 
	pplacer --pretend \\
	-c \$REFDB/\${protein}.refpkg \\
	\${protein}.\${minID}.fa > \${protein}.\${minID}.pplacer.log

	### check if pplacer worked
	if grep -Fxq \"everything looks OK.\" \${protein}.\${minID}.pplacer.log
	then
		echo -e \"\${protein}\\t\${minID}\" >> \$OUTDIR/hmmerfiltered.txt
		break

	else 
		rm -f \$REFDIR/blast_filter/total_\${protein}.\${minID}.hmm.faa
		rm -f \$OUTDIR/\${protein}.\${minID}.fa
		continue

	fi


done

cd \$OUTDIR

rm -f \${protein}.*.pplacer.log

### now can run pplacer 
pplacer -c \$REFDB/\${protein}.refpkg \\
\${protein}.\${minID}.fa \\
-p --keep-at-most 20 

guppy to_csv --point-mass --pp \${protein}.\${minID}.jplace > \${protein}.\${minID}.csv
guppy fat --node-numbers --point-mass --pp \${protein}.*.jplace


	" > $OUTDIR/${protein}.pplacer.sh

	
	qsub $OUTDIR/${protein}.pplacer.sh
	sleep 5


done




