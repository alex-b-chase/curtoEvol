#!/bin/bash

# take a metagenomic file, parse it for marker genes in the database

BASEDIR=/bio/abchase/MG/elevation_transplant
MGDIR1=$BASEDIR/initialtime
MGDIR2=$BASEDIR/finaltime

OUTDIR=$BASEDIR/comm_markers

# using a new method to reduce computation time, substituting BLAT for BLAST
# with this, don't need to partition the reads into smaller chunks
# BLAT should be able to parse this in a few hours for each library from benchmarking on local machine

cd $MGDIR1

for f in *.filter.total.faa
do

	output=${f%.filter.total.faa}

	echo "#!/bin/bash
#$ -N ${output}.blat
#$ -m a
#$ -q pub8i,bio,mic

module load enthought_python/7.3.2
module load blat/36
module load BBMap/37.50

MGDIR=$MGDIR1
BLASTDB=/bio/abchase/refDB/total_markers.faa
OUTDIR=$OUTDIR

cd \${MGDIR}

blat -prot -fastMap -minIdentity=20 -out=blast8 \\
\$BLASTDB $f temp.20.${output}.txt

# subset the giant output file for only relevant information (i.e. query sequence and protein match)
cut -f1-2 temp.20.${output}.txt | \\
awk 'BEGIN{FS=\"\t\"; OFS=\"\t\"} {gsub(/^[^_]*_/, \"\", \$2); print}' | \\
sort -u > \$OUTDIR/${output}.blat.total.txt

rm -f temp.20.${output}.txt

# now subet the marker gene reads from the MG library
cut -f1 \$OUTDIR/${output}.blat.total.txt | sort -u | sed 's/[ ]*$//' > \$OUTDIR/${output}.blat.temp.txt

# filter each MG by marker gene using BBMap software (WAY faster than my own scripts)
filterbyname.sh \\
in=$f \\
out=\$OUTDIR/${output}.blat.markers.faa \\
names=\$OUTDIR/${output}.blat.temp.txt ow=t include=t 2>/dev/null

rm -f \$OUTDIR/${output}.blat.temp.txt



	" > $output.sh

	qsub $output.sh

	sleep 2s
done

