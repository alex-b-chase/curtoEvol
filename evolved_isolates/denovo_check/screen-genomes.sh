#!/bin/bash

BASE=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-evolution/isolate_genomics_final/spades

REFDB=/Volumes/Bennett_BACKUP/Research/reference_db
ESSEN=$REFDB/essen_singlecopy

rm -rf $BASE/curated_assembly
mkdir -p $BASE/curated_assembly

rm -f $BASE/blobplot/*.fix.*
rm -f $BASE/blobplot/*.png


cd $BASE/genomes


# curate the genomes by filtering and renaming to finalize genomes
# first need to rename the contig IDs since prokka doesn't like long names <20 characters
for fas in *.scaffolds.fasta
do
	genome=${fas%.scaffolds.fasta}

	# filter out contigs that had <30 depth of coverage in the assembly
	# we are going to use the BLOBPLOT output to do so - column4
	cat $BASE/blobplot/$genome.bowtie.blobplot.txt | awk '(NR>1) && ($4 > 30 )' | cut -f1 > $genome'_temp.txt'

	search-fasta.py -i $fas -m $genome'_temp.txt' -o $genome.temp.fasta

	# remove contigs with less than 500bp and less than 55% GC%
	removesmalls.py -i $genome.temp.fasta -l 500 -o $genome.fix.fasta -g 55 > /dev/null 2>&1

	print-fasta-id.py $genome.fix.fasta

	# create new BLOBPLOT with fixed contigs
	cd $BASE/blobplot/

	# need to keep the header information from the file
	cat $genome.bowtie.blobplot.txt | head -n1 > $genome.fix.bowtie.blobplot.txt

	while read line
	do
		cat $genome.bowtie.blobplot.txt | grep $line >> $genome.fix.bowtie.blobplot.txt
	done < $BASE/genomes/$genome.fix_ids.txt

	$BASE/makeblobplot.R $genome.bowtie.blobplot.txt 0.01 taxlevel_genus 1 > /dev/null 2>&1
	$BASE/makeblobplot.R $genome.fix.bowtie.blobplot.txt 0.01 taxlevel_genus 1 > /dev/null 2>&1

	# now simply rename the files to match everything and be good to go for Prodigal
	cd $BASE/genomes

	cat $genome.fix_ids.txt | sed -e "s/^/$genome./" | cut -f1-2 -d '_' | sed 's/-//g' | sed 's/_//g' > $genome'_temp.txt'
	paste $genome.fix_ids.txt $genome'_temp.txt' > $genome'_newnames.txt'

	# now simply rename the files with homemade python script
	fasta-rename.py $genome.fix.fasta $genome'_newnames.txt' $genome.final.fasta
	

	# calculate some genomic stats
	n50_calc.py $fas > /dev/null 2>&1
	n50_calc.py $genome.final.fasta > /dev/null 2>&1

	mv $genome.final.fasta $BASE/curated_assembly

	rm $genome'_temp.txt'
	rm $genome.fix_ids.txt
	rm $genome.temp.fasta
	rm $genome'_newnames.txt'
	rm $genome.fix_rejected.fasta
	rm $genome.fix.fasta

	echo "done with ${genome}"


done

cd $BASE/genomes

# combine all the genomic stats together
rm -rf $BASE/total.genome.stats.txt
cat *_n50.txt | awk '!seen[$0]++' > $BASE/total.genome.stats.txt
rm *_n50.txt


cd $BASE/curated_assembly

cp /Volumes/Bennett_BACKUP/Research/curtobacterium/curto-evolution/reference_genome/MMLR14002_pacbio-corrected.fasta $BASE/curated_assembly/

ani-aai-matrix.sh -i .fasta -m ani -t 8


