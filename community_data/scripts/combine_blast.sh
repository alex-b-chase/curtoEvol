#!/bin/bash

BASEDIR=/bio/abchase/MG/elevation_transplant


for timepoint in initialtime finaltime
do
	cd $BASEDIR/$timepoint

	ls *.filter.total.faa | sed 's/.filter.total.faa//g' > temp.txt

	while read line
	do
	
		cat ${line}*.blast.txt > ${line}.blast.total.txt
		wc -l ${line}.blast.total.txt

	done < temp.txt

	rm temp.txt
	rm -f *bP.e*
	rm -f *bP.o*

done
