#! /bin/bash
#$ -q bio,abio

module load wgs/8.2

PBcR -length 500 -partitions 200 -l mmlr15002 -s pacbio.spec -fastq R088-C01.1-Reads.fastq genomeSize=3700000

# for reference - http://wgs-assembler.sourceforge.net/wiki/index.php/PBcR


# this crashes on the server - doesn;t fully run all the way through
# more info here - https://sourceforge.net/p/wgs-assembler/bugs/338/
# doesn't seem like they support this anymore for cluster use...

# reduce the memory on the local machine to 12 mb ram needed - include the pabcio.spec file
######
## limit to 32GB. By default the pipeline will auto-detect memory and try to use maximum. This allow limiting it
#ovlMemory = 12
#ovlStoreMemory= 12000
#merylMemory = 12000

# need to install QUIVER
# reference - https://github.com/PacificBiosciences/GenomicConsensus

# first, install boost
# brew install boost
# brew install swig

# download the dependencies for QUIVER
# https://github.com/PacificBiosciences/ConsensusCore
# ./configure  [see ./configure --help for details]
# make csharp
# sudo python setup.py install

# https://github.com/PacificBiosciences/ConsensusCore2
# ???? nothing in README

# https://github.com/PacificBiosciences/pbcore
# sudo pip install -r requirements.txt
# sudo python setup.py install

# cd ../GenomicConsensus-master
# sudo python setup.py install

# QIVER is ready to go!!!

# need pbalign now...
