#!/bin/bash
#PBS -l walltime=23:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

date
# nohups are dying on qlogin nodes, data nodes are better for long data installation runs
# data2 has modules, data7 does not.
######################################################################
# fresh install of new bcbio instance:
# 1. Don't mix with old environments
# mv ~/.conda/environments.txt ~/.conda/environments.default.txt - move back
# export PATH=/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/bin
# export PYTHONPATH=
# wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
# echo "Installing to " $1
# module load python/2.7.12 - not working, just plain python from the system is better
# python bcbio_nextgen_install.py $1 --tooldir $1 --genomes GRCh37 --aligners bwa
 #--isolate --nodata
######################################################################
# 2. Use the new environment: 
# create a .test_profile:
# export PATH=$HOME/cre:$HOME/crt:$HOME/crg:/hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/anaconda/bin:$HOME/tools/mc-4.8.16/bin:$HOME/jkent_tools:$HOME/bioscripts:.:/usr/local/bin:/opt/moab/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin
# export PYTHONPATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/anaconda/lib/python3.6
. /hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/.profile115
which python
echo $PYTHONPATH
######################################################################
# 3. Upgrade tools. If tooldir was set before, no need to specify it again
which bcbio_nextgen.py
# bcbio_nextgen.py upgrade -u skip --tools
#--tooldir $1
######################################################################
# 4. Install indices
# genomes = {GRCh37, hg38}
# bcbio_nextgen.py upgrade -u skip --genomes hg38 --aligners bwa --cores 10
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --aligners star --cores 10
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --aligners hisat2 --cores 10
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --aligners rtg --cores 10
#########################################################################
# 5. upgrade code to the latest stable version
# bcbio_nextgen.py upgrade -u stable
# upgrade code to development
# bcbio_nextgen.py upgrade -u development
#########################################################################
# 6. data installation/upgrade 
# GRCh37, hg38, mm10
# data installation takes a lot of time (gnomad, dbnsfp) it is better to have data and just upgrade bcbio code
# 1.1.5 - a huge update to python3, installed from scratch

# upgrades data installed before (gemini, cadd) for all references
# bcbio_nextgen.py upgrade --data --genomes GRCh37

# VEP is upgraded quite often ~2-3 months - when upgrading tools it looks for new VEP cache
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget vep

# gemini ~3h for GRCh37
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget gemini
# --genomes hg38

# cadd is in dbnsfp
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget cadd
# --genomes hg38

# gnomad 14h
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget gnomad
# bcbio_nextgen.py upgrade -u skip --genomes hg38 --datatarget gnomad

# dbnsfp
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget dbnsfp

# rnaseq
bcbio_nextgen.py upgrade -u skip --genomes hg38 --datatarget rnaseq

######################################################################
# fresh installation for Sam with human and mouse genome

# to check what enviroments were picked up during the installation
# conda info --envs --json
# check file ~/.conda/environments.txt - if it has environments from all installations they could interfere
# wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
# export PYTHONPATH=/hpf/largeprojects/lauryl/bcbio110/anaconda/lib/python2.7

# PATH=/hpf/largeprojects/lauryl/bcbio110/anaconda/bin
# PATH=${PATH}:/usr/local/bin:/opt/moab/bin:/home/naumenko/cre:/home/naumenko/crt:/home/naumenko/crg:/home/naumenko/tools/mc-4.8.16/bin:/home/naumenko/jkent_tools
# PATH=${PATH}:/home/naumenko/bioscripts:.:/home/naumenko/.aspera/connect/bin:/usr/local/bin:/usr/lib64/qt-3.3/bin:/opt/moab/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin
# PATH=${PATH}:/opt/ibutils/bin:/sbin:/usr/sbin:/sbin:/usr/sbin
# export PATH

# python bcbio_nextgen_install.py /hpf/largeprojects/lauryl/bcbio110 --tooldir=/hpf/largeprojects/lauryl/bcbio110 --genomes mm10 --aligners bwa --isolate
# bcbio_nextgen.py upgrade -u skip --tools --tooldir /hpf/largeprojects/lauryl/bcbio110
# bcbio_nextgen.py upgrade -u skip --data --genomes mm10 --datatarget variation --datatarget vep

date
