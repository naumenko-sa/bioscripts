#!/bin/bash

# multiz: download, make, and add to PATH
# wget https://www.bx.psu.edu/miller_lab/dist/multiz-tba.012109.tar.gz
#
# Kent utilities - add to PATH
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/axtToMaf
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faOneRecord

http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
http://hgdownload.soe.ucsc.edu/goldenPath/dm6/vsDroSim1/dm6.droSim1.net.axt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/droSim1/bigZips/droSim1.chrom.sizes
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes

axtToMaf -tPrefix=dm6. -qPrefix=droSim1. -tSplit dm6.droSim1.net.axt dm6.chrom.sizes droSim1.chrom.sizes dm6.droSim1.maf_dir

cat dm6.chrom.sizes  | grep -v chrUn | grep -v random | awk '{print $1}' > dm6.chrom.list
cat dm6.fa | sed s/"chr"/"dm6.chr"/ > dm6.chr.fa

for chr in `cat dm6.chrom.list`
do
    faOneRecord dm6.chr.fa dm6.$chr > dm6.$chr.fa
done

cd dm6.droSim1.maf_dir

rm *Un*
rm *random*

for f in *.maf;
do
    fasta=`echo $f | sed s/maf/fa/`
    maf2fasta ../dm6.$fasta $f fasta2 iupac2n > dm6.droSim1.$f.fa
done

cd ..