#!/bin/bash

# download a wgs run from sftp://tcagfts.ccm.sickkids.ca
# $1 = credentials: user:password
# $2 = RUNID, i.e. HIR9438
# $3 = dir to download in RUNID, i.e. 157639_files

# run on data2
# curl sftp is required
module load curl

mkdir $3
cd $3

url="sftp://tcagfts.ccm.sickkids.ca/inbox"

# download files in the top dir

# for debugging use -r 0-100 - download first chunk
# curl_cmd="curl -r 0-10"
curl_cmd="curl"

for f in `curl -u $1 -l -s $url/$2/$3/`
do
    $curl_cmd -u $1 -O -s $url/$2/$3/$f
done

flowcell=`curl -u $1 -l -s $url/$2/$3/fastq/ | tail -n1`
echo $flowcell

mkdir -p fastq/$flowcell
cd fastq/$flowcell
for f in `curl -u $1 -l -s $url/$2/$3/fastq/$flowcell/`
do
    $curl_cmd -u $1 -O -s $url/$2/$3/fastq/$flowcell/$f
done
cd ../../

for d in {annot,gatk,sv,cnv}
do
    mkdir $d
    cd $d
    for f in `curl -u $1 -l -s $url/$2/$3/$d/`
    do
	$curl_cmd -u $1 -O -s $url/$2/$3/$d/$f
    done
    cd ..
done

cd ..
