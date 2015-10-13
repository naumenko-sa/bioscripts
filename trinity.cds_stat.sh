#!/bin/bash
#makes cds stat table using gam species after trinity_to_cds.sh
#for f in {gam0,gam1,gam1.1,gam2,gam2quest,gam3,gam3.1,gam4,gam4.1,gam4.2,gam5,gam5.6,gam6,gam7,gam7.3,gam8,gam8.1,gam8.4,gam8.5,gam8.7,gam9,gam9.4,gam9.5,gam9.9,gam10,gam10.1,gam10.2,gam11.1,gam11.4,gam12.1,gam13.4,gam14.1,gam15.1,gam16.4,gam19.2,gamk1,gam22.1,gam27.1,gam27.2,gam27.3,gam32.1,gam34.2,gam35.3,gam35.6,gam39.1,gam43.1,phaw,dpulex};
#do
f=gam35.6.newbler
    num_cds=`grep -c '>' $f.fasta.cds`;
    cds_avr_len=`cat $f.fasta.cds | grep -v '>' | awk '{sum+=length($0);}END{print sum/NR}'`
    cds_total=`cat $f.fasta.cds | grep -v '>' | awk '{sum+=length($0);}END{print sum}'`;
    cds_complete=`grep -c 'complete' $f.fasta.cds`;
    cds_complete_avr=`cat $f.fasta.cds | awk '{name=$0;getline; if (name ~ /complete/) {cds_num++;sum+=length($0);}}END{print sum/cds_num}'`;
    cds_complete_len=`cat $f.fasta.cds | awk '{name=$0;getline; if (name ~ /complete/) sum+=length($0);}END{print sum}'`;
    echo -e $f"\t"${num_cds}"\t"$cds_avr_len"\t"$cds_total"\t"$cds_complete"\t"$cds_complete_avr"\t"$cds_complete_len;
#done;