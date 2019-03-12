#!/bin/bash

#percent of columns with gaps in an alignment

a=`cat $1 | grep -v '>' | sed s/[A-Z]//g | awk '{sum+=length($0)}END{print sum}'`
b=`cat $1 | grep -v '>' |  awk '{sum+=length($0)}END{print sum}'`

echo $((a*100/b))

