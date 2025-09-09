#!/bin/bash

for f in *kraken.txt;do echo $f; cat $f | awk '{if ($1>0.01) print $0}';done;