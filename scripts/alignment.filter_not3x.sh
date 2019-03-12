#!/bin/bash

cat $1 |  awk '{name=$0; getline; if (length($0)%3 == 0){print name;print $0;}}'