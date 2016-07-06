#!/bin/bash
#only sites with variants because variants are usually the same, only few like A - ATTT in one sample  A - A/ATTTT in other

gunzip -c $1 | grep -v "^#" |  awk '{print $1"-"$2}'