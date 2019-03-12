#!/bin/bash

bgzip -cd $1 | grep -v "^#" |  awk '{print $1"-"$2"-"$4"-"$5}' | sort 