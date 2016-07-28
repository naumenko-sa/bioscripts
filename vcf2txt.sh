#!/bin/bash

cat $1 | grep -v "^#" |  awk '{print $1"-"$2"-"$4"-"$5}'