#!/bin/bash

#get gam species list from tree

cat $1 | sed -r s/":"/"\n"/g  | grep gam | grep -o gam.*@.*$

