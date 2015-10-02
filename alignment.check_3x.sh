#!/bin/bash

cat $1 | grep  -v '>' | awk '{print length($0)/3}' | grep '\.' | wc -l