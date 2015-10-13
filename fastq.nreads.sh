#!/bin/bash

cat $1 | wc -l | awk '{print $0/4}' > $1.nreads