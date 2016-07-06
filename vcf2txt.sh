#!/bin/bash

gunzip -c $1 | grep -v "^#" |  awk '{print $1"-"$2"-"$4"-"$5}'