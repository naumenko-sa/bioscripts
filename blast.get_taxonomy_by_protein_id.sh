#!/bin/bash

#get taxonomy by protein gi

efetch -db protein -id "$1" -format gpc | xtract -element "INSDSeq_taxonomy"