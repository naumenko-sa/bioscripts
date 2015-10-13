#!/bin/bash

#get protein name by protein gi

efetch -db protein -id "$1" -format gpc | xtract -element "INSDSeq_definition"