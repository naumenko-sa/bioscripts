#!/bin/bash

bs list biosamples list -F BioSampleName -F DefaultProject.Name -f csv > biosamples_runs.csv
