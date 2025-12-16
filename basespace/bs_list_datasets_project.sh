#!/bin/bash

/data/tools/bs list datasets -F Name -F Id -F DataSetType.Id -F TotalSize -F Project.Name --project-name=$1