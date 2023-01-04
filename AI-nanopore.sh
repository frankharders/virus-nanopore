#!/bin/bash

##  list of different scripts for assembly of a draft genome


## create a sample file and concatenate all raw output files into a new file with the sample name for downstream processing
#./00_structure.sh

## adapter trimming with porechop 
#./01_porechop.sh

## construct a draft genome
./02_assembly.sh



exit 1
