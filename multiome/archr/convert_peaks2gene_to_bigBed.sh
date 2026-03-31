#! /bin/bash

bedSort $1 $1
bedToBigBed $1 ~/nest/assembly/lonStrDom2/ucsc/chrom.sizes.ucsc "${1%.bed}.bb" -as=interact.as -type=bed5+13
