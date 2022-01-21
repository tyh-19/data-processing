#!/bin/bash

taxolevel=$2
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_microbe/$1/microbe_${taxolevel}_level.txt"
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_microbe/$1/report"

python /data/taoyuhuan/projects/exOmics_RNA/level_3_microbe/summarize-kraken.py -i ${input} -l ${taxolevel} -o ${output}
 
