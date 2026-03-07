#!/bin/bash

repo=/home/neuro/Documents/NMD_analysis/other_repos/REMBRANDTS
file=/home/neuro/Documents/NMD_analysis/Analysis/NMD-analysis/data/data.txt
inputdir=/home/neuro/Documents/NMD_analysis/Analysis/Results/EISA

bash ${repo}/REMBRANDTS.sh NMD_data $file $inputdir .99 linear

