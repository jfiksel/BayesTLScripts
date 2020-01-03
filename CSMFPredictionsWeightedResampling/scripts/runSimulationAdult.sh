#!/bin/bash
#$ -l mem_free=3G
#$ -l h_vmem=3G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-8000
Rscript runSimulationAdult.R $SGE_TASK_ID
