#!/bin/bash
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-3000
module load java/1.8.0
export JAVA_HOME=/usr/lib/jvm/java-1.8.0/jre
export PATH=$PATH:$JAVA_HOME/bin
export JAVA_OPTS="-Xmx2048m -XX:CompressedClassSpaceSize=256m"
export _JAVA_OPTIONS="-Xmx2048m -XX:CompressedClassSpaceSize=256m"
export _JAVA_OPTIONS="-XX:+UseSerialGC"
Rscript SimulationChildExpandedCauses.R $SGE_TASK_ID
