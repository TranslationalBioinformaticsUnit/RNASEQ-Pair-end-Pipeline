#!/bin/bash
#$ -V
#$ -N QualControl
#$ -m ea
#$ -o /workdir/jobs/fastqc.stdout
#$ -e /workdir/jobs/fastqc.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=10

# Perform quality control of sequencing reads

abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

set -e
trap 'abort' 0

# Get the fastq id file
FastqID=$1
# Get the option
Option=$2
#Get input direction
Input_Dir=$3
#Get output direction
Ouput_Dir=$4

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2


if (( $Option == 1 )) 
then
  # Run quality control of the original fastq using FastQC
  /opt/FastQC/fastqc -o $Ouput_Dir/ -t 10 $Input_Dir/${FastqID}_R1_001.fastq.gz $Input_Dir/${FastqID}_R2_001.fastq.gz
else if (( $Option == 2 ))
then
  #Quality Control del trimmed fastq
  /opt/FastQC/fastqc -o $Output_Dir/ -t 10 $Input_Dir/${FastqID}_1_trimmed.fastq $Input_Dir/${FastqID}_2_trimmed.fastq
fi


printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'