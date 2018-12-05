#!/bin/bash
#$ -V
#$ -N MappingTophat
#$ -m ea
#$ -o /workdir/jobs/mapping.stdout
#$ -e /workdir/jobs/mapping.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20


# Mapping with Mapping Tophat program

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
#Get input direction
Input_Dir=$2
#Get output direction
Output_Dir=$3
#Get reference direction
Reference_Dir=$4

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

#Create folder inside bam folder called 'FastqID'
#Where saved all the outputs of tophat
if [ ! -d $Input_Dir/$FastqID ]; then
  mkdir $Input_Dir/$FastqID
  chmod 777 $Input_Dir/$FastqID
fi

# Run Mapping tophat (reference_genome.fa)
/opt/tophat-2.1.1.Linux_x86_64/tophat2 --no-coverage-search  -G $Reference_Dir/Mus_musculus38.gtf -o $Output_Dir/$FastqID $Reference_Dir/Mus_musculus38 $Input_Dir/${FastqID}_1_trimmed.fastq $Input_Dir/${FastqID}_2_trimmed.fastq


printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'