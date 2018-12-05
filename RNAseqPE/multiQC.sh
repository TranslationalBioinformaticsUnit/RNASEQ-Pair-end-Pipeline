#!/bin/bash
#$ -V
#$ -N MultiQC
#$ -m ea
#$ -o /workdir/jobs/multiQC.stdout
#$ -e /workdir/jobs/multiQC.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20


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

#Get input direction
Input_Dir=$1
#Get output direction
Output_Dir=$2


hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2


/opt/python/bin/multiqc $Input_Dir/ -o $Output_Dir/

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'