#!/bin/bash
#$ -V
#$ -N IndexReferenceGENOME
#$ -m ea
#$ -o /work/jobs/indexReferenceGenome.stdout
#$ -e /work/jobs/indexReferenceGenome.stderr
#$ -w e
##$ -l h_vmem=3G
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

#Get reference direction
Reference_Dir=$1


#.fa file and .gtf file
/opt/bowtie2-2.2.3/bowtie2-build $Reference_Dir/Mus_musculus38.fa $Reference_Dir/Mus_musculus38
#index with samtools
samtools faidx $Reference_Dir/Mus_musculus38.fa

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'