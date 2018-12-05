#!/bin/bash
#$ -V
#$ -N Htseq
#$ -m ea
#$ -o /workdir/jobs/htseq.stdout
#$ -e /workdir/jobs/htseq.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20


# Count the number of reads mapping to each feature
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
#Option
Option=$2
#Get input direction
Input_Dir=$3
#Get output direction
Output_Dir=$4
#Get reference direction
Reference_Dir=$5




hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2


#Htseq-count
if (( $Option == 1 )) 
then
  /opt/python/bin/htseq-count -i gene_id --mode=intersection-nonempty --nonunique=none --format=bam $Input_Dir/${FastqID}/accepted_hits.bam $Reference_Dir/Mus_musculus38.gtf > $Output_Dir/${FastqID}_nonempty_count_table.txt
else if (( $Option == 2 )) 
then
  /opt/python/bin/htseq-count -i gene_id --mode=union --nonunique=none --format=bam $Input_Dir/${FastqID}/accepted_hits.bam $Reference_Dir/Mus_musculus38.gtf > $Output_Dir/${FastqID}_nonempty_count_table.txt
else if (( $Option == 3 )) 
then
  /opt/python/bin/htseq-count -i gene_id --mode=intersection-strict --nonunique=none --format=bam $Input_Dir/${FastqID}/accepted_hits.bam $Reference_Dir/Mus_musculus38.gtf > $Output_Dir/${FastqID}_nonempty_count_table.txt
fi

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'