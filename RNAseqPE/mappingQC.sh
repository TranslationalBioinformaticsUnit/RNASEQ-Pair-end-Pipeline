#!/bin/bash
#$ -V
#$ -N MappingQC
#$ -m ea
#$ -o /workdir/jobs/mappingQC.stdout
#$ -e /workdir/jobs/mappingQC.stderr
#$ -w e
##$ -l h_vmem=2G
#$ -l m_core=20


# Quality control of the mapping
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
# Get the fastq id file
FastqID=$1
#Get input direction
Input_Dir=$2
#Get output direction
Output_Dir=$3
#Get reference direction
Reference_Dir=$4


/opt/jdk1.8.0_161/bin/java -Xmx18000m -Djava.io.tmpdir=`pwd`/tmp -jar /opt/picard/build/libs/picard.jar CollectAlignmentSummaryMetrics \
		I= $Input_Dir/${FastqID}/accepted_hits.bam \
		O= $Output_Dir/${FastqID}/${FastqID}.metrics.txt \
		REFERENCE_SEQUENCE= $Reference_Dir/Mus_musculus38.fa \
		VALIDATION_STRINGENCY=LENIENT


printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'