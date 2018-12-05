#!/bin/bash

#/example/fileId.txt
SamplesNames_Dir=$1
#/example/fastq -> where are the fastq files
Fastq_Dir=$2
#/example/workdir
Work_Dir=$3
#/example/scripts
Code_Dir=$4
#/example/reference (where is the reference genome and gtf file)
Reference_Dir=$5

#INDEX REFERENCE GENOME (mus_musculus38)
#Parameters(1): ReferenceDir
#qsub -N Index_Genome $Code_Dir/indexReferenceGenome.sh $Reference_Dir


#for each fastq
for FastqID in $(cat $SamplesNames_dir); 
  do
    echo "$FastqID"
    ###############
    #1-Quality control of the original fastq
    ###############
    #Create folder fastqQC
    if [ ! -d $Work_Dir/fastqQC ]; then
      mkdir $Work_Dir/fastqQC
      chmod 777 $Work_Dir/fastqQC
    fi
    #Parameters(4): FastqID Option(1) InputDir OutputDir
    ###Option=1 -> fastqc of original fastq
    ###Option=2 -> fastqc of trimmed fastq
    qsub -N QC_${FastqID} $Code_Dir/fastqQC.sh $FastqID 1 $Fastq_Dir $Work_Dir/fastqQC


    
    ##############
    #2-Trimmed Reads
    ##############
    #Create folder trimmed_reads
    if [ ! -d $Work_Dir/trimmed_reads ]; then
      mkdir $Work_Dir/trimmed_reads
      chmod 777 $Work_Dir/trimmed_reads
    fi
    #Parameters(3): FastqID InputDir OutputDir
    #Don't start this job until fastqc is not finish 
    qsub -hold_jid QC_${FastqID} -N trim_${FastqID} $Code_Dir/trimmedReads.sh $FastqID $Fastq_Dir $Work_Dir/trimmed_reads
 
 
    
    ###############
    #3-Quality control of the TRIMMED fastq
    ###############
    #Create folder trimmed_fastqQC
    if [ ! -d $Work_Dir/trimmed_reads/trimmed_fastqQC ]; then
      mkdir $Work_Dir/trimmed_reads/trimmed_fastqQC
      chmod 777 $Work_Dir/trimmed_reads/trimmed_fastqQC
    fi
    #Parameters(4): FastqID Option(2) InputDir OutputDir
    #Don't start this job until trimming is not finish 
    qsub -hold_jid trim_${FastqID} -N QC_trim_${FastqID} $Code_Dir/fastqQC.sh $FastqID 2 $Work_Dir/trimmed_reads $Work_Dir/trimmed_reads/trimmed_fastqQC
 
 
  
    ################
    #4-Mapping TopHAT
    ################
    #Create folder bam
    if [ ! -d $Work_Dir/bam ]; then
      mkdir $Work_Dir/bam
      chmod 777 $Work_Dir/bam
    fi
    #Parameters(4): FastqID InputDir OutputDir ReferenceDir
    #Don't start this job until fastqc trimming is not finish 
    qsub -hold_jid QC_trim_${FastqID} -N mapping_${FastqID} $Code_Dir/mappingTophat.sh $FastqID $Work_Dir/trimmed_reads $Work_Dir/bam $Reference_Dir



    ###############
    #5-Quality control of the mapping
    ###############
    #Parameters(4): FastqID InputDir OutputDir (the inputdir and oputputdir is the same) ReferenceDir
    #Don't start this job until mapping is not finish 
    qsub -hold_jid mapping_${FastqID} -N mappingQC_${FastqID} $Code_Dir/mappingQC.sh $FastqID $Work_Dir/bam $Work_Dir/bam $Reference_Dir



    ################
    #6-HTSEQ Count Tables
    ################
    #Create folder count_tables
    if [ ! -d $Work_Dir/count_tables ]; then
      mkdir $Work_Dir/count_tables
      chmod 777 $Work_Dir/count_tables
    fi
    #Parameters(5): FastqID Mode_option InputDir OutputDir ReferenceDir
    ###Mode_option=1 -> intersection-nonempty(by default)
    ###Mode_option=2 -> union
    ###Mode_option=3 -> intersection_strict
    #Don't start this job until quality control of mapping is not finish 
    qsub -hold_jid mappingQC_${FastqID} -N count_${FastqID} $Code_Dir/htseq.sh $FastqID 1 $Work_Dir/bam $Work_Dir/count_tables $Reference_Dir
    
    
    
    ################
    #7-MultiQC Plots
    ################
    if [ ! -d $Work_Dir/multiQCPlots ]; then
      mkdir $Work_Dir/multiQCPlots
      chmod 777 $Work_Dir/multiQCPlots
    fi
    #MultiQC plots of the mapping quality control results
    #Parameters(2): InputDir OutputDir
    #Don't start this job until the generating of count tables is not finish 
    qsub -hold_jid count_${FastqID} -N multiQC_${FastqID} $Code_Dir/multiQC.sh $Work_Dir/bam $Work_Dir/multiQCPlots
    
  done 