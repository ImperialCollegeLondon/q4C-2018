#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#@ script for alignment using GSNAP. 


# aliases

# software
GSNAPPATH=<PATH/TO/GSNAP>

# files
RSPATHALIGN=<PATH/TO/GSNAP/OUTPUT/>
RSPATHRAW=<PATH/TO/FASTQ/FILES>


# files to process (passed in from the job index in a serial job)
FILENAME=<PATH/TO/READ1/FILE_1.fastq.gz>
fbname=$(basename "$f" _1.fastq.gz)

# define read1, read2
READ1=$RSPATHRAW/$fbname"_1.fastq.gz"
READ2=$RSPATHRAW/$fbname"_2.fastq.gz"

# run gsnap
COMMAND="$GSNAPPATH --split-output $RSPATHALIGN/$fbname \
-m 10 \
-t 8 \
-A sam \
--gunzip \
-d AB513134_hg19 \
-D ~/data/Ref \
-s AB513134_hg19.splicesites \
-N 1 \
-v snp138_gsnap \
--use-sarray=0 \
$READ1 $READ2 "

eval $COMMAND
