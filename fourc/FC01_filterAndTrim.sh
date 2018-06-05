#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#!/bin/bash

# define path
FASTQPATH= <PATH/TO/FASTQ/FILES>

# files to process (passed in from the job index in a serial job)
FILENAME=<PATH/TO/INPUT/FILE.fastq.gz>

f=$FILENAME
f2=`echo $f  | sed 's#_R1#_R2#g' - ` 

filename=$(basename "$f" .gz)
filename2=$(basename "$f2" .gz)

cp $FASTQPATH/$f .
cp $FASTQPATH/$f2 .

gunzip *.gz


# search and filter for first 4b , those will be trimmed out
echo "filter1..."
~/Applications/cutadapt-1.8.3/bin/cutadapt \
-g ^CGGG \
-e 0 --overlap 4 \
--discard-untrimmed \
-o filter1_1.fastq -p filter1_2.fastq $filename $filename2
 
# search and filter for CATG - not trimmed. 
echo "filter2..."
~/Applications/cutadapt-1.8.3/bin/cutadapt \
-g ^CATG \
-e 0 --overlap 4 \
--no-trim \
--discard-untrimmed \
-o $filename -p $filename2 filter1_1.fastq filter1_2.fastq 

# trim
echo "trimming..."
trim_galore --paired -q 20 -a AGATCGGAAGAGC -a2 CCCGGAGGACCTG --no_report_file $filename $filename2

gzip *_val_*

cp *_val_* $FASTQPATH

