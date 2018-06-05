#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#!/bin/bash

# This script is submitted from an external loop, 
# once for each trimmed read 1 file name 
# (outputdir and filename variable passed from outside)


filename2=`echo $filename  | sed 's#R1_001#R2_001#g; s#_val_1#_val_2#g' - `
echo "File2 is "$filename2

output1=$outputdir"NlaIIIsummary1.tsv"
output2=$outputdir"NlaIIIsummary2.tsv"

echo $outputdir
echo $filename
echo $filename2

date
cp $filename .
cp $filename2 .

gunzip *.gz

ls -l -h -R


# digest read 1 and read 2 separately
for read1 in *_R1_*
do

read2=`echo $read1  | sed 's#R1_001#R2_001#g; s#_val_1#_val_2#g' - `


echo "Read1 is "$read1
echo "Read2 is "$read2

firstline=`head -1 $read1 | awk -F ":" '{print $1}'`

perl FC02_digestseqRead1advanced.pl $read1 $output1 $firstline
perl FC02_digestseqRead2.pl $read2 $output2 $firstline

date
echo "  "
done

ls


# compress and copy back to working directory
gzip *NlaIII.fq

cp *NlaIII* $outputdir

	