#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#@ script for alignment using GSNAP. 
#@ Note: gsnap output is separated by type of alignment result - use here the concordant uniq 
#@ Note2: if samples were sequenced across multiple lanes, merge and rename appropriately. 


# alias

# software
PICARDPATH=<PATH/TO/PICARD.jar>
PICARD="java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PICARDPATH "

# files
RSPATHBAM=<PATH/TO/BAM/OUTPUT/>

# files to process (passed in from the job index in a serial job)
FILENAME=<PATH/TO/ALIGEND/FILE.concordant_uniq>
fbname=`echo $FILENAME  | sed 's#_S.*$##g' - ` 

outfile=$fbname.bam

# use picard to 
COMMAND="$PICARD MergeSamFiles \
I=$FILENAME \
O=$RSPATHBAM/$outfile"

eval $COMMAND


