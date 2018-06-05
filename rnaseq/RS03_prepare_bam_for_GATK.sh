#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#@ use picard to prepare a bam file for processing by GATK (e.g. indexing, add group information to header)

# aliases

# software
PICARDPATH=<PATH/TO/PICARD.jar>
PICARD="java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PICARDPATH "

# files
RSPATHBAM=<PATH/TO/BAM/FILES/>
SAMPLELIST=<PATH/TO/METADATA/FILE> # list of clones. see suppl file 1 as example. 

# files to process (passed in from the job index in a serial job)
FILENAME=<PATH/TO/INPUT/FILE.bam>
fbname=$(basename "$FILENAME" concordant_uniq.bam)


###########################################################################
# remove duplicates using picard
MARKDUPS="MarkDuplicates \
I=$RSPATHBAM/$fbname.concordant_uniq.bam \
O=$RSPATHBAM/$fbname.dd.bam \
REMOVE_DUPLICATES="true" \
M=$RSPATHBAM/$fbname.marked_dup_metrics.txt"

COMMAND=$PICARD$MARKDUPS

eval "$COMMAND"

###########################################################################
# rename file
COMMAND="mv $RSPATHBAM/$fbname.dd.bam $RSPATHBAM/$fbname.bam"
eval "$COMMAND"

###########################################################################
# organize RG line
ADDREPLACEREADGREOUPS="$PICARD AddOrReplaceReadGroups "

smname=`awk -v items="$fbname" -F"[\t]" '$1 ~ items {print $2}' $SAMPLELIST \
| uniq | sed s#' '#_#g -`

COMMAND="$ADDREPLACEREADGREOUPS I=$RSPATHBAM/$fbname.bam "
COMMAND="$COMMAND O=$RSPATHBAM/$fbname.t.bam "
COMMAND="$COMMAND RGID=1 "
COMMAND="$COMMAND RGLB=HMMVHBBXX "
COMMAND="$COMMAND RGPL=illumina "
COMMAND="$COMMAND RGPU=$fbname "
COMMAND="$COMMAND RGSM=$smname"

eval "$COMMAND"


###########################################################################
# rename file
COMMAND="mv $RSPATHBAM/$fbname.t.bam $RSPATHBAM/$fbname.bam"
eval "$COMMAND"


###########################################################################
# create index for bam file
COMMAND="$PICARD BuildBamIndex "
COMMAND="$COMMAND I=$RSPATHBAM/$fbname.bam"
eval "$COMMAND"








