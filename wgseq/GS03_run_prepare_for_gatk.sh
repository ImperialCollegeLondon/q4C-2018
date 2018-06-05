#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

# using picard to prepare bam files for downstream gatk analysis (including indexing and adding group data)

# aliases

# software
PICARDPATH=<PATH/TO/PICARD.jar>
PICARD="java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PICARDPATH SortSam SO=coordinate "

# files
BAMPATH=<PATH/TO/BAM/OUTPUT>

# files to process (passed in from the job index in a serial job)
FILENAME=<PATH/TO/INPUT/FILE.sam>
fbname=$(basename "$FILENAME" .sam)

###########################################################################
# sort sam and convert to bam
COMMAND="$PICARD INPUT=$FILENAME OUTPUT=$BAMPATH/$fbname.bam \
VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=false"

eval "$COMMAND"

###########################################################################
# remove duplicates using picard
PICARD="java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PICARDPATH "
MARKDUPS="MarkDuplicates \
I=$BAMPATH/$fbname.bam \
O=$BAMPATH/$fbname.dd.bam \
REMOVE_DUPLICATES="true" \
M=$BAMPATH/$fbname.marked_dup_metrics.txt"

COMMAND=$PICARD$MARKDUPS

eval "$COMMAND"

###########################################################################
# rename file
COMMAND="mv $BAMPATH/$fbname.dd.bam $BAMPATH/$fbname.bam"
eval "$COMMAND"

###########################################################################
# organize RG line
ADDREPLACEREADGREOUPS="java -jar $PICARDPATH \
AddOrReplaceReadGroups "

IFS='# ' read -r -a name_items <<< "$f"
IFS='_ ' read -r -a items <<< "${name_items[0]}"

COMMAND="$ADDREPLACEREADGREOUPS I=$BAMPATH/$fbname.bam "
COMMAND="$COMMAND O=$BAMPATH/$fbname.t.bam "
COMMAND="$COMMAND RGID=1 "
COMMAND="$COMMAND RGLB=DN475270M "
COMMAND="$COMMAND RGPL=illumina "
COMMAND="$COMMAND RGPU=${items[0]} "
COMMAND="$COMMAND RGSM=$fbname"

eval "$COMMAND"

###########################################################################
# rename file
COMMAND="mv $BAMPATH/$fbname.t.bam $BAMPATH/$fbname.bam"
eval "$COMMAND"

###########################################################################
# create index for bam file
COMMAND="java -jar $PICARDPATH BuildBamIndex "
COMMAND="$COMMAND I=$BAMPATH/$fbname.bam"
eval "$COMMAND"

