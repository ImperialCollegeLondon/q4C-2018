#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018


###########################################################################
# base recalibration using GATK BQSR
ALIGNPATH=<PATH/TO/ALIGNED/FILES>
BAMPATH=<PATH/TO/BAM/FILES>
GATKPATH=<PATH/TO/GenomeAnalysisTK.jar>

REF=<PATH/TO/FASTA/REFERENCE> # using a concatanated reference hg19 + HTLV-1
VAR=<PATH/TO/DBSNP/FILE> 

GATK="java -jar $GATKPATH "

# files to process (passed in from the job index in a serial job)
FILENAME=<PATH/TO/INPUT/FILE.bam>
fbname=$(basename "$FILENAME" .bam)

###########################################################################
# base calibration - create calibrator

COMMAND="$GATK -T BaseRecalibrator \
    -R $REF \
    -knownSites $VAR \
    -I $FILENAME \
    -o $BAMPATH/$fbname.recal"

eval "$COMMAND"

###########################################################################
# base calibration
echo "*************Base Calibration**************"

COMMAND="$GATK -T PrintReads \
   -R $REF \
   -I $FILENAME \
   -BQSR $BAMPATH/$fbname.recal \
   -o $BAMPATH/$fbname.c.bam"

eval "$COMMAND"

###########################################################################
# rename files

echo "*************Renaming files**************"

COMMAND="mv $BAMPATH/$fbname.c.bam $BAMPATH/$fbname.bam"
eval "$COMMAND"

COMMAND="mv $BAMPATH/$fbname.c.bai $BAMPATH/$fbname.bai"
eval "$COMMAND"