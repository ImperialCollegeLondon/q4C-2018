#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

# script for alignment using BWA. the -p flag is required for 
# interleaved fastq, the output from cramtools cram to fastq conversion
### See http://bio-bwa.sourceforge.net/bwa.shtml


# aliases

# files
BWAPATH=<PATH/TO/BWA>
REFPATH=<PATH/TO/REFERENCE/GENOME>
ALIGNPATH=<PATH/TO/ALIGHMENT/OUTPUT>

# software
BWA="$BWAPATH mem -p -M -t 8 $REFPATH "

# files to process (passed in from the job index in a serial job)
FILENAME=<PATH/TO/INPUT/FILE.fastq>
fbname=$(basename "$FILENAME" .fastq)

# align reads using bwa
COMMAND="$BWA $FILENAME > $ALIGNPATH/$fbname.sam"
eval "$COMMAND"



