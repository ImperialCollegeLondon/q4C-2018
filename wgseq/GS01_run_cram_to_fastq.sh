#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#@ convert cram data to fastq files using cramtools

# aliases

# files
FASTQPATH=<PATH/TO/FASTQ/OUTPUT>

# software
CRAMTOOLSPATH=<PATH/TO/CRAMTOOLS.jar>
CRAMTOFASTQ="java -jar $CRAMTOOLSPATH fastq -I "

# files to process (passed in from the job index in a serial job)
FILENAME=<PATH/TO/INPUT/FILE.cram>
fbname=$(basename "$FILENAME" .cram)

# convert cram to fastq using cramtools
COMMAND="$CRAMTOFASTQ $FILENAME > $FASTQPATH/$fbname.fastq"

eval "$COMMAND"

