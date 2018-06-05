#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018


###########################################################################
# base recalibration using GATK BQSR

# aliases

# files
REF=<PATH/TO/FASTA/REFERENCE.fa> # using a concatanated reference hg19 + HTLV-1
GATKPATH=<PATH/TO/GenomeAnalysisTK.jar>
OUTPATH=<PATH/TO/GVCF/OUTPUT/>

# software
GATK="java -jar $GATKPATH "

# files to process (passed in from the job index in a serial job)
FILENAME=<PATH/TO/INPUT/FILE.bam>
fbname=$(basename "$FILENAME" .bam)

echo fbname is $fbname

# run base recalibration using gatk
COMMAND="$GATK -T HaplotypeCaller  "
COMMAND="$COMMAND -R $REF "
COMMAND="$COMMAND -nct 8 "
COMMAND="$COMMAND -I $FILENAME "
COMMAND="$COMMAND -o $OUTPATH/$fbname.g.vcf"
COMMAND="$COMMAND --emitRefConfidence GVCF"
COMMAND="$COMMAND --variant_index_type LINEAR"
COMMAND="$COMMAND --variant_index_parameter 128000"

eval "$COMMAND"

