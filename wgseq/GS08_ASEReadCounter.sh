#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018


# aliases

# files
REF=<PATH/TO/FASTA/REFERENCE.fa> # using a concatanated reference hg19 + HTLV-1
VAR=<PATH/TO/DBSNP/FILE> 
CSVPATH=<PATH/TO/CSV/OUTPUT/>
VCF=<PATH/TO/VCF/>
FILENAME=<PATH/TO/INPUT/FILE.bam>

# software
GATKPATH=<PATH/TO/GenomeAnalysisTK.jar>
GATK="java -jar $GATKPATH "

# files to process (passed in from the job index in a serial job)
fbname=$(basename "$FILENAME" .bam)

# count alleles using gatk
COMMAND="$GATK -T ASEReadCounter \
    -R $REF \
    -o $CSVPATH/$fbname.csv \
    -I $FILENAME \
    -sites $VCF/htlv-all-biallelic.vcf \
    -U ALLOW_N_CIGAR_READS \
    --minMappingQuality 10 \
    --minBaseQuality 2 \
    -minDepth 10"

eval $COMMAND