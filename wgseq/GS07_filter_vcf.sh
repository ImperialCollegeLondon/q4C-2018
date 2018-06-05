#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018


# filter VCF to only biallelic variants to enable allele counting
############################################################################

# aliases

# files
VCF=<PATH/TO/VCF/>

# software
GATKPATH=<PATH/TO/GenomeAnalysisTK.jar>
GATK="java -jar $GATKPATH "

# reference
REF=<PATH/TO/FASTA/REFERENCE.fa> # using a concatanated reference hg19 + HTLV-1

# filter vcf using gatk
COMMAND="$GATK -T SelectVariants \
   -R $REF \
   -V $VCF/htlv-all.vcf \
   -o $VCF/htlv-all-biallelic.vcf \
   -restrictAllelesTo BIALLELIC"

eval $COMMAND




