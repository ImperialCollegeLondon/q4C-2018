#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018




# aliases

GVCF=<PATH/TO/GVCF/>
VCF=<PATH/TO/VCF/OUTPUT/>
GATKPATH=<PATH/TO/GenomeAnalysisTK.jar>

GATK="java -jar $GATKPATH "
REF=<PATH/TO/FASTA/REFERENCE.fa> # using a concatanated reference hg19 + HTLV-1


# run genotypeGVCFs to generate a merged list of variants across samples 

COMMAND="$GATK -T GenotypeGVCFs \
   -R $REF \
   -V $GVCF/gvcf.list \
   -nt 10 \
   -o $VCF/htlv-all-T.vcf"

eval $COMMAND

