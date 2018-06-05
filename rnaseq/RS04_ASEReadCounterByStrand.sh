#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#@ example script for allele specific analysis using GATK. 

# aliases

# software
GATKPATH=</PATH/TO/GATK.jar>
GATK="java -jar $GATKPATH "

# files
OUTPATH=<PATH/TO/STRANDED/BAM>
REF=<PATH/TO/FASTA/REFERENCE.fa> # using a concatanated reference hg19 + HTLV-1
CSVPATH=PATH/TO/CSV/OUTPUT/>

# files to process (passed in from the job index in a serial job)
FILENAME=<PATH/TO/INPUT/FILE.bam>
fbname=$(basename "$FILENAME" .bam)


# split bam files based on mapping direction

# using samtools to identify mapped reads by flag (plus or minus strand)
samtools view -f 80 $FILENAME > $OUTPATH/$fbname-80.sam
samtools view -f 160 $FILENAME > $OUTPATH/$fbname-160.sam
samtools view -f 96 $FILENAME > $OUTPATH/$fbname-96.sam
samtools view -f 144 $FILENAME > $OUTPATH/$fbname-144.sam

# using samtools to extract sam header only
samtools view -H $FILENAME > $OUTPATH/$fbname-h.sam

# concatanate based on orientation and convert back to bam 
cat $OUTPATH/$fbname-h.sam $OUTPATH/$fbname-80.sam $OUTPATH/$fbname-160.sam | samtools view -bS - > $OUTPATH/$fbname-same.bam 
cat $OUTPATH/$fbname-h.sam $OUTPATH/$fbname-96.sam $OUTPATH/$fbname-144.sam | samtools view -bS - > $OUTPATH/$fbname-anti.bam 

# sort and index
samtools sort $OUTPATH/$fbname-same.bam $OUTPATH/$fbname-same.s 
samtools sort $OUTPATH/$fbname-anti.bam $OUTPATH/$fbname-anti.s

mv $OUTPATH/$fbname-same.s.bam $OUTPATH/$fbname-same.bam
mv $OUTPATH/$fbname-anti.s.bam $OUTPATH/$fbname-anti.bam

samtools index $OUTPATH/$fbname-same.bam 
samtools index $OUTPATH/$fbname-anti.bam 

rm $OUTPATH/$fbname-*.sam


# using GATK to count alleles separately for the expression from the plus / minus stand. 
COMMAND="$GATK -T ASEReadCounter \
    -R $REF \
    -o $CSVPATH/$fbname-same.csv \
    -I $OUTPATH/$fbname-same.bam \
    -sites ~/htlv/vcf/htlv-all-biallelic.vcf \
    -U ALLOW_N_CIGAR_READS \
    --minMappingQuality 10 \
    --minBaseQuality 2"

eval $COMMAND


COMMAND="$GATK -T ASEReadCounter \
    -R $REF \
    -o $CSVPATH/$fbname-anti.csv \
    -I $OUTPATH/$fbname-anti.bam \
    -sites ~/htlv/vcf/htlv-all-biallelic.vcf \
    -U ALLOW_N_CIGAR_READS \
    --minMappingQuality 10 \
    --minBaseQuality 2"

eval $COMMAND


