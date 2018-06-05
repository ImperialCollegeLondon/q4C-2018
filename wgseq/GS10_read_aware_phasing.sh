#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018


# based on steps here https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#readaware
# haplotype reference data from https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html (already in NCBI37 - hg19 so should correspond to my data)


# choose chromosome - chromosomes passed in from job index in a serial job
CHROM=<CHROM_NUM>


# aliases

# software paths
BINPATH=<SOFTWARE/PATH/>
EXPIR=$BINPATH"/extractPIRs"
SHAPEIT=$BINPATH"shapeit2/bin"

# reference paths
REFDIR=<PATH/TO/1000GP/REFERENCE/PATH/>
REFDATA=$REFDIR"/1000GP_Phase3_chr"$CHROM

# file paths
OUTDIR=<PATH/TO/PIR/OUTPUT/>
BAMPATH=<PATH/TO/BAM/FILES/>
BAMLIST="$OUTDIR/bamlist_"$CHROM
VCF=<PATH/TO/VCF/FILES/>
VCFDATA=$VCF"bychr/htlv-all-biallelic-chr"$CHROM".vcf"
OUTFILE="$OUTDIR/PIRlist_chr$CHROM"


# create bam list

touch $BAMLIST
for f in $BAMDIR/*.bam
do
	fbname=$(basename "$f" .bam)
	echo -e "$fbname\t$f\tchr$CHROM" >> $BAMLIST
done


# extract PIRs

COMMAND="$EXPIR/extractPIRs --bam $BAMLIST \
            --base-quality 20 \
            --vcf $VCFDATA \
            --out $OUTFILE"

eval $COMMAND


# check alignment of reference to variants

COMMAND="$SHAPEIT/shapeit -check \
        --input-vcf $VCFDATA \
        -R $REFDATA.hap.gz $REFDATA.legend.gz $REFDIR/1000GP_Phase3.sample \
        --output-log $OUTDIR/myAlignmentChecks_chr$CHROM.log"

eval $COMMAND


# assemble haplotypes using shapeit. 

COMMAND="$SHAPEIT/shapeit -assemble \
        --input-vcf $VCFDATA \
        --input-pir $OUTFILE \
        -R $REFDATA.hap.gz $REFDATA.legend.gz $REFDIR/1000GP_Phase3.sample \
        -O $OUTDIR/HTLV_chr$CHROM \
        --exclude-snp $OUTDIR/myAlignmentChecks_chr$CHROM.snp.strand.exclude \
        --no-mcmc"

eval $COMMAND


