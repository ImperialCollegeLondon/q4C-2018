#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#@ example script for coverage analysis using bedtools. 
#@ the same procedure should be done using multiple sam flags, 
#@ to account for the different combination of R1/R2 plus sense / minus sense.
#@ here example shown is using 80 as flag. use also 160, 144, 96.  

# aliases

# software
BEDTOOLS=<PATH/TO/BEDTOOLS/BIN>

# files
OUTPATH=<PATH/TO/SORTED/BAM/FILE/INTERMEDIATE/>
WINDOWS=<PATH/TO/windows.bed> # 1 kb window from 5Mb upstream to 5Mb downstream to each integration site. bed format

# files to process (passed in from the job index in a serial job)
FILENAME=<PATH/TO/INPUT/FILE.bam>
fbname=$(basename "$FILENAME" .bam)

# run bedtools coverage on file
COMMAND="samtools view -f 80 -b $FILENAME | $BEDTOOLS coverage -b - -a $WINDOWS > $OUTPATH/$fbname.bed"
eval $COMMAND

