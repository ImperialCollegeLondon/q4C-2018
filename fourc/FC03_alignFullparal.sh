#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#@ FC03_alignFullparal - alignment of single read data (R1 and R2 aligned separately) against
#                       mixed reference of human + virus. 

# input variables:  $outputdir - <path/to/output/> 
#                   $filename - <file/path/to/align>.fq.gz

# require adding samtools to path (e.g. version 0.1.19)

# requires creating an index from concatanated hg19, HTLV-1 genomes (AB513134 used here)
INDEX_PATH=<PATH/TO/BOWTIE/INDEX>

date
echo " "

echo $outputdir
echo $filename

## copy file to location
cp $filename .
gunzip *.gz

ls -l -h -R

#full reference ONLY, BOWTIE2 ONLY
for mate in *.fq
do
    bamName=`echo $mate  | sed 's#.fq#_full.bam#g' - ` 
    samName=`echo $mate  | sed 's#.fq#_full.sam#g' - ` 
    echo "mate is "$mate   
    echo "output is "$bamName

    if [ ! -f $outputdir$bamName".gz" ]; then
        echo "using bowtie2"
        bowtie2 \
        -x $INDEX_PATH \
        -k 100 \
        -L 20 \
        -p 20 \
        --reorder \
        -U $mate \
        | samtools view -bS - > $bamName
        
        samtools view -F 256 $bamName > $samName
        gzip $samName
        gzip $bamName
        
    else
        echo $bamName" already exists, skipping."
    fi
done


date
ls -l -h -R

date
cp *.bam.gz $outputdir
cp *.sam.gz $outputdir

echo "   done!"