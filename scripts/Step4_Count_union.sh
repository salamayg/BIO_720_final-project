#!/bin/bash
#Step 4 of RNA-seq analysis: Counting reads to genes using union mode.

annotation=/home/yasser/bio720/final_project/data/esal_reference/Esalsugineum_173_v1.0.gene_exons.gff3;
dir=htseq_counts_union_newest;
mkdir $dir;

for bam_file in */*.bam ; do
    touch $dir/`basename $bam_file .Aligned.sortedByCoord.out.bam`.htseq-count.log 
    /usr/bin/htseq-count -f bam -r pos -s no -t exon -i Parent $bam_file $annotation 2> $dir/`basename $bam_file .Aligned.sortedByCoord.out.bam`.htseq-count.log > $dir/`basename $bam_file .Aligned.sortedByCoord.out.bam`.htseq-count.txt &
    done;

