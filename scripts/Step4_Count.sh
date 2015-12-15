#!/bin/bash
#Step 4 of RNA-seq analysis: Counting reads to genes

annotation=/home/yasser/bio720/final_project/data/esal_reference/Esalsugineum_173_v1.0.gene_exons.gff3;
mkdir htseq_counts;

for bam_file in */*.bam ; do
    /usr/bin/htseq-count -f bam -r pos -s no -t exon -i Parent -q $bam_file $annotation > htseq_counts/`basename $bam_file .Aligned.sortedByCoord.out.bam`.htseq-count.txt &
    done;

