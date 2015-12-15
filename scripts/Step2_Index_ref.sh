#!/usr/bin/perl 
#use STAR to generate genome reference index

#pathways
my $out_dir="/home/yasser/bio720/final_project/data/esal_reference/STAR_reference";
my $genome_file="/home/yasser/bio720/final_project/data/esal_reference/Esals173.fa";
my $gff3_file="/home/yasser/bio720/final_project/data/esal_reference/Esalsugineum_173_v1.0.gene_exons.gff3";

#Command to generate genome for STAR. --sjdbGTFfile specifies the gff3 annotation file. --sjdbGTFfeatureExon specifies what feature in the gff3 file to care about. --sjdbGTFtagExonParentTranscript Parent specifies what field to look for when grouping features from the same transcript/gene
exec " /usr/local/STAR_2.4.2/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --genomeDir $out_dir --genomeFastaFiles $genome_file --sjdbGTFfile $gff3_file --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --runThreadN 4";



