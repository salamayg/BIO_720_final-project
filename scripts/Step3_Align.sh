#!/bin/bash
#Step 3 of RNA-seq analysis: Read mapping to reference

#The reads for each sample are stored in separate directories. This script will go through each directory and map reads against the reference built by STAR using STAR. 
#If the sample was split over two lanes, it will have four files, so the if statements will differentiate these two scenarios and output accordingly.
#Notable flags in the alignment command are: --readFilesCommand zcat to read gzipped files as input, and --outSAMtype BAM SortedByCoordinate which outputs the files as BAM and already sorted.
for sample in Sample_* ; do
    cd $sample;
    FILES=(*.gz);

    if [ ${#FILES[@]} = 2 ]; then
	R1=${FILES[0]}
	R2=${FILES[1]}
	/usr/local/STAR_2.4.2/bin/Linux_x86_64_static/STAR --runMode alignReads --runThreadN 14 --genomeDir /home/yasser/bio720/final_project/data/esal_reference/STAR_reference --readFilesIn $R1 $R2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix ./$sample.

    elif [ ${#FILES[@]} = 4 ]; then
	A_R1=${FILES[0]}
	A_R2=${FILES[1]}
	B_R1=${FILES[2]}
	B_R2=${FILES[3]}
	/usr/local/STAR_2.4.2/bin/Linux_x86_64_static/STAR --runMode alignReads --runThreadN 6 --genomeDir /home/yasser/bio720/final_project/data/esal_reference/STAR_reference --readFilesIn $A_R1 $A_R2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix ./"$sample.L1." &
	/usr/local/STAR_2.4.2/bin/Linux_x86_64_static/STAR --runMode alignReads --runThreadN 6 --genomeDir /home/yasser/bio720/final_project/data/esal_reference/STAR_reference --readFilesIn $B_R1 $B_R2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix ./"$sample.L2."  
    fi
    cd ../;
done;
