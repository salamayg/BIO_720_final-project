#!/bin/tcsh 
#This script will get the number of uniquely mapped reads for each sample and output a csv with sample name and the number of reads. Useful for data visualization later

touch mapped_statistics.txt

foreach x ( /home/yasser/bio720/final_project/data/processed_data/*/*Log.final* )
   echo "`basename $x`," >> mapped_statistics.txt
   grep Uniquely $x | cut -f 2 | grep -v "%" >> mapped_statistics.txt
end

perl -pi -e 's/,\n/,/g' mapped_statistics.txt
perl -pi -e 's/Sample_//g' mapped_statistics.txt
perl -pi -e 's/.Log.final.out//g' mapped_statistics.txt
