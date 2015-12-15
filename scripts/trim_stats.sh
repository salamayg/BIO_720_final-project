#!/bin/tcsh 
#This script will read the BBDuk trimming output and give values for number of reads prior and post trimming. Useful for data visualization later.

touch pretrim_statistics.txt
touch posttrim_statistics.txt

foreach x (/home/yasser/bio720/final_project/data/processed_data/*/*bbduk*txt)
     grep Input: $x | cut -f 2 >> pretrim_statistics.txt
     grep Result: $x | cut -f 2 >> posttrim_statistics.txt
end

perl -pi -e 's/ reads//' pretrim_statistics.txt
perl -pi -e 's/ reads .*//' posttrim_statistics.txt
