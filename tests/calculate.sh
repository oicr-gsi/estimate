#!/bin/bash
cd $1

# - .txt files with information produced with RSEM and STAR
# - .gct file made with estimate

# Therefore:
# - Check md5sums for all types of files, sort

echo ".txt files:"
for v in *.txt;do echo $v;md5sum $v;done 

echo ".gct files:"
find . -name "*.gct" | xargs md5sum | sort

