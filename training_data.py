from setup import *

!wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr21.fa.gz
!wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr22.fa.gz
!wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz
### Schema at https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.sql

#Extraction of the files
!gunzip -f chr21.fa.gz
!bgzip -l 9 chr21.fa
!samtools faidx chr21.fa.gz
!gunzip -f chr22.fa.gz
!bgzip -l 9 chr22.fa
!samtools faidx chr22.fa.gz