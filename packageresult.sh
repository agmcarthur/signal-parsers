#!/bin/sh

tar cvf /home/mcarthua/$1.tar bioproject_fastq results results_mutations_freq.tsv results_mutations.tsv results.tsv   
gzip /home/mcarthua/$1.tar
