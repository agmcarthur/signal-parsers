#!/bin/sh

rm -rf /home/mcarthua/$1_*
mkdir /home/mcarthua/$1_fastq
cp bioproject_fastq/* /home/mcarthua/$1_fastq/.
tar cvf /home/mcarthua/$1_results.tar results results_mutations_freq.tsv results_mutations.tsv results.tsv pangolin summary*
gzip /home/mcarthua/$1_results.tar
chgrp mcarthua /home/mcarthua/$1_*
chown mcarthua /home/mcarthua/$1_*
