#!/bin/sh

tar cvf /home/mcarthua/$1_fastq.tar bioproject_fastq
tar cvf /home/mcarthua/$1_results.tar results results_mutations_freq.tsv results_mutations.tsv results.tsv summary*
gzip /home/mcarthua/$1.tar
chgrp mcarthua /home/mcarthua/$1*
chown mcarthua /home/mcarthua/$1*
