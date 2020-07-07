#!/bin/sh

tar cvf /home/mcarthua/$1.tar bioproject_fastq results results_mutations_freq.tsv results_mutations.tsv results.tsv summary*
gzip /home/mcarthua/$1.tar
chgrp mcarthua /home/mcarthua/$1.tar.gz
chown mcarthua /home/mcarthua/$1.tar.gz
