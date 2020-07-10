#!/bin/sh

tar cvf /home/mcarthua/$1_ncov.tar plots qc_analysis
gzip /home/mcarthua/$1_ncov.tar
chgrp mcarthua /home/mcarthua/$1*
chown mcarthua /home/mcarthua/$1*
