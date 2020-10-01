#!/bin/sh

rm -rf /home/mcarthua/$1_ncov.tar.gz
tar cvf /home/mcarthua/$1_ncov.tar *
gzip /home/mcarthua/$1_ncov.tar
chgrp mcarthua /home/mcarthua/$1*
chown mcarthua /home/mcarthua/$1*
