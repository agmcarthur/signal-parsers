#!/bin/sh

echo "assembly_report.pl"
rm -rf result*
rm -rf bioproject_fastq
for file in `ls | grep -v summary | grep -v config`; do cd $file; assembly_report.pl; cd ..; done
echo "assembly_summarize.pl"
assembly_summarize.pl
echo "assembly_mutations.pl"
assembly_mutations.pl
echo "prepbioproject.pl"
prepbioproject.pl
