#!/usr/bin/perl -w

# run in the qc_analysis folder, lists samples with ambiguities found in more than one sample (possible contamination)
# treat this the same was as parsencovqc.pl

chomp(@files=<*alleles.tsv>);

foreach $file (@files) {
	open (INPUT,"< $file");
	while (defined($line=<INPUT>)) {
		chomp($line);
		if ($line =~ /samples_with_allele/) {
			next;
		}
		@temp = split(/\t/,$line);
		unless ($temp[3] =~ /[ACGTN]/) {
			if ($temp[4] > 2) {
				($name) = $temp[0] =~ /Consensus_(.*)\.consensus/;
				print "$name\tFREQ_AMBIGUITY\n";
			}
		}
	}
	close (INPUT);
}

