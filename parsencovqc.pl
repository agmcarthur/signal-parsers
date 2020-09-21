#!/usr/bin/perl -w

chomp(@files=<*.summary.qc.tsv>);

foreach $file (@files) {
	open (INPUT,"< $file");
	while (defined($line=<INPUT>)) {
		chomp($line);
		if ($line =~ /num_consensus_iupac/) {
			chomp($line=<INPUT>);
			@temp = split(/\t/,$line);
			$result{$temp[0]}=$temp[15];
			$stats{$temp[15]}++;	
		}
	}
	close (INPUT);
}

foreach $entry (sort keys %stats) {
	print STDERR "$stats{$entry}\t$entry\n";
}

foreach $entry (sort keys %result) {
	print "$entry\t$result{$entry}\n";
}

