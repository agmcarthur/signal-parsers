#!/usr/bin/perl -w

chomp(@directories=<*>);
system("mkdir results");

foreach $directory (@directories) {
	system("cp $directory/${directory}* results/.");
	open (INPUT,"< $directory/${directory}_statistics.tsv");
	chomp($header=<INPUT>);
	chomp($data=<INPUT>);
	$results{$directory}=$data;
	close (INPUT);
}
open (OUTPUT,"> results.tsv");
print OUTPUT "Sample\t$header\n";
foreach $entry (sort keys %results) {
	print OUTPUT "$entry\t$results{$entry}\n";
}
close (OUTPUT);
