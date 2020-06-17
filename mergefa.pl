#!/usr/bin/perl -w

chomp(@files=<*.fa>);
foreach $file (@files) {
	($name) = $file =~ /^(.*)\.fa/;
	open (INPUT,"< $file");
	chomp($line=<INPUT>);
	chomp($line=<INPUT>);
	print ">HCoV-19/Canada/Toronto/2020/$name\n";
	print "$line\n";
	close (INPUT);
}
