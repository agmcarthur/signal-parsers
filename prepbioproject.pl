#!/usr/bin/perl -w

chomp(@directories=<*>);
system("mkdir bioproject_fastq");

foreach $directory (@directories) {
	unless ($directory =~ /results/ or $directory =~ /bioproject_fastq/) {
		system("cp $directory/host_removal/R1.fastq.gz bioproject_fastq/${directory}-R1.fastq.gz");
		system("cp $directory/host_removal/R2.fastq.gz bioproject_fastq/${directory}-R2.fastq.gz");
	}
}
