#!/usr/bin/perl -w

chomp(@directories=<*>);
system("mkdir bioproject_fastq");

foreach $directory (@directories) {
	unless ($directory =~ /results/ or $directory =~ /bioproject_fastq/) {
		system("cp $directory/host_removal/${directory}_R1.fastq.gz bioproject_fastq/${directory}_R1.fastq.gz");
		system("cp $directory/host_removal/${directory}_R2.fastq.gz bioproject_fastq/${directory}_R2.fastq.gz");
	}
}
