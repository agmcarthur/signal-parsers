#!/usr/bin/perl -w

unless (@ARGV) {
	print STDERR "\tFirst argument is tab-delimited versin of All Sequencing Results.xlsx\n";
	print STDERR "\tSecond argument is tab-delimited isolate metadata, columns as:\n";
	print STDERR "\t\tSample Identifier\n";
	print STDERR "\t\tDe-Identified Biosample Source ID\n";
	print STDERR "\t\tSubmit to NCBI/GSIAID (YES OR NO)\n";
	print STDERR "\t\tOriginal Patient Sample (e.g. from nasal swab)? YES OR NO\n";
	print STDERR "\t\tDate Sample Collected (date, month, or year)\n";
	print STDERR "\t\tHealthcare Facility that Collected the Sample\n";
	print STDERR "\t\t1st 3 or 4 Digits of Patient’s Postal Code \n";
	print STDERR "\t\tYear of Birth of Patient \n";
	print STDERR "\tSTDOUT is a annotated version of All Sequencing Results\n";
	exit;
}

unless (-e "$ARGV[0]") {
	print STDERR "\tCannot find: $ARGV[0]\n";
	exit;
}

unless (-e "$ARGV[1]") {
	print STDERR "\tCannot find: $ARGV[1]\n";
	exit;
}

open (INPUT,"< $ARGV[0]");
while (defined($line=<INPUT>)) {
	chomp($line);
	if ($line =~ /De-Identified/) {
		$header = $line;
		next;
	}
	@temp=split(/\t/,$line);
	$count = @temp;
	if (defined($data{$temp[4]})) {
		print STDERR "Repeat in $ARGV[0]: Skipping $temp[0] $temp[4]\n";
	}
	for ($i=0; $i <= ($count -1); $i++) {
		$data{$temp[4]}{$i}=$temp[$i];
	}
	$i = 13;	# QC: Genome Fraction less than 90%
	if ($data{$temp[4]}{$i} eq "FAIL") {
		$i=6;
		$data{$temp[4]}{$i}="FAIL";	# Submit to NCBI/GSIAID (YES OR NO)
	}
	$i = 18;	# QC: At least 90% of positions have better than 100x Coverage
	if ($data{$temp[4]}{$i} eq "FAIL") {
		$i=6;
		$data{$temp[4]}{$i}="FAIL";	# Submit to NCBI/GSIAID (YES OR NO)
	}
}
close (INPUT);

open (INPUT,"< $ARGV[1]");
while (defined($line=<INPUT>)) {
	chomp($line);
	if ($line =~ /De-Identified/) {
		next;
	}
	@temp = split(/\t/,$line);
	if (defined($data{$temp[0]})) {
		if (defined($temp[1])) {
			$i=5;
			$data{$temp[0]}{$i}=$temp[1];	# De-Identified Biosample Source ID
		}
		if (defined($temp[2])) {
			$i=6;
			$data{$temp[0]}{$i}=$temp[2];	# Submit to NCBI/GSIAID (YES OR NO)
		}
		if (defined($temp[3])) {
			$i=7;
			$data{$temp[0]}{$i}=$temp[3];	# Original Patient Sample (e.g. from nasal swab)? YES OR NO
		}
		if (defined($temp[4])) {
			$i=8;
			$data{$temp[0]}{$i}=$temp[4];	# Date Sample Collected (date, month, or year)
		}
		if (defined($temp[5])) {
			$i=9;
			$data{$temp[0]}{$i}=$temp[5];	# Healthcare Facility that Collected the Sample
		}
		if (defined($temp[6])) {
			$i=10;
			$data{$temp[0]}{$i}=$temp[6];	# 1st 3 or 4 Digits of Patient’s Postal Code
		}
		if (defined($temp[7])) {
			$i=11;
			$data{$temp[0]}{$i}=$temp[7];	# Year of Birth of Patient
		}
		$i = 13;	# QC: Genome Fraction less than 90%
		if ($data{$temp[0]}{$i} eq "FAIL") {
			$i=6;
			$data{$temp[0]}{$i}="FAIL";	# Submit to NCBI/GSIAID (YES OR NO)
		}
		$i = 18;	# QC: At least 90% of positions have better than 100x Coverage
		if ($data{$temp[0]}{$i} eq "FAIL") {
			$i=6;
			$data{$temp[0]}{$i}="FAIL";	# Submit to NCBI/GSIAID (YES OR NO)
		}
	}
}
close (INPUT);

print "$header\n";
foreach $entry (sort keys %data) {
	$count = keys %{$data{$entry}};
	for ($i=0; $i <= ($count -1); $i++) {
		print "$data{$entry}{$i}\t";
	}
	print "\n";
}
