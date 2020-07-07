#!/usr/bin/perl -w

open (OUTPUT1,"> results_mutations.tsv");

open (INPUT,"< results.tsv");
while (defined($line=<INPUT>)) {
	chomp($line);
	if ($line =~ /Nucleotides/) {
		next;
	}
	@temp = split(/\t/,$line);
	if ($temp[40] eq "n.d.") {
		next;
	}
	$fraction{$temp[0]} = $temp[23];
	$coverage{$temp[0]} = $temp[30];
	$temp[40] =~ s/^ //;
	$temp[40] =~ s/,//g;
	$temp[40] =~ s/&nbsp;/ /g;
	$temp[40] =~ s/&#8209;/-/g;
	$temp[40] =~ s/&#8231;/-/g;
	$temp[40] =~ s/&#8211;/-/g;
	$temp[40] =~ s/&Delta;/delta/g;
	@mutations = split(/;/,$temp[40]);
	foreach $mutant (@mutations) {
		($position,$description,$gene) = $mutant =~ /^[A-Z]+ ([0-9]+)(.*of reads\)) (.*)/;
		if ($mutant =~ /^MC INDEL/) {
			($position,$gene) = $mutant =~ /^MC INDEL ([0-9]+).*of reads\) (.*)/;
			($description) = $mutant =~ /(INDEL.*of reads\))/;
		}
		if (defined($position)) {
			$annotation{$position}=$gene;
			$data{$position}{$temp[0]}=$description;
			($desc2) = $description =~ /^(.*) \([0-9]+\.[0-9]+\% of reads\)/;
			unless ($desc2) {
				($desc2) = $description =~ /^(.*) \([0-9]+\% of reads\)/;
			}
			$freqdesc = "$position\t$gene\t$desc2";
			$frequency{$freqdesc}++;
		} else {
			print STDERR "\tSkipping $temp[0]: $mutant\n";	
		}
	}
}
close (INPUT);

print OUTPUT1 "\t";
foreach $sample (sort keys %fraction) {
	print OUTPUT1 "\t$sample";
}
print OUTPUT1 "\n";

print OUTPUT1 "MN908947.3\tGenome Fraction (%)";
foreach $sample (sort keys %fraction) {
	print OUTPUT1 "\t$fraction{$sample}";
}
print OUTPUT1 "\n";

print OUTPUT1 "Genome Position\tAverageDepth of Coverage (fold)";
foreach $sample (sort keys %fraction) {
	print OUTPUT1 "\t$coverage{$sample}";
}
print OUTPUT1 "\n";

foreach $location (sort {$a <=> $b} keys %data) {
	print OUTPUT1 "$location\t$annotation{$location}";
	foreach $sample (sort keys %fraction) {
		if (defined($data{$location}{$sample})) {
			print OUTPUT1 "\t$data{$location}{$sample}";
		} else {
			print OUTPUT1 "\t";
		}
	}
	print OUTPUT1 "\n";
}

close (OUTPUT1);

open (OUTPUT2,"> results_mutations_freq.tsv");
print OUTPUT2 "# Samples\tPosition\tGenee\tMutation\n";
foreach $entry (sort keys %frequency) {
	print OUTPUT2 "$frequency{$entry}\t$entry\n";
}
close (OUTPUT2);
