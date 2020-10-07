#!/usr/bin/perl -w

print STDERR "Scatterplots: https://scatterplot.online\n";

system("rm -rf results_plots");
system("mkdir results_plots");

$fractionpass = 0;
$pass100x = 0;
$pass1000x = 0;
$passindel = 0;
$passframeshift = 0;

if ($ARGV[0]) {							# load Ct values
	open (INPUT,"< $ARGV[0]");
	while (defined($line=<INPUT>)) {
		chomp($line);
		@temp = split(/\t/,$line);
		if (defined($temp[1])) {
			if ($temp[1] =~ /[0-9]/) {
				$ctvalue{$temp[0]}=$temp[1];
			}
		}
	}
	close (INPUT);
}

open (INPUT,"< results.tsv");			# load data
while (defined($line=<INPUT>)) {
	chomp($line);
	if ($line =~ /Pipeline Version/) {
		next;
	}
	@temp = split(/\t/,$line);
	$genomefraction{$temp[0]}=$temp[23];
	$sarscontent{$temp[0]}=$temp[20];
	$avecoverage{$temp[0]}=$temp[30];

	$ambiquity{$temp[0]}=$temp[25];
	$mismatch{$temp[0]}=$temp[26];
	$indels{$temp[0]}=$temp[27];

	$value = $temp[34] + $temp[35] + $temp[36] + $temp[37];
	$coverage100plus{$temp[0]} = $value;

	$count++;
	if ($temp[2] eq "PASS") {
		$fractionpass++;
	}
	if ($temp[7] eq "PASS") {
		$pass100x++;
	}
	if ($temp[8] eq "PASS") {
		$pass1000x++;
	}
	if ($temp[9] eq "PASS") {
		$passindel++;
	}
	if ($temp[10] eq "PASS") {
		$passframeshift++;
	}
	
}
close (INPUT);

print STDERR "$fractionpass of $count PASS Genome Fraction of at least 90\%\n";
print STDERR "$pass100x of $count PASS 90\% of positions have at least 100x coverage\n";
print STDERR "$pass1000x of $count PASS 90\% of positions have at least 1000x coverage\n";
print STDERR "$passindel of $count PASS No indels detected (maximum length 85bp)\n";
print STDERR "$passframeshift of $count PASS No frameshifts in SARS-CoV-2 open reading frames\n";

# overview plot

open (OUTPUT,"> results_plots/overview.tsv");
print OUTPUT "Sample\t\%SARS-CoV-2\tGenome Fraction (\%)\tColour\n";
foreach $sample (keys %genomefraction) {
	$label = "$sample ${avecoverage{$sample}}x";
	print OUTPUT "$label\t$sarscontent{$sample}\t$genomefraction{$sample}";
	if ($coverage100plus{$sample} > .90) {
		if ($ctvalue{$sample}) {
			if ($ctvalue{$sample} > 30) {
				print OUTPUT "\t>90% of positions have at least 100x coverage & Ct > 30\n";
			} else {
				print OUTPUT "\t>90% of positions have at least 100x coverage & Ct <= 30\n";
			}		
		} else {
			print OUTPUT "\t>90% of positions have at least 100x coverage (no Ct data)\n";
		}
	} else {
		if ($ctvalue{$sample}) {
			if ($ctvalue{$sample} > 30) {
				print OUTPUT "\t<90% of positions have at least 100x coverage & Ct > 30\n";
			} else {
				print OUTPUT "\t<90% of positions have at least 100x coverage & Ct <= 30\n";
			}		
		} else {
			print OUTPUT "\t<90% of positions have at least 100x coverage (no Ct data)\n";
		}
	}
}
close (OUTPUT);

# Ct plot

if (%ctvalue) {

	open (OUTPUT,"> results_plots/ctplot.tsv");
	print OUTPUT "Sample\tCt Value\tGenome Fraction (\%)\tColour\n";
	foreach $sample (keys %genomefraction) {
		$label = "$sample ${avecoverage{$sample}}x";
		if (defined($ctvalue{$sample})) {
			print OUTPUT "$label\t$ctvalue{$sample}\t$genomefraction{$sample}";
			if ($coverage100plus{$sample} > .90) {
				print OUTPUT "\t>90% of positions have at least 100x coverage\n";
			} else {
				print OUTPUT "\t<90% of positions have at least 100x coverage\n";
			}
		}
	}
	close (OUTPUT);	

	open (OUTPUT,"> results_plots/ct_ambiquity.tsv");
	print OUTPUT "Sample\tCt Value\tN's per 100 kbp\tColour\n";
	foreach $sample (keys %genomefraction) {
		$label = "$sample ${avecoverage{$sample}}x";
		if (defined($ctvalue{$sample})) {
			print OUTPUT "$label\t$ctvalue{$sample}\t$ambiquity{$sample}";
			if ($coverage100plus{$sample} > .90) {
				print OUTPUT "\t>90% of positions have at least 100x coverage\n";
			} else {
				print OUTPUT "\t<90% of positions have at least 100x coverage\n";
			}
		}
	}
	close (OUTPUT);

	open (OUTPUT,"> results_plots/ct_mismatch.tsv");
	print OUTPUT "Sample\tCt Value\tMismatches per 100 kbp\tColour\n";
	foreach $sample (keys %genomefraction) {
		$label = "$sample ${avecoverage{$sample}}x";
		if (defined($ctvalue{$sample})) {
			print OUTPUT "$label\t$ctvalue{$sample}\t$mismatch{$sample}";
			if ($coverage100plus{$sample} > .90) {
				print OUTPUT "\t>90% of positions have at least 100x coverage\n";
			} else {
				print OUTPUT "\t<90% of positions have at least 100x coverage\n";
			}
		}
	}
	close (OUTPUT);

	open (OUTPUT,"> results_plots/ct_indels.tsv");
	print OUTPUT "Sample\tCt Value\tIndels per 100 kbp\tColour\n";
	foreach $sample (keys %genomefraction) {
		$label = "$sample ${avecoverage{$sample}}x";
		if (defined($ctvalue{$sample})) {
			print OUTPUT "$label\t$ctvalue{$sample}\t$indels{$sample}";
			if ($coverage100plus{$sample} > .90) {
				print OUTPUT "\t>90% of positions have at least 100x coverage\n";
			} else {
				print OUTPUT "\t<90% of positions have at least 100x coverage\n";
			}
		}
	}
	close (OUTPUT);

}

# coverage plots

open (OUTPUT,"> results_plots/ave_coverage.tsv");
print OUTPUT "Sample\tAverage Depth of Coverage\tGenome Fraction (\%)\tColour\n";
foreach $sample (keys %genomefraction) {
	print OUTPUT "$sample\t$avecoverage{$sample}\t$genomefraction{$sample}";
	if ($ctvalue{$sample}) {
		if ($ctvalue{$sample} > 30) {
			print OUTPUT "\tCt > 30\n";
		} else {
			print OUTPUT "\tCt <= 30\n";
		}		
	} else {
		print OUTPUT "\tNo Ct Data)\n";
	}
}
close (OUTPUT);

open (OUTPUT,"> results_plots/100x_coverage.tsv");
print OUTPUT "Sample\tFraction of Genome with >100x coverage\tGenome Fraction (\%)\tColour\n";
foreach $sample (keys %genomefraction) {
	print OUTPUT "$sample\t$coverage100plus{$sample}\t$genomefraction{$sample}";
	if ($ctvalue{$sample}) {
		if ($ctvalue{$sample} > 30) {
			print OUTPUT "\tCt > 30\n";
		} else {
			print OUTPUT "\tCt <= 30\n";
		}		
	} else {
		print OUTPUT "\tNo Ct Data)\n";
	}
}
close (OUTPUT);

open (OUTPUT,"> results_plots/100x_ambiquity.tsv");
print OUTPUT "Sample\tFraction of Genome with >100x coverage\tN's per 100 kbp\tColour\n";
foreach $sample (keys %genomefraction) {
	$label = "$sample ${avecoverage{$sample}}x";
	print OUTPUT "$label\t$coverage100plus{$sample}\t$ambiquity{$sample}";
	if ($ctvalue{$sample}) {
		if ($ctvalue{$sample} > 30) {
			print OUTPUT "\tCt > 30\n";
		} else {
			print OUTPUT "\tCt <= 30\n";
		}		
	} else {
		print OUTPUT "\tNo Ct Data)\n";
	}
}
close (OUTPUT);

open (OUTPUT,"> results_plots/100x_mismatch.tsv");
print OUTPUT "Sample\tFraction of Genome with >100x coverage\tMismatches per 100 kbp\tColour\n";
foreach $sample (keys %genomefraction) {
	$label = "$sample ${avecoverage{$sample}}x";
	print OUTPUT "$label\t$coverage100plus{$sample}\t$mismatch{$sample}";
	if ($ctvalue{$sample}) {
		if ($ctvalue{$sample} > 30) {
			print OUTPUT "\tCt > 30\n";
		} else {
			print OUTPUT "\tCt <= 30\n";
		}		
	} else {
		print OUTPUT "\tNo Ct Data)\n";
	}
}
close (OUTPUT);

open (OUTPUT,"> results_plots/100x_indels.tsv");
print OUTPUT "Sample\tFraction of Genome with >100x coverage\tIndels per 100 kbp\tColour\n";
foreach $sample (keys %genomefraction) {
	$label = "$sample ${avecoverage{$sample}}x";
	print OUTPUT "$label\t$coverage100plus{$sample}\t$indels{$sample}";
	if ($ctvalue{$sample}) {
		if ($ctvalue{$sample} > 30) {
			print OUTPUT "\tCt > 30\n";
		} else {
			print OUTPUT "\tCt <= 30\n";
		}		
	} else {
		print OUTPUT "\tNo Ct Data)\n";
	}
}
close (OUTPUT);

# seq error

open (OUTPUT,"> results_plots/ambiquity.tsv");
print OUTPUT "Sample\tN's per 100 kbp\tGenome Fraction (\%)\tColour\n";
foreach $sample (keys %genomefraction) {
	$label = "$sample ${avecoverage{$sample}}x";
	print OUTPUT "$label\t$ambiquity{$sample}\t$genomefraction{$sample}";
	if ($coverage100plus{$sample} > .90) {
		if ($ctvalue{$sample}) {
			if ($ctvalue{$sample} > 30) {
				print OUTPUT "\t>90% of positions have at least 100x coverage & Ct > 30\n";
			} else {
				print OUTPUT "\t>90% of positions have at least 100x coverage & Ct <= 30\n";
			}		
		} else {
			print OUTPUT "\t>90% of positions have at least 100x coverage (no Ct data)\n";
		}
	} else {
		if ($ctvalue{$sample}) {
			if ($ctvalue{$sample} > 30) {
				print OUTPUT "\t<90% of positions have at least 100x coverage & Ct > 30\n";
			} else {
				print OUTPUT "\t<90% of positions have at least 100x coverage & Ct <= 30\n";
			}		
		} else {
			print OUTPUT "\t<90% of positions have at least 100x coverage (no Ct data)\n";
		}
	}
}
close (OUTPUT);

open (OUTPUT,"> results_plots/mismatch.tsv");
print OUTPUT "Sample\tMismatches per 100 kbp\tGenome Fraction (\%)\tColour\n";
foreach $sample (keys %genomefraction) {
	$label = "$sample ${avecoverage{$sample}}x";
	print OUTPUT "$label\t$mismatch{$sample}\t$genomefraction{$sample}";
	if ($coverage100plus{$sample} > .90) {
		if ($ctvalue{$sample}) {
			if ($ctvalue{$sample} > 30) {
				print OUTPUT "\t>90% of positions have at least 100x coverage & Ct > 30\n";
			} else {
				print OUTPUT "\t>90% of positions have at least 100x coverage & Ct <= 30\n";
			}		
		} else {
			print OUTPUT "\t>90% of positions have at least 100x coverage (no Ct data)\n";
		}
	} else {
		if ($ctvalue{$sample}) {
			if ($ctvalue{$sample} > 30) {
				print OUTPUT "\t<90% of positions have at least 100x coverage & Ct > 30\n";
			} else {
				print OUTPUT "\t<90% of positions have at least 100x coverage & Ct <= 30\n";
			}		
		} else {
			print OUTPUT "\t<90% of positions have at least 100x coverage (no Ct data)\n";
		}
	}
}
close (OUTPUT);

open (OUTPUT,"> results_plots/indels.tsv");
print OUTPUT "Sample\tIndels per 100 kbp\tGenome Fraction (\%)\tColour\n";
foreach $sample (keys %genomefraction) {
	$label = "$sample ${avecoverage{$sample}}x";
	print OUTPUT "$label\t$indels{$sample}\t$genomefraction{$sample}";
	if ($coverage100plus{$sample} > .90) {
		if ($ctvalue{$sample}) {
			if ($ctvalue{$sample} > 30) {
				print OUTPUT "\t>90% of positions have at least 100x coverage & Ct > 30\n";
			} else {
				print OUTPUT "\t>90% of positions have at least 100x coverage & Ct <= 30\n";
			}		
		} else {
			print OUTPUT "\t>90% of positions have at least 100x coverage (no Ct data)\n";
		}
	} else {
		if ($ctvalue{$sample}) {
			if ($ctvalue{$sample} > 30) {
				print OUTPUT "\t<90% of positions have at least 100x coverage & Ct > 30\n";
			} else {
				print OUTPUT "\t<90% of positions have at least 100x coverage & Ct <= 30\n";
			}		
		} else {
			print OUTPUT "\t<90% of positions have at least 100x coverage (no Ct data)\n";
		}
	}
}
close (OUTPUT);


