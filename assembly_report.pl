#!/usr/bin/perl -w

# filename

use Cwd;
$directory = getcwd;
@temp = split(/\//,$directory);
$filebase = $temp[-1];
print STDERR "$directory\n";

# defaults

$rawreadpairs = 0;
$rawbp = 0;
$R1withprimer = 0;
$R2withprimer = 0;
$postprimerreadpairs = 0;
$postprimerbp = 0;
$finalR1 = 0;
$length = 0;
$fraction = 0;
$features = 0;
$unresolved = 0;
$mismatches = 0;
$indels = 0;
$prime5 = 0;
$prime3 = 0;
$average = 0;
$mutations = " ";
$fold1++;
$fold10++;
$fold100++;
$fold1000++;
$fold2000++;
$fold10000++;
$bigfold++;

open (OUTPUT1,"> ${filebase}_statistics.tsv");

# header

print OUTPUT1 "Pipeline Version";
print OUTPUT1 "\tQC: Genome Fraction less than 90%\tQC: Depth of Coverage less than 2000 fold\tQC: Variant with less than 90% frequency among reads\tQC: Per base sequence quality\tQC: Adapter Content";
print OUTPUT1 "\tQC: At least 90% of positions have better than 100x Coverage\tQC: At least 90% of positions have better than 1000x Coverage";

print OUTPUT1 "\tQC: No indels detected (maximum length 85bp)\tQC: Frameshifts in SARS-CoV-2 open reading frames";

print OUTPUT1 "\tRaw Read Pairs\tRaw Nucleotides (bp)\tR1 with PCR primer\tR2 with PCR primer\tPost Primer Trim Read Pairs\tPost Primer Trim Nucelotides (bp)";
print OUTPUT1 "\t# Paired Sequences Used\tFASTQC FAIL\tFASTQC WARN";
print OUTPUT1 "\tSARS-COV-2 in Trimmed FASTQ (\%)";
print OUTPUT1 "\tSARS-COV-2 in Host Removed FASTQ (\%)";
print OUTPUT1 "\tGenome Length (bp)\tGenome Fraction (\%)\tGenomic Features\tN's per 100 kbp\tMismatches per 100 kbp\tindels per 100 kbp";
print OUTPUT1 "\t5' Ns\t3' Ns";
print OUTPUT1 "\tAverage Depth of Coverage\tFraction with 0 Coverage\tFraction 1x-10x Coverage\tFraction 11x-100x Coverage\tFraction 101x-1000x Coverage\tFraction 1000x-2000x Coverage\tFraction 2000x-10000x Coverage\tFraction with greater than 10000x Coverage";
print OUTPUT1 "\tTaxonomic Composition of Assembly";
print OUTPUT1 "\tiVar Variants";
print OUTPUT1 "\tBreSeq Variants";
print OUTPUT1 "\n";

open (INPUT,"< ${filebase}_sample.txt");
while (defined($line=<INPUT>)) {
	chomp($line);
	if ($line =~ /version/) {
		($version) = $line =~ /version (.*)/;
	}
	if ($line =~ /Raw Data \(read pairs\)/) {
		($rawreadpairs) = $line =~ /Raw Data \(read pairs\): ([0-9]+)/;
	}
	if ($line =~ /Raw Data \(base pairs\)/) {
		($rawbp) = $line =~ /Raw Data \(base pairs\): ([0-9]+)/;
	}
	if ($line =~ /Primer Removal \(read pairs\)/) {
		($postprimerreadpairs) = $line =~ /Primer Removal \(read pairs\): ([0-9]+)/;
	}
	if ($line =~ /Primer Removal \(base pairs\)/) {
		($postprimerbp) = $line =~ /Primer Removal \(base pairs\): ([0-9]+)/;
	}
	if ($line =~ /Post Trim \(read pairs\)/) {
		($finalR1) = $line =~ /Post Trim \(read pairs\): ([0-9]+)/;
	}
	if ($line =~ /FASTQC Flags/) {
		while ($line =~ /[a-zA-Z]/) {
			chomp($line=<INPUT>);
			if ($line =~ /FAIL/) {
				($entry) = $line =~ /FAIL  (.*)/;
				$fastqcfail .= "$entry; ";
			}
			if ($line =~ /WARN/) {
				($entry) = $line =~ /WARN  (.*)/;
				$fastqcwarn .= "$entry; ";
			}
		}
	}
	if ($line =~ /Reads SARS-CoV-2 \(\%\): /) {
		($sarsload) = $line =~ /Reads SARS-CoV-2 \(\%\): (.*)/;
	}
	if ($line =~ /Genome Length \(bp\): /) {
		($length) = $line =~ /Genome Length \(bp\): (.*)/;
	}
	if ($line =~ /Genome Fraction \(%\): /) {
		($fraction) = $line =~ /Genome Fraction \(%\): (.*)/;
	}
	if ($line =~ /Genomic Features: [0-9]/) {
		($features) = $line =~ /Genomic Features: (.*)/;
	}
	if ($line =~ /N's per 100 kbp: /) {
		($unresolved) = $line =~ /N's per 100 kbp: (.*)/;
	}
	if ($line =~ /Mismatches per 100 kbp: /) {
		($mismatches) = $line =~ /Mismatches per 100 kbp: (.*)/;
	}
	if ($line =~ /Indels per 100 kbp: /) {
		($indels) = $line =~ /Indels per 100 kbp: (.*)/;
	}
	if ($line =~ /5' Ns: /) {
		($prime5) = $line =~ /5' Ns: (.*)/;
	}
	if ($line =~ /3' Ns: /) {
		($prime3) = $line =~ /3' Ns: (.*)/;
	}
	if ($line =~ /Average Depth of Coverage: /) {
		($average) = $line =~ /Average Depth of Coverage: (.*)/;
	}
	if ($line =~ /Fraction with 0 coverage: /) {
		($fold1) = $line =~ /Fraction with 0 coverage: (.*)/;
	}
	if ($line =~ /Fraction with 1x-10x coverage: /) {
		($fold10) = $line =~ /Fraction with 1x-10x coverage: (.*)/;
	}
	if ($line =~ /Fraction with 11x-100x coverage: /) {
		($fold100) = $line =~ /Fraction with 11x-100x coverage: (.*)/;
	}
	if ($line =~ /Fraction with 101x-1000x coverage: /) {
		($fold1000) = $line =~ /Fraction with 101x-1000x coverage: (.*)/;
	}
	if ($line =~ /Fraction with 1001x-2000x coverage: /) {
		($fold2000) = $line =~ /Fraction with 1001x-2000x coverage: (.*)/;
	}
	if ($line =~ /Fraction with 2001x-10000x coverage: /) {
		($fold10000) = $line =~ /Fraction with 2001x-10000x coverage: (.*)/;
	}
	if ($line =~ /Fraction with > 10000x coverage: /) {
		($bigfold) = $line =~ /Fraction with > 10000x coverage: (.*)/;
	}

	if ($line =~ /Variants in Consensus Genome/) {
		while ($line =~ /[a-zA-Z]/) {
			chomp($line=<INPUT>);
			if ($line =~ /[A-Z]/) {
				$mutations .= $line;
			}
		}
		$mutations =~ s/\s+/ /g;
		$mutations =~ s/^\s//g;
	}

#	if ($line =~ /Variants in Read Alignment/) {
#		while ($line =~ /[a-zA-Z]/) {
#			chomp($line=<INPUT>);
#			if ($line =~ /[A-Z]/) {
#				$temp = "$line";
#				$temp =~ s/^\s+//g;
#				$breseq .= "$temp; ";
#			}
#		}
#	}

	if ($line =~ /Quality Control Flags/) {
		while ($line =~ /[a-zA-Z]/) {
			chomp($line=<INPUT>);
			if ($line =~ /Genome Fraction greater/) {
				($qcfraction) = $line =~ /(....)  Genome Fraction greater/
			}
			if ($line =~ /Depth of coverage/) {
				($qccoverage) = $line =~ /(....)  Depth of coverage/
			}
			if ($line =~ /All variants with at least/) {
				($qcvariants) = $line =~ /(....)  All variants with at least/
			}
			if ($line =~ /Reads per base sequence quality/) {
				($qcperbase) = $line =~ /(....)  Reads per base sequence quality/
			}
			if ($line =~ /Sequencing adapter removed/) {
				($qcadapter) = $line =~ /(....)  Sequencing adapter removed/
			}
			if ($line =~ /100x/) {
				($qc100) = $line =~ /(....)  At least/
			}
			if ($line =~ /1000x/) {
				($qc1000) = $line =~ /(....)  At least/
			}
			if ($line =~ /No indels detected/) {
				($indelcheck) = $line =~ /(....)  No indels detected/
			}
			if ($line =~ /Frameshifts in SARS-CoV-2 open reading frames/) {
				($frameshifts) = $line =~ /(....)  Frameshifts in SARS-CoV-2 open reading frames/
			}
		}
	}

}
close (INPUT);

$R1withprimer = "n.d.";
$R2withprimer = "n.d.";
$purgeload = "n.d.";

# breseq

if (-e "breseq/${filebase}_output/index.html") {
	open (INPUT,"< breseq/${filebase}_output/index.html");
	while (defined($line=<INPUT>)) {
		chomp($line);
		@temp = split(/\t/,$line);
		if ($line =~ /-- Evidence --/) {
			($evidence) = $line =~ /center\">(.*)<\/td/;
			$evidence =~ s/RA/SNP/g;
			$evidence =~ s/JC/INDEL/g;
			chomp($line=<INPUT>);
			($position) = $line =~ /right\">(.*)<\/td/;
			chomp($line=<INPUT>);
			$line =~ s/\&rarr;/->/g;
			($mutation) = $line =~ /center\">(.*)<\/td>/;
			chomp($line=<INPUT>);
			($frequency) = $line =~ /right\">(.*)\%<\/td/;
			chomp($line=<INPUT>);
			($description) = $line =~ /center\">(.*)<\/td/;
			@temp2 = split(/<\/font>/,$description);
			if (defined($temp2[0])) {
				if ($temp2[0] =~ /[A-Z][0-9]+[A-Z]$/) {
					($description) = $temp2[0] =~ /([A-Z][0-9]+[A-Z])$/;
				} else {
					$temp2[0] =~ s/<br>/ - /g;
					$temp2[0] =~ s/&nbsp;/ /g;
					$description = $temp2[0];
				}
			} else {
				$description = " ";
			}
			chomp($line=<INPUT>);
			chomp($line=<INPUT>);
			$line =~ s/&#8209;/-/g;
			$line =~ s/<br>/, /g;
			($annotation) = $line =~ /left\">(.*)<\/td/;
			$breseq .= "$evidence $position $mutation $description ($frequency\% of reads) $annotation;";
		}
	}
}

unless ($breseq) {
	$breseq = "n.d.";
}
unless ($breseq =~ /;/) {
	$breseq = "n.d.";
}

# output

print OUTPUT1 "$version";
print OUTPUT1 "\t$qcfraction\t$qccoverage\t$qcvariants\t$qcperbase\t$qcadapter\t$qc100\t$qc1000";
print OUTPUT1 "\t$indelcheck\t$frameshifts";
print OUTPUT1 "\t$rawreadpairs\t$rawbp\t$R1withprimer\t$R2withprimer\t$postprimerreadpairs\t$postprimerbp";
print OUTPUT1 "\t$finalR1";
print OUTPUT1 "\t$fastqcfail\t$fastqcwarn";
print OUTPUT1 "\t$sarsload";
print OUTPUT1 "\t$purgeload";
print OUTPUT1 "\t$length\t$fraction\t$features\t$unresolved\t$mismatches\t$indels";
print OUTPUT1 "\t$prime5\t$prime3";
print OUTPUT1 "\t$average\t$fold1\t$fold10\t$fold100\t$fold1000\t$fold2000\t$fold10000\t$bigfold";
print OUTPUT1 "\tn.d.";
print OUTPUT1 "\t$mutations";
print OUTPUT1 "\t$breseq";
print OUTPUT1 "\n";

close (OUTPUT1);

# exit

system("cp coverage/${filebase}_coverage_plot.png ${filebase}_coverage.png");
system("cp core/${filebase}.consensus.fa ${filebase}.fa");
system("cp ${filebase}_sample.txt ${filebase}_report.txt");

if (-e "breseq/${filebase}_output/index.html") {
	system("cp breseq/${filebase}_output/index.html ${filebase}_mutations.html");
}

open (OUTPUT,"> ${filebase}_report.html");
open (INPUT,"< ${filebase}_sample.html");
while (defined($line=<INPUT>)) {
	chomp($line);
	if ($line =~ /coverage\/${filebase}_coverage_plot.png/) {
		$line =~ s/coverage\/${filebase}_coverage_plot.png/${filebase}_coverage.png/g;
	} elsif ($line =~ /breseq\/${filebase}_output\/index.html/) {
		$line =~ s/breseq\/${filebase}_output\/index.html/${filebase}_mutations.html/g;
	}
	print OUTPUT "$line\n";
}
close (INPUT);
close (OUTPUT);
exit;


 