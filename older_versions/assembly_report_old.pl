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

open (OUTPUT1,"> ${filebase}_statistics.tsv");
open (OUTPUT2,"> ${filebase}_report.txt");

# header

print OUTPUT1 "Pipeline Version";
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
print OUTPUT1 "\tQC: Genome Fraction less than 90%\tQC: Depth of Coverage less than 2000 fold\tQC: Variant with less than 90% frequency among reads\tQC: Per base sequence quality\tQC: Adapter Content";
print OUTPUT1 "\tQC: At least 90% of positions have better than 100x Coverage\tQC: At least 90% of positions have better than 1000x Coverage";
print OUTPUT1 "\n";

print OUTPUT2 "SARS-CoV-2 Genome Sequencing Consensus & Variants\n";
print OUTPUT2 "https://github.com/jaleezyy/covid-19-sequencing\n";
print OUTPUT2 "\n";

# pipeline version
print OUTPUT1 "468d152";

# input data and primer removal

open (INPUT,"< fastq_trimmed/cutadapt.log");
while (defined($line=<INPUT>)) {
	chomp($line);
	if ($line =~ /Total read pairs processed/) {
		($rawreadpairs) = $line =~ /([0-9].*[0-9])/;
		$rawreadpairs =~ s/,//g;
	}
	if ($line =~ /Read 1 with adapter/) {
		($R1withprimer) = $line =~ /([0-9]+,[0-9]+)/;
		$R1withprimer =~ s/,//g;
	}
	if ($line =~ /Read 2 with adapter/) {
		($R2withprimer) = $line =~ /([0-9]+,[0-9]+)/;
		$R2withprimer =~ s/,//g;
	}
	if ($line =~ /Pairs written/) {
		($postprimerreadpairs) = $line =~ /([0-9]+,[0-9]+)/;
		$postprimerreadpairs =~ s/,//g;
	}
	if ($line =~ /Total basepairs processed/) {
		($rawbp) = $line =~ /([0-9].*[0-9])/;
		$rawbp =~ s/,//g;
	}
	if ($line =~ /Total written/) {
		if ($line =~ /[0-9]+,[0-9]+,[0-9]+/) {
			($postprimerbp) = $line =~ /([0-9]+,[0-9]+,[0-9]+)/;
			$postprimerbp =~ s/,//g;		
		} else {
			($postprimerbp) = $line =~ /([0-9]+,[0-9]+)/;
			$postprimerbp =~ s/,//g;		
		}
	}
}
close (INPUT);

print OUTPUT1 "\t$rawreadpairs\t$rawbp\t$R1withprimer\t$R2withprimer\t$postprimerreadpairs\t$postprimerbp";

print OUTPUT2 "Raw Data: $rawreadpairs read pairs, $rawbp bp\n";

# after Trimmomatic

open (INPUT,"< fastq_trimmed/R1_paired_fastqc/fastqc_data.txt");
while (defined($line=<INPUT>)) {
	chomp($line);
	if ($line =~ /Total Sequences/) {
		($finalR1) = $line =~ /([0-9]+)/;
	}
	if ($line =~ /Sequences flagged as poor quality/) {
		($flaggedR1) = $line =~ /([0-9]+)/;
	}
}
close (INPUT);

open (INPUT,"< fastq_trimmed/R2_paired_fastqc/fastqc_data.txt");
while (defined($line=<INPUT>)) {
	chomp($line);
	if ($line =~ /Total Sequences/) {
		($finalR2) = $line =~ /([0-9]+)/;
	}
	if ($line =~ /Sequences flagged as poor quality/) {
		($flaggedR2) = $line =~ /([0-9]+)/;
	}
}
close (INPUT);

open (INPUT,"< fastq_trimmed/R1_paired_fastqc/summary.txt");
while (defined($line=<INPUT>)) {
	chomp($line);
	@temp = split(/\t/,$line);
	unless ($temp[0] eq "PASS") {
		$fastqc{$temp[0]}{$temp[1]}++;
	}
}
close (INPUT);

$FLAG4 = 0;
$FLAG5 = 0;
open (INPUT,"< fastq_trimmed/R2_paired_fastqc/summary.txt");
while (defined($line=<INPUT>)) {
	chomp($line);
	@temp = split(/\t/,$line);
	unless ($temp[0] eq "PASS") {
		$fastqc{$temp[0]}{$temp[1]}++;
		if ($temp[1] eq "Per base sequence quality" and $temp[0] eq "FAIL") {
			$FLAG4 = 1;
		} elsif ($temp[1] eq "Adapter Content" and $temp[0] eq "FAIL") {
			$FLAG5 = 1;
		}
	}
}
close (INPUT);

if ($flaggedR1 > 0) {
	print STDERR "WARNING flaggedR1\n";
}
if ($flaggedR2 > 0) {
	print STDERR "WARNING flaggedR2\n";
}
unless ($finalR1 == $finalR2) {
	print STDERR "WARNING finalR1 finalR2\n";
}

print OUTPUT1 "\t$finalR1";
print OUTPUT2 "Post-Trim: $postprimerreadpairs read pairs used in consensus generation\n";

foreach $entry (sort keys %fastqc) {
	print OUTPUT1 "\t";
	print OUTPUT2 "\nFASTQC $entry:\n";
	foreach $entry2 (sort keys %{$fastqc{$entry}}) {
		print OUTPUT1 "$entry2; ";
		print OUTPUT2 "\t$entry2\n";
	}
}
print OUTPUT2 "\n";

# Kraken

$sarsload = 0;
open (INPUT,"< kraken2/report");
while (defined($line=<INPUT>)) {
	chomp($line);
	@temp = split(/\t/,$line);
	if ($line =~ /Severe acute respiratory syndrome coronavirus 2/) {
		($sarsload) = $temp[0] =~ /([0-9].*)/;
	}
}
close (INPUT);

print OUTPUT1 "\t$sarsload";
print OUTPUT2 "Reads SARS-CoV-2 (Kraken2): ${sarsload}\%\n";
print OUTPUT2 "\n";

# Human purge

$purgeload = 0;
open (INPUT,"< host_removed/hisat2.log");
while (defined($line=<INPUT>)) {
	chomp($line);
	if ($line =~ /overall alignment rate/) {
		($purgeload) = $line =~ /(.*)\%/;
	}
}
close (INPUT);

print OUTPUT1 "\t$purgeload";

# QUAST

open (INPUT,"< quast/report.tsv");
while (defined($line=<INPUT>)) {
	chomp($line);
	@temp = split(/\t/,$line);
	if ($line =~ /Total length \(>= 0 bp\)/) {
		$length = $temp[1];
	}
	if ($line =~ /genomic features/) {
		$features = $temp[1];
	}
	if ($line =~ /Genome fraction/) {
		$fraction = $temp[1];
	}
	if ($line =~ /N's per 100 kbp/) {
		$unresolved = $temp[1];
	}
	if ($line =~ /mismatches per 100 kbp/) {
		$mismatches = $temp[1];
	}
	if ($line =~ /indels per 100 kbp/) {
		$indels = $temp[1];
	}
}
close (INPUT);

print OUTPUT1 "\t$length\t$fraction\t$features\t$unresolved\t$mismatches\t$indels";
print OUTPUT2 "Genome Statistics (QUAST & HiSAT2):\n";
print OUTPUT2 "\tGenome Length: $length bp\n";
print OUTPUT2 "\tGenome Fraction: ${fraction}\%\n";
print OUTPUT2 "\tGenomic Features: $features\n";
print OUTPUT2 "\tN's per 100 kbp: $unresolved\n";
print OUTPUT2 "\tMismatches per 100 kbp: $mismatches\n";
print OUTPUT2 "\tindels per 100 kbp: $indels\n";

$FLAG1 = 0;
if ($fraction < 90) {
	$FLAG1 = 1;
}

# consensus

open (INPUT,"< consensus/virus.consensus.fa");
chomp($line=<INPUT>);
chomp($line=<INPUT>);
close (INPUT);
($temp1) = $line =~ /^(N+)/;
($temp2) = $line =~ /(N+)$/;
$prime5 = length($temp1);
$prime3 = length($temp2);

unless (defined($prime5)) {
	$prime5 = 0;
}
unless (defined($prime3)) {
	$prime3 = 0;
}

print OUTPUT1 "\t$prime5\t$prime3";

# coverage

$fold1++;
$fold10++;
$fold100++;
$fold1000++;
$fold2000++;
$fold10000++;
$bigfold++;
if (-e "coverage/depth.txt") {
	open (INPUT,"< coverage/depth.txt");
	while (defined($line=<INPUT>)) {
		chomp($line);
		@temp = split(/\t/,$line);
		$count++;
		$sum += $temp[2];
		if ($temp[2] < 1) {
			$fold1++;
		} elsif ($temp[2] <= 10) {
			$fold10++;		
		} elsif ($temp[2] <= 100) {
			$fold100++;		
		} elsif ($temp[2] <= 1000) {
			$fold1000++;		
		} elsif ($temp[2] <= 2000) {
			$fold2000++;		
		} elsif ($temp[2] <= 10000) {
			$fold10000++;		
		} else {
			$bigfold++
		}
	}
	close (INPUT);
	$average = int(($sum / $count)*100) / 100;
	$fold1= int($fold1 / $count * 10000)/10000;
	$fold10= int($fold10 / $count * 10000)/10000;
	$fold100= int($fold100 / $count * 10000)/10000;
	$fold1000= int($fold1000 / $count * 10000)/10000;
	$fold2000= int($fold2000 / $count * 10000)/10000;
	$fold10000= int($fold10000 / $count * 10000)/10000;
	$bigfold= int($bigfold / $count * 10000)/10000;
} else {
	$average = 0;
	$fold1= 0;
	$fold10= 0;
	$fold100= 0;
	$fold1000= 0;
	$fold2000= 0;
	$fold10000= 0;
	$bigfold = 0;
}

print OUTPUT1 "\t$average\t$fold1\t$fold10\t$fold100\t$fold1000\t$fold2000\t$fold10000\t$bigfold";

print OUTPUT2 "\tAverage Depth of Coverage: $average\n";
print OUTPUT2 "\t\tFraction with 0 Coverage: $fold1\n";
print OUTPUT2 "\t\tFraction 1x-10x Coverage: $fold10\n";
print OUTPUT2 "\t\tFraction 11x-100x Coverage: $fold100\n";
print OUTPUT2 "\t\tFraction 101x-1000x Coverage: $fold1000\n";
print OUTPUT2 "\t\tFraction 1000x-2000x Coverage: $fold2000\n";
print OUTPUT2 "\t\tFraction 2000x-10000x Coverage: $fold10000\n";
print OUTPUT2 "\t\tFraction with greater than 10000x Coverage: $bigfold\n";
print OUTPUT2 "\t5' Ns: $prime5\n";
print OUTPUT2 "\t3' Ns: $prime3\n";

$FLAG2 = 0;
if ($average < 2000) {
	$FLAG2 = 1;
}

$FLAG6 = 0;
$FLAG7 = 0;
unless (($bigfold + $fold10000 + $fold2000 + $fold1000) >= .9) {
	$FLAG6 = 1;	# At least 90% of positions have better than 100x Coverage
}
unless (($bigfold + $fold10000 + $fold2000) >= .9) {
	$FLAG7 = 1;	# At least 90% of positions have better than 1000x Coverage
}

# LMAT

open (INPUT,"< lmat/parseLMAT_output.txt");
chomp($line=<INPUT>);
chomp($line=<INPUT>);
close (INPUT);
@temp = split(/\t/,$line);

if (defined($temp[5])) {
	print OUTPUT1 "\t$temp[5]";
} else {
	print OUTPUT1 "\t";
}

# iVAR

open (INPUT,"< ivar_variants/ivar_variants.tsv");
while (defined($line=<INPUT>)) {
	chomp($line);
	@temp = split(/\t/,$line);
	if ($line =~ /ALT_FREQ/) {
		next;
	}
	if ($temp[3] =~ /[A-Z0-9]/) {
		$mutation = "${temp[2]}${temp[1]}${temp[3]}";
		$mutations .= "$mutation ";
	}
}

print OUTPUT1 "\t$mutations";
print OUTPUT2 "\n";
print OUTPUT2 "Variants in Consensus Genome (iVar):\n";
print OUTPUT2 "\t$mutations\n";
print OUTPUT2 "\n";

# BreSeq

print OUTPUT2 "Variants in Read Alignment (BreSeq):\n";

$FLAG3=0;
$breseq = " ";
open (INPUT,"< breseq/output/index.html");
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
		print OUTPUT2 "\t$evidence $position $mutation $description ($frequency\% of reads) $annotation\n";
		if ($frequency < 90) {
			$FLAG3 = 1;
		}
	}
}

print OUTPUT1 "\t$breseq";

# QC

print OUTPUT1 "\t$FLAG1\t$FLAG2\t$FLAG3\t$FLAG4\t$FLAG5\t$FLAG6\t$FLAG7";

print OUTPUT2 "\n";
print OUTPUT2 "Quality Control Flags:\n";
if ($FLAG1 == 0) {
	print OUTPUT2 "\tPASS\tGenome Fraction greater than 90%\n";
} else {
	print OUTPUT2 "\tFAIL\tGenome Fraction greater than 90%\n";
}
if ($FLAG2 == 0) {
	print OUTPUT2 "\tPASS\tDepth of Coverage at least 2000 fold\n";
} else {
	print OUTPUT2 "\tFAIL\tDepth of Coverage at least 2000 fold\n";
}
if ($FLAG3 == 0) {
	print OUTPUT2 "\tPASS\tAll variants with at least 90% frequency among reads\n";
} else {
	print OUTPUT2 "\tWARN\tAll variants with at least 90% frequency among reads\n";
}
if ($FLAG4 == 0) {
	print OUTPUT2 "\tPASS\tReads per base sequence quality\n";
} else {
	print OUTPUT2 "\tFAIL\tReads per base sequence quality\n";
}
if ($FLAG5 == 0) {
	print OUTPUT2 "\tPASS\tSequencing adapater removed\n";
} else {
	print OUTPUT2 "\tFAIL\tSequencing adapater removed\n";
}
if ($FLAG6 == 0) {
	print OUTPUT2 "\tPASS\tAt least 90% of positions have better than 100x Coverage\n";
} else {
	print OUTPUT2 "\tFAIL\tAt least 90% of positions have better than 100x Coverage\n";
}
if ($FLAG7 == 0) {
	print OUTPUT2 "\tPASS\tAt least 90% of positions have better than 1000x Coverage\n";
} else {
	print OUTPUT2 "\tWARN\tAt least 90% of positions have better than 1000x Coverage\n";
}

# done reports

print OUTPUT2 "\n";
print OUTPUT2 "Provided Read Coverage Plot (BreSeq):\n";
print OUTPUT2 "\tUnique: reads with only one best match to the SARS-CoV-2 genome\n";
print OUTPUT2 "\tRepeat: reads with multiple equally good matches to repeat sequences\n";
print OUTPUT2" \tTop & Bottom: mapping to the forward and reverse complement matches, respectively\n";

print OUTPUT1 "\n";
close (OUTPUT1);
close (OUTPUT2);

# images

system("cp breseq/output/evidence/MN908947.3.overview.png ${filebase}_coverage.png");
system("cp consensus/virus.consensus.fa ${filebase}.fa");
system("cp breseq/output/index.html ${filebase}_mutations.html");


 