#!/usr/bin/perl -w

# NCoV reference is MN908947.3
# Homoplastic sites: https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
# ARTIC3 amplicons: https://github.com/artic-network/artic-ncov2019/blob/master/primer_schemes/nCoV-2019/V3/nCoV-2019.tsv

@homoplastic=("187","1059","2094","3037","3130","6990","8022","10323","10741","11074","13408","14786","19684","20148","21137","24034","24378","25563","26144","26461","26681","28077","28826","28854","29700","4050","13402","11083","15324","21575");

foreach $entry (@homoplastic) {
	$homoplasy{$entry}++;
}

while (defined($line=<DATA>)) {
	chomp($line);
	if ($line =~ /MN908947.3/) {
		@temp = split(/\t/,$line);
		for ($i=($temp[1]+1); $i <= $temp[2]; $i++) {			# BED zero-based start and a one-based end
			$amplicon{$i} .= "amplicon_$temp[3] ";
		}
	}
}

chomp(@files=<*.tsv>);

foreach $file (@files) {
	open (INPUT,"< $file");
	while (defined($line=<INPUT>)) {
		chomp($line);
		if ($line =~ /ref_allele/) {
			next;
		}
		@temp = split(/\t/,$line);
		if ($temp[3] =~ /[ACGT]/) {
			$posfreq{$temp[1]}{$temp[0]}++;			# first key is position, second key is isolate
		} else {
			if ($temp[3] =~ /N/) {
				$Ncount++;
			} else {
				$heterofreq{$temp[1]}{$temp[0]}++;			# first key is position, second key is isolate	
				$heterocount++;	
				$isolateheteropos{$temp[0]}{$temp[1]}++;	# first key is isolate, second key is position	
			}
		}
		$isolates{$temp[0]}++;						# key is isolate
		$basecount++;
	}
	close (INPUT);
}

$count1 = keys %isolates;
print "$count1 isolates in the NCOV-Tools tree\n";
print "Amplicons are from ARTIC version 3, https://artic.network/resources/ncov/ncov-amplicon-v3.pdf\n";
print "Homoplastic sites are from https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473\n\n";

print "Distribution of multiple nucleotides (X) among isolates\n";

foreach $isolate (keys %isolateheteropos) {
	$count = keys %{$isolateheteropos{$isolate}};
	$freqdist{$count}++;
}

foreach $entry (sort {$a <=> $b} keys %freqdist) {
	print "$freqdist{$entry}\tisolates with\t$entry\tX nucelotide(s)\n";
}
print "\n";

foreach $position (keys %posfreq) {
	$count = keys %{$posfreq{$position}};
	if ($count > 1) {
		$notnoise1++;
	}
}

foreach $position (keys %posfreq) {
	$count = keys %{$posfreq{$position}};
	if ($count > 9) {
		$notnoise2++;
	}
}

foreach $position (keys %posfreq) {
	$count = keys %{$posfreq{$position}};
	if ($count > 99) {
		$notnoise3++;
	}
}


print "Positions in the alignment having evidence of multiple nucleotides (X) in 2 or more isolates\n";
foreach $position (sort {$a <=> $b} keys %heterofreq) {
	$icount = keys %{$heterofreq{$position}};
	if ($icount > 1) {
		$nothetnoise1++;
		$note = " ";
		if ($position < 56) {
			$note .= "5' sequencing/mapping error region; ";
		}
		if ($position > 29803) {
			$note .= "3' sequencing/mapping error region; ";
		}
		if ($homoplasy{$position}) {
			$note .= "homoplastic site; ";
		}
		if ($posfreq{$position}) {
			$count = keys %{$posfreq{$position}};
			$note .= "$count isolates with resolved variants";
		} else {
			$note .= "0 isolates with resolved variants";
		}
		print "\t$icount\tisolates \@ position $position  \t$amplicon{$position}\t$note\n";
		foreach $isolate (keys %{$heterofreq{$position}}) {
			$badamplicon{$amplicon{$position}}{$isolate}++;
		}
	}
}

foreach $amplicon (keys %badamplicon) {
	open (OUTPUT,"> $amplicon.txt");
	foreach $isolate (keys %{$badamplicon{$amplicon}}) {
		($name) = $isolate =~ /Consensus_(.*)\.consensus_threshold/;
		print OUTPUT "$name\n";
	}
	close (OUTPUT);
}

print "\nSummary\n";
$count2 = keys %posfreq;
$count3 = keys %heterofreq;
print "\t$Ncount of $basecount variant bases across all isolates have insufficient coverage to make a call (N; <10x coverage)\n";
print "\t$count2 positions have a resolved (ACGT) sequence variant in at least 1 isolate\n";
print "\t$notnoise1 positions have a resolved (ACGT) sequence variant in 2 or more isolates\n";
print "\t$notnoise2 positions have a resolved (ACGT) sequence variant in 10 or more isolates\n";
print "\t$notnoise3 positions have a resolved (ACGT) sequence variant in 100 or more isolates\n";
print "\t$heterocount of $basecount variant bases across all isolates have evidence of multiple nucleotides (X), i.e. at least 75% of the reads support more than one base\n";
print "\t$count3 positions in the alignment have evidence of multiple nucleotides (X) in at least one isolate\n";
print "\t$nothetnoise1 positions in the alignment have evidence of multiple nucleotides (X) in 2 or more isolates (see above)\n";


__DATA__
MN908947.3	54	385	1	nCoV-2019_1	+
MN908947.3	342	704	2	nCoV-2019_2	+
MN908947.3	664	1004	3	nCoV-2019_1	+
MN908947.3	965	1312	4	nCoV-2019_2	+
MN908947.3	1264	1623	5	nCoV-2019_1	+
MN908947.3	1595	1942	6	nCoV-2019_2	+
MN908947.3	1897	2242	7	nCoV-2019_1	+
MN908947.3	2205	2568	8	nCoV-2019_2	+
MN908947.3	2529	2880	9	nCoV-2019_1	+
MN908947.3	2850	3183	10	nCoV-2019_2	+
MN908947.3	3166	3507	11	nCoV-2019_1	+
MN908947.3	3482	3826	12	nCoV-2019_2	+
MN908947.3	3795	4142	13	nCoV-2019_1	+
MN908947.3	4077	4402	14	nCoV-2019_2	+
MN908947.3	4322	4666	15	nCoV-2019_1	+
MN908947.3	4658	4995	16	nCoV-2019_2	+
MN908947.3	4966	5296	17	nCoV-2019_1	+
MN908947.3	5287	5620	18	nCoV-2019_2	+
MN908947.3	5586	5932	19	nCoV-2019_1	+
MN908947.3	5894	6247	20	nCoV-2019_2	+
MN908947.3	6197	6526	21	nCoV-2019_1	+
MN908947.3	6495	6846	22	nCoV-2019_2	+
MN908947.3	6745	7092	23	nCoV-2019_1	+
MN908947.3	7058	7389	24	nCoV-2019_2	+
MN908947.3	7332	7671	25	nCoV-2019_1	+
MN908947.3	7651	7997	26	nCoV-2019_2	+
MN908947.3	7968	8319	27	nCoV-2019_1	+
MN908947.3	8275	8635	28	nCoV-2019_2	+
MN908947.3	8619	8954	29	nCoV-2019_1	+
MN908947.3	8913	9245	30	nCoV-2019_2	+
MN908947.3	9226	9557	31	nCoV-2019_1	+
MN908947.3	9502	9834	32	nCoV-2019_2	+
MN908947.3	9806	10146	33	nCoV-2019_1	+
MN908947.3	10099	10437	34	nCoV-2019_2	+
MN908947.3	10384	10737	35	nCoV-2019_1	+
MN908947.3	10688	11048	36	nCoV-2019_2	+
MN908947.3	11022	11372	37	nCoV-2019_1	+
MN908947.3	11331	11668	38	nCoV-2019_2	+
MN908947.3	11584	11927	39	nCoV-2019_1	+
MN908947.3	11889	12234	40	nCoV-2019_2	+
MN908947.3	12133	12465	41	nCoV-2019_1	+
MN908947.3	12439	12779	42	nCoV-2019_2	+
MN908947.3	12732	13074	43	nCoV-2019_1	+
MN908947.3	13029	13363	44	nCoV-2019_2	+
MN908947.3	13344	13660	45	nCoV-2019_1	+
MN908947.3	13625	13961	46	nCoV-2019_2	+
MN908947.3	13946	14271	47	nCoV-2019_1	+
MN908947.3	14232	14579	48	nCoV-2019_2	+
MN908947.3	14570	14898	49	nCoV-2019_1	+
MN908947.3	14895	15224	50	nCoV-2019_2	+
MN908947.3	15193	15538	51	nCoV-2019_1	+
MN908947.3	15503	15861	52	nCoV-2019_2	+
MN908947.3	15851	16186	53	nCoV-2019_1	+
MN908947.3	16144	16485	54	nCoV-2019_2	+
MN908947.3	16444	16804	55	nCoV-2019_1	+
MN908947.3	16770	17130	56	nCoV-2019_2	+
MN908947.3	17087	17430	57	nCoV-2019_1	+
MN908947.3	17406	17738	58	nCoV-2019_2	+
MN908947.3	17697	18036	59	nCoV-2019_1	+
MN908947.3	17993	18324	60	nCoV-2019_2	+
MN908947.3	18275	18650	61	nCoV-2019_1	+
MN908947.3	18618	18957	62	nCoV-2019_2	+
MN908947.3	18918	19275	63	nCoV-2019_1	+
MN908947.3	19232	19591	64	nCoV-2019_2	+
MN908947.3	19570	19911	65	nCoV-2019_1	+
MN908947.3	19866	20231	66	nCoV-2019_2	+
MN908947.3	20200	20542	67	nCoV-2019_1	+
MN908947.3	20496	20867	68	nCoV-2019_2	+
MN908947.3	20813	21146	69	nCoV-2019_1	+
MN908947.3	21104	21427	70	nCoV-2019_2	+
MN908947.3	21386	21716	71	nCoV-2019_1	+
MN908947.3	21682	22013	72	nCoV-2019_2	+
MN908947.3	21990	22324	73	nCoV-2019_1	+
MN908947.3	22290	22626	74	nCoV-2019_2	+
MN908947.3	22542	22877	75	nCoV-2019_1	+
MN908947.3	22821	23189	76	nCoV-2019_2	+
MN908947.3	23144	23500	77	nCoV-2019_1	+
MN908947.3	23466	23822	78	nCoV-2019_2	+
MN908947.3	23812	24145	79	nCoV-2019_1	+
MN908947.3	24100	24443	80	nCoV-2019_2	+
MN908947.3	24416	24765	81	nCoV-2019_1	+
MN908947.3	24721	25052	82	nCoV-2019_2	+
MN908947.3	25003	25347	83	nCoV-2019_1	+
MN908947.3	25301	25646	84	nCoV-2019_2	+
MN908947.3	25623	25969	85	nCoV-2019_1	+
MN908947.3	25924	26290	86	nCoV-2019_2	+
MN908947.3	26219	26566	87	nCoV-2019_1	+
MN908947.3	26542	26890	88	nCoV-2019_2	+
MN908947.3	26860	27190	89	nCoV-2019_1	+
MN908947.3	27164	27511	90	nCoV-2019_2	+
MN908947.3	27471	27825	91	nCoV-2019_1	+
MN908947.3	27808	28145	92	nCoV-2019_2	+
MN908947.3	28104	28442	93	nCoV-2019_1	+
MN908947.3	28416	28756	94	nCoV-2019_2	+
MN908947.3	28699	29041	95	nCoV-2019_1	+
MN908947.3	29007	29356	96	nCoV-2019_2	+
MN908947.3	29316	29665	97	nCoV-2019_1	+
MN908947.3	29510	29836	98	nCoV-2019_2	+
