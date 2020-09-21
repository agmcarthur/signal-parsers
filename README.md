## SIGNAL Processing Scripts

Processing of SIGNAL (https://github.com/jaleezyy/covid-19-signal) output.

Inside a directory of a single isolate's analysis:

* assembly_report.pl, generate report files for a single isolate

Inside a folder containing directories of many isolate results:

* assembly_summarize.pl, summarize the assembly_report.pl results of all the isolates
* assembly_mutations.pl, summarize the mutations of all the isolates
* prepbioproject.pl, create a directory of all release ready FASTQ for a set of isolates

Run all scripts in bulk inside a folder containing directories of many isolate results:

* sarsreportbulk.sh, run all of the above on a set of isolates

Generate plots using the following script, once you have run assembly_summarize.pl. You can include a TSV file as the single argument that contains sampleID and Ct value to generate Ct plots:

* sarsplots.pl, run all of the above on a set of isolates

## Output

The only output to the screen is some summary statistics from sarsplots.pl, e.g.

> 9 of 18 PASS Genome Fraction of at least 90%

> 1 of 18 PASS 90% of positions have at least 100x coverage

> 1 of 18 PASS 90% of positions have at least 1000x coverage

Within the folder of a single isolate, assembly_report.pl will generate, e.g.

> S1012.fa (consensus genome sequence)

> S1012_coverage.png (plot of read coverage along the genome)

> S1012_mutations.html (HTML version of BreSeq read mapping mutation analysis)

> S1012_report.html (HTML version of full report for isolate, which include coverage plot and BreSeq results)

> S1012_report.txt (text version of full report for isolate)

> S1012_statistics.tsv (tab-delimited full statistical report for isolate)

Summarizing across many isolate folders will produce:

> results.tsv (concatenated statistics for all isolates)

> results_mutations.tsv (tab-delimited report of all mutations across all isolates)

> results_mutations_freq.tsv (tab-delimited report of frequency of mutations across all isolates; does not report % of reads supporting the mutations)

> bioproject_fastq (folder with host-removed FASTQ files)

> results (copies of all the above single isolate files in one directory)

> results_plots (tab-delimited files for graphical reports, see below)

## Plots

Plots can be created from tab-delimited files at https://scatterplot.online. Note that many plots use final consensus genome fraction, for which our QC PASS is at least 90%.

# overview.tsv

Plot of %SARS-CoV-2 in sequencing reads against final consensus genome fraction, labeled with isolate name and average fold coverage of the genome. Axes should run from 0 to 100.

![overview.tsv](/images/overview.jpg)

# ctplot.tsv

Plot of RT-PCR Ct value against final consensus genome fraction, labeled with isolate name and average fold coverage of the genome. Y-axis should run from 0 to 100, X-axis should start at 0.

![overview.tsv](/images/ctplot.jpg)

# ave_coverage.tsv

Plot of average genome coverage against final consensus genome fraction, labeled with isolate name. Y-axis should run from 0 to 100, X-axis should start at 0. While average depth of coverage >2000x was an early cut-off based on MiniSeq over-sequencing, we now consider >= 100x coverage for at least 90% of positions a more reliable QC metric.

![overview.tsv](/images/ave_coverage.jpg)

# 100x_coverage.tsv

Plot of fraction of genome with >100x coverage against final consensus genome fraction, labeled with isolate name. Y-axis should run from 0 to 100, X-axis should run from 0 to 1. Sequences passing coverage QC should have >90% of bases with at least 100x coverage.

![overview.tsv](/images/100x_coverage.jpg)

## SIGNAL Supplementary Scripts

* ncovalleles.pl, summarizes mutations in the an NCoV-Tools tree, reading the output tsv files. Can accept a tab-delimited files of Ct values as a single argument

* packagencov.sh, a custom script for use on the ascension server to package up NCoV-Tools results for download

* packageresult.sh, a custom script for use on the ascension server to package up SIGNAL results for download
