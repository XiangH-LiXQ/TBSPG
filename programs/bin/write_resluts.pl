#!/usr/bin/perl

my %hash;
#my $file_name = $ARGV[0];
my $one;

open IN;
while (<STDIN>) {
	chomp;
	my $first=(split(/\t/))[1];
	my $name=(split(/\t/))[2];
	$hash{$first}="$hash{$first}"."\t$name";
}
print "Chr_name\t\tReads_No\tCovered_ref_size\tDepth_sum\n";
close IN;

foreach $one (sort keys %hash) {
	print "$one	$hash{$one}\n";
	}
