#!/usr/bin/perl -w

#print "Hello, World...\n";

my %hash;
my %sum;
my $file_name = $ARGV[0];
my $one;
my $total;


open (IN, "$file_name");
while (<IN>) {
	chomp;
	$_=~s/:/\t/g;
	my $name=(split(/\t/))[0];
	my $num=(split(/\t/))[2];
if ($num>0) {
	$sum{$name}=$sum{$name}+$num;
	$hash{$name}++;
}
}

close IN;

foreach $one (sort keys %hash) {
	print "Total_locus	$one	$hash{$one}\n";
	}
foreach $total (sort keys %sum) {
	print "Total_Depth	$total	$sum{$total}\n";
	}



