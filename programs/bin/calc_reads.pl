#!/usr/bin/perl -w

#print "Hello, World...\n";

my %hash;
my $file_name = $ARGV[0];
my $one;

open (IN, "$file_name");
while (<IN>) {
	chomp;
	#my $name=$_;
my $first=(split(/\t/))[0];
if (!($first=~/@/g)){
	my $name=(split(/\t/))[2];
	#print "$name\n";
	$hash{$name}++;
}
}
close IN;

foreach $one (sort keys %hash) {
	print "reads_no	$one	$hash{$one}\n";
	}
