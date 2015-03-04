#!/usr/bin/perl
open IN;
while (<STDIN>) {
	chomp;
	$_=~s/\t\t/\t/g;
	my @love = split /\t/ ;
if ($love[0]=~/Chr_name/){
	print "$_\tCoverage\n";
}elsif ($love[0]=~/Locus/){
}else{
	if ($love[2]>0) {
	my $cov=$love[3]/$love[2];
	print "$_\t$cov\n";}
		else{
		print "$_\t0\n";} 
}
}
close IN;

