#!/usr/bin/perl


open (SEQ, $ARGV[0]) || die "Can't open $ARGV[0]\n";
$/="\@HWI";<SEQ>;
while (<SEQ>){
	chomp;
	my ($name,$seq,$add,$quality)=split("\n", $_, 4);
	my $n1=(split(" ",$name))[0]."_part1";
	my $n2=(split(" ",$name))[0]."_part2";	
	my $num=length($seq);
if ($num%2==0){
	my $s1=substr($seq,0,$num/2);
	my $s2=substr($seq,$num/2,$num/2);
	my $q1=substr($quality,0,$num/2);
	my $q2=substr($quality,$num/2,$num/2);
	print"\@HWI$n1\n$s1\n+\n$q1\n";
	print"\@HWI$n2\n$s2\n+\n$q2\n";
}else{
	my $s1=substr($seq,0,$num/2+1);
	my $s2=substr($seq,$num/2+1,$num/2+1);
	my $q1=substr($quality,0,$num/2+1);
	my $q2=substr($quality,$num/2+1,$num/2+1);
	print"\@HWI$n1\n$s1\n+\n$q1\n";
	print"\@HWI$n2\n$s2\n+\n$q2\n";
}
}
$/="\n";
close(SEQ);

