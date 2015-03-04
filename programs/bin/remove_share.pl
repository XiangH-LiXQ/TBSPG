#!/usr/bin/perl
use strict;
unless(@ARGV==6)
{
    print "Usage:$0 chr_name mit_name chl_name out.sam share noshare\n";
    exit(0);
}


open IN,"$ARGV[0]" or die "$!";
my @chrname=<IN>;
chomp @chrname;
close IN;

open IN,"$ARGV[1]"  or die "$!";
my @mitname=<IN>;
chomp @mitname;
close IN;

open IN,"$ARGV[2]" or die "$!";
my @chlname=<IN>;
chomp @chlname;
close IN;

open IN,"$ARGV[3]" or die "$!";
my @sam=<IN>;
chomp @sam;
close IN;


my %chr;
my %mit;
my %chl;
my $chr;
my $chl;
my $mit;

for $chr(@chrname)
{
    for (@sam)
    {
        if (/$chr/)
        {
            $chr{$_}=1;
        }
    }
}

for $chl(@chlname)
{
    for (@sam)
    {
        if (/$chl/)
        {
            $chl{$_}=1;
        }
    }
}

for $mit(@mitname)
{
    for (@sam)
    {
        if (/$mit/)
        {
            $mit{$_}=1;
        }
    }
}


my @NOS; # for noshare
my @SHA; # for share
my $sam; 

for $sam(@sam)
{


my $head = (split/\t/,$sam)[0]; 
my $ID = (split/\t/,$sam)[12];
    if (($head eq '@SQ')|| ($head eq '@RG') || ($head eq '@PG') || ($ID eq 'XT:A:U'))  
    {
	push @NOS,$sam,"\n";
    }
    elsif ($mit{$sam} && $chl{$sam} && $chr{$sam})
    {
        push @SHA, "mit&chl&chr"."\t".$sam."\n";
    }
    elsif ($mit{$sam} && $chr{$sam})
    {
        push @SHA,  "mit&chr"."\t".$sam."\n";
    }
    elsif($mit{$sam} && $chl{$sam})
    {
        push @SHA,  "chl&mit"."\t".$sam."\n";
    }
    elsif($chl{$sam} && $chr{$sam})
    {
       push @SHA, "chl&chr"."\t".$sam."\n";
    }
    else
    {
        push @NOS,$sam,"\n";
    }
}




open SHA,">$ARGV[4]"||die"Can't open sequence file $ARGV[4]";
print SHA @SHA;
close SHA;

open NOS,">$ARGV[5]"||die"Can't open sequence file $ARGV[5]";
print NOS @NOS;
close NOS;
