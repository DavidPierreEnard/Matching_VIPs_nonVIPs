#!/usr/bin/perl

#use warnings;

use strict;

my ($valid_file,$segment_file,$out_file,$cutoff) = @ARGV;

my %ens;

open(DATA,$valid_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    $ens{$splitter_line[0]} = "0";
}
close DATA;

open(DATA,$segment_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split("\t",$_);
    if(($splitter_line[1]>=$cutoff)and($splitter_line[1]<=500000)){
	$ens{$splitter_line[0]} = "1";
    } 

}
close DATA;

open OUT, ">".$out_file;
my $key_01;
foreach $key_01 (sort keys %ens){
    print OUT $key_01."\t".$ens{$key_01}."\n";
}
close OUT;
