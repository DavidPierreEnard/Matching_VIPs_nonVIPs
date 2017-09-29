#!/usr/bin/perl

#use warnings;

use strict;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my ($coord_file,$vip_file,$nonvip_file,$signal_file,$cutoff) = @ARGV;

my %coords;
my %cluster;
my %used;
open(DATA,$coord_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    my $center = int(($splitter_line[2]+$splitter_line[3])/2);
    $coords{$splitter_line[0]} = $splitter_line[1]." ".$center;
    my $int_center = int($center/100000)*100000;
    $cluster{$splitter_line[0]} = $splitter_line[1]." ".$int_center;
    $used{$splitter_line[1]." ".$int_center} = "no";
}
close DATA;

my %signals;
my %vip_scores;
open(DATA,$signal_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    my $chain = "";
    my $sc_sl = scalar(@splitter_line);
    for(my $i=1;$i<=$sc_sl-1;$i++){
	my $ind = $i-1;
	$vip_scores{$ind} = 0;
	$chain .= $splitter_line[$i]." ";
    }
    chop $chain;
    $signals{$splitter_line[0]} .= $chain;
}
close DATA;

my %nonvip_av;
my %nonvip_den;
my $nonvip_lowci = 0;
my $nonvip_highci = 0;
my %pval;

open(DATA,$vip_file);
while(<DATA>){
    chomp $_;
    my $gene = $_;
    my $coord = $cluster{$gene};
    if($used{$coord} eq "no"){
	my $chain = $signals{$gene};
	my @splitter_chain = split(" ",$chain);
	my $sc_sc = scalar(@splitter_chain);
	for(my $i=0;$i<=$sc_sc-1;$i++){
	    if($splitter_chain[$i]>=$cutoff){
		$vip_scores{$i} += 1;
		my @splitter_co = split(" ",$coord);
		my $center = $splitter_co[1];
		for(my $g=0;$g<=5;$g++){
		    my $co_left = $center-100000*$g;
		    my $co_right = $center+100000*($g-1);
		    $used{$splitter_co[0]." ".$co_left} = "yes";
		    $used{$splitter_co[0]." ".$co_right} = "yes";
		}
	    }
	}
    }
}
close DATA;

my @rat = ();
my $total_pval = 0;
my $av_rat = 0;
my $total_den = 0;
my @nonvip = ();

open(DATA,$nonvip_file);
while(<DATA>){
    chomp $_;
    my $key_used;
    foreach $key_used (sort keys %used){
	$used{$key_used} = "no";
    }
    my @splitter_line = split(" ");
    my $sc_sl = scalar(@splitter_line);
    my %nonvip_scores;
    my %reps;
    for(my $i=1;$i<=$sc_sl-1;$i++){
	my $chain = $signals{$splitter_line[$i]};
	my @splitter_chain = split(" ",$chain);
	my $sc_sc = scalar(@splitter_chain);
	my $sign = "no";
	for(my $j=0;$j<=$sc_sc-1;$j++){
	    if($splitter_chain[$j]>=$cutoff){	    
		$sign = "yes";
	    }
	}
	if($sign eq "yes"){
	    $reps{$splitter_line[$i]} += 1;
	}
    }
    my $key_rep;
    foreach $key_rep (shuffle keys %reps){
	my $coord = $cluster{$key_rep};
	if($used{$coord} eq "no"){
	    my $chain = $signals{$key_rep};
	    my @splitter_chain = split(" ",$chain);
	    my $sc_sc = scalar(@splitter_chain);
	    for(my $j=0;$j<=$sc_sc-1;$j++){
		if($splitter_chain[$j]>=$cutoff){
		    $nonvip_scores{$j} += $reps{$key_rep};                                                                                                                                             
		}
	    }
	    my @splitter_co = split(" ",$coord);
	    my $center = $splitter_co[1];
	    for(my $g=0;$g<=5;$g++){
		my $co_left = $center-100000*$g;
		my $co_right = $center+100000*($g-1);
		$used{$splitter_co[0]." ".$co_left} = "yes";
		$used{$splitter_co[0]." ".$co_right} = "yes";
	    }
	}
    }
    my $rat = 0;
    my $rat_den  = 0;
    $total_den += 1;
    my $key_001;
    foreach $key_001 (sort keys %vip_scores){
	if($nonvip_scores{$key_001}>0){
	    $rat += ($nonvip_scores{$key_001});
	}
	else{
	    $rat += 0.1;
	}
	$rat_den += 1;   
	$nonvip_av{$key_001} += $nonvip_scores{$key_001};
	$nonvip_den{$key_001} += 1;
	if($nonvip_scores{$key_001}>=$vip_scores{$key_001}){
	    $pval{$key_001} += 1;
	}
    }
    $rat = $rat/$rat_den;
    $av_rat += $rat;
    push(@rat,$rat);
    if($rat<=1){
	$total_pval += 1;
    }
}
close DATA;

my @sorted_rat = sort{$a<=>$b}@rat;
my $sc_sr = scalar(@sorted_rat);
my $inf_ind = int(0.05*$sc_sr);
my $sup_ind = int(0.95*$sc_sr);
my $inf = $sorted_rat[$inf_ind];
my $sup = $sorted_rat[$sup_ind];

$total_pval = $total_pval/$total_den;
$av_rat = $av_rat/$total_den;
my $line = "";
my $num  = 0;
my $den = 0;
my $ratio = 0;
my $ratio_den = 0;

my $key_01;
foreach $key_01 (sort{$a<=>$b} keys %vip_scores){
    $nonvip_av{$key_01} = $nonvip_av{$key_01}/$nonvip_den{$key_01};
    $pval{$key_01}= $pval{$key_01}/$nonvip_den{$key_01};
    $num += $vip_scores{$key_01};
    $den += $nonvip_av{$key_01};
    if($nonvip_av{$key_01}==0){
	$nonvip_av{$key_01} = 0.1;
    }
    $ratio_den += 1;
    my $excess = $vip_scores{$key_01}/$nonvip_av{$key_01};
    $excess = int($excess*1000)/1000;
    $ratio += $excess;
    $line .= $excess."\t".$vip_scores{$key_01}."\t".$nonvip_av{$key_01}."\t".$inf."\t".$sup."\t".$pval{$key_01};
}

$ratio = $ratio/$ratio_den;
my @sorted_rat = sort{$a<=>$b}@rat;
my $sc_sr = scalar(@sorted_rat);
my $inf_ind = int(0.05*$sc_sr);
my $sup_ind = int(0.95*$sc_sr);
my $inf = $sorted_rat[$inf_ind];
my $sup = $sorted_rat[$sup_ind];
print $line."\n";

