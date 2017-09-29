#!/usr/bin/perl

#use warnings;

use strict;

my ($indir,$valid_file,$vip_file,$rec_file,$factors_table_file,$distance_file,$intervals_file,$iterations,$out_file_1,$out_file_2,$dist,$cutrec) = @ARGV;

my %valid;

open(DATA,$indir.$valid_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    $valid{$splitter_line[0]} = "yes";
}
close DATA;

my %rec;
open(DATA,$indir.$rec_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    $rec{$splitter_line[0]} = $splitter_line[1];
}
close DATA;

my $factor_number = 0;
my %factors;
open(DATA,$indir.$factors_table_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    my $gene = $splitter_line[0];
    shift @splitter_line;
    $factor_number = scalar(@splitter_line);
    $factors{$gene} = "@splitter_line";
}
close DATA;

my %include;
my $key_001;
foreach $key_001 (sort keys %valid){
    if(($rec{$key_001}>=$cutrec)and($factors{$key_001} ne "")){
	$include{$key_001} = "yes";
    }
}

my %distance;
open(DATA,$indir.$distance_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    $distance{$splitter_line[0]} = $splitter_line[1];
}
close DATA;

my %vips;
my %nonvips;
my %avail;
my $vip_number = 0;
my @nonvips = ();
my $nvip_number = 0;
open(DATA,$indir.$vip_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    if($include{$splitter_line[0]} eq "yes"){
	if($splitter_line[1] eq "yes"){
	    $vips{$splitter_line[0]} = $factors{$splitter_line[0]};
	    $vip_number += 1;
	}
	if(($splitter_line[1] eq "no")and($distance{$splitter_line[0]}>=$dist)){
	    $nonvips{$splitter_line[0]} = $factors{$splitter_line[0]};
	    $avail{$splitter_line[0]} = "yes";
	    push(@nonvips,$splitter_line[0]);
	    $nvip_number += 1;
	}
    }
}
close DATA;

#print $vip_number."\t".$nvip_number."\n";

my $fact_counter = 0;
my %lower_bound;
my %upper_bound;
open(DATA,$intervals_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split("\t",$_);
    $fact_counter += 1;
    $lower_bound{$fact_counter} = $splitter_line[1];
    $upper_bound{$fact_counter} = $splitter_line[2];
}
close DATA;


my $num=0;
my $den=0;

my %matches;
my %controls;
my %kept_vips;

my $key_01;
foreach $key_01 (sort keys %vips){
    my $vip = $key_01;
    my $vip_factors_chain = $vips{$vip};
    #for the specific case of Neanderthal to modern human introgressions, the factor chain includes seven factors in this specific order:
    #CDS density, DNASEI density, FUNSEQ, GC content, recombination, Tajima's D and PhastCons conserved elements density.
    #The definition of tolerated intervals (see methods) is given by $intervals_file. You need the same number of factors in $intervals_file
    #as you have in the factors table. Line 1 in $intervals_file must correspond to the factor in column 1 of factors table, and so on.
    
    my %vip_factors;
    my $counter_001 = 0;
    my @splitter_chain = split(" ",$vip_factors_chain);
    my $sc_sc = scalar(@splitter_chain);
    for(my $i=0;$i<=$sc_sc-1;$i++){
	$counter_001 += 1;
	my $factor_value = $splitter_chain[$i];
	$vip_factors{$counter_001} = $factor_value;
    }


    my $key_02;
    foreach $key_02 (sort keys %nonvips){
	my $nonvip = $key_02;
	my $nonvip_factors_chain = $nonvips{$nonvip};
	my %nonvip_factors;
	my $counter_002 = 0;
	my @splitter_chain = split(" ",$nonvip_factors_chain);
	my $sc_sc = scalar(@splitter_chain);
	for(my $i=0;$i<=$sc_sc-1;$i++){
	    $counter_002 +=1;
	    my $factor_value = $splitter_chain[$i];
	    $nonvip_factors{$counter_002} = $factor_value;
	}

	my $matching_number = 0;
	my $key_003;
	foreach $key_003 (sort keys %vip_factors){
	    if(($nonvip_factors{$key_003}>=(1-$lower_bound{$key_003})*$vip_factors{$key_003})and($nonvip_factors{$key_003}<=(1+$upper_bound{$key_003})*$vip_factors{$key_003})){
		$matching_number += 1;
	    }
	}

	if($matching_number==$factor_number){
	    $matches{$vip} .= $nonvip. " ";
	    $controls{$nonvip} += 1;
	}
    }

    chop $matches{$vip};
    my @splitter_chain2 = split(" ",$matches{$vip});
    my $sc_ch = scalar(@splitter_chain2);
    if($sc_ch>=3){
	$num += 1;
	$kept_vips{$vip} = "yes";
    }
    $den += 1;

    #print $num."\t".$den."\n";
}

open OUT, ">".$out_file_1;
my $key_04;
foreach $key_04 (sort keys %kept_vips){
    print OUT $key_04."\n";
}
close OUT;

open OUT, ">".$out_file_2;
for(my $r=1;$r<=$iterations;$r++){

    my %new_avail;
    my $key_av;
    foreach $key_av (sort keys %avail){
	$new_avail{$key_av} = $avail{$key_av};
    }

    my $line = "sample_".$r;
    my %rep;
    my $repeats = 0;
    my $key_001;
    foreach $key_001 (sort keys %kept_vips){
        my $matches = $matches{$key_001};
        my @splitter = split(" ",$matches);
        my $sc_sp = scalar(@splitter);
        for(my $i=1;$i<=200;$i++){
            my $rand_ind = int(rand($sc_sp));
            my $rand_gene = $splitter[$rand_ind];
            if($new_avail{$rand_gene} eq "yes"){
                $line .= " ".$rand_gene;
                $rep{$rand_gene} += 1;
                #$new_avail{$rand_gene} = "no";
                last;
            }
            if($i==200){
                $repeats += 1;
                $line .= " ".$rand_gene;
                #$new_avail{$rand_gene} = "no";
            }
        }
    }
    my $used = 0;
    my $key_rep;
    foreach $key_rep (sort keys %rep){
        if($rep{$key_rep}>=1){
            $used += 1;
        }
    }
    print $r."\t".$used."\n";
    print OUT $line."\n";
}
