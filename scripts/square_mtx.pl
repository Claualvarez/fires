#!/usr/bin/perl

use strict ;

my $infile  = $ARGV[0];

open IN, "$infile" or die "\n$!" ;
my (@X,@Y,@MAT);
while (my $lin = <IN>){
	if ($lin =~/^\s$/){
	}else{
		chomp $lin;
		my ($x,$y,undef) = split (" ",$lin);
		$MAT[$x][$y] = 1;
		$MAT[$y][$x] = 1;
		push @X, $x ;
		push @Y, $y ;
	}
}

my $size_x = scalar @X;
my $size_y = scalar @Y;

for (my $x = 0; $x < $size_x+$size_y+1 ; $x ++){
	for (my $y = 0 ; $y < $size_y+$size_x+1 ; $y ++){
		if ($MAT[$y][$x] == 1){
			print "$x $y\n";
		}
	}
}
