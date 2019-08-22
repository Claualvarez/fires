#!/usr/bin/perl

use strict; use List::Util qw[min max] ;

my (%DIAGONAL,@MAT);
my $in;               #infile name
my $reward      = 2 ;  #rewarding points for a diagonal extension
my $gap_penalty = 1 ;
my $threshold   = 10; #minimum score for a diagonal
my $options     = "";

for (my $i = 0; $i < @ARGV ; $i++) {$options .= " $ARGV[$i] " }
if ($options =~ / -i\s+(\S+) /){$in = $1}
if ($options =~ / -ext\s+(\S+) /){$reward = $1}
if ($options =~ / -p\s+(\S+) /){$gap_penalty = $1}
if ($options =~ / -b\s+(\S+) /){$threshold = $1}

open IN, "$in" or die ;
while (my $lin = <IN>){
        if ($lin =~/^\s$/) {}
	else{
                chomp $lin;
                my ($i, $j)   = split (" ", $lin);
		$MAT[$i][$j]  = 1;
                my $ind       = maximal($i, $j, $gap_penalty, $reward, $threshold);

		store($i,$j);

		if($ind ==1){
                	my @coord ;
                	$coord[0] = $i;
                	$coord[1] = $j;
                	my $name  = "$i:$j";
                	push ( @{$DIAGONAL{$name}}, @coord);
                	push (@{$DIAGONAL{$name}}, $MAT[$i][$j]);
		}
        }
}

sub store{
        my $x = shift (@_);
        my $y = shift (@_);
	foreach my $k (sort {$a<=>$b} keys %DIAGONAL){
                my $anteriorx = ${$DIAGONAL{$k}}[0];
                my $anteriory = ${$DIAGONAL{$k}}[1];
                if ((($x-1 == $anteriorx) or ($x == $anteriorx)) and (($y-1 == $anteriory) or ($y == $anteriory))){
                        my @coord ;
                        $coord[0] = $x ;
                        $coord[1] = $y ;
                        unshift( @{$DIAGONAL{$k}}, @coord );
                        pop @{$DIAGONAL{$k}};
                        push (@{$DIAGONAL{$k}}, $MAT[$x][$y] );
                }
        }
}

sub maximal{
        my $x            = shift (@_);
        my $y            = shift (@_);
	my $gap_penalty  = shift (@_);
        my $reward       = shift (@_);
	my $threshold    = shift (@_);
        my ($ai, $bi, $c)  = 0;

        if ($MAT[$x-1][$y-1] > 0)  {$ai = $MAT[$x-1][$y-1] + $reward}
	elsif ($MAT[$x][$y-1] > 0) {$bi = $MAT[$x][$y-1] - $gap_penalty}
	elsif ($MAT[$x-1][$y] > 0) {$c  = $MAT[$x-1][$y] - $gap_penalty}
        
	$MAT[$x][$y] = max($MAT[$x][$y], $ai, $bi, $c);
        
	return($MAT[$x][$y]);
}
	

foreach my $k (sort {$a<=>$b} keys %DIAGONAL){
        my $baseline = pop @{$DIAGONAL{$k}};
        
        if ($baseline >= $threshold){
                my $ai = 0;
                my $bi = 1;

                if (${$DIAGONAL{$k}}[$ai] ne "" and ${$DIAGONAL{$k}}[$bi] ne ""){
                	my ($vector1_ini,$vector2_ini)= split(":",$k);
                        print "$vector1_ini:${$DIAGONAL{$k}}[$ai]\t$vector2_ini:${$DIAGONAL{$k}}[$bi]\n";
                        $ai ++; $ai ++; $bi ++; $bi ++;
                }
        }
}
