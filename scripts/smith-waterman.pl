#!/usr/bin/perl
# Author: Claudia Alvarez-Carreno
use strict; use List::Util qw[min max] ;

my $in = $ARGV[0];
my (@X,@Y,%DIAGONAL,%DIAG,$id,@MAT,$delta,$size,$ai,$bi,$c,$cx,$cy,$sumax,$sumay);
$delta  = 2;
my $gap = 1;
my $gap_size = 10;

my ($last_x,$last_y)="FIRST";
open IN, "$in" or die "\n$!" ;

while (my $lin = <IN>){
        if ($lin =~/^\s$/){
        }else{
                chomp $lin;
                my ($i,$j,undef) = split (" ",$lin,3);
		$MAT[$i][$j]=1;
                my $ind = maximal($i,$j);
		store($i,$j);
		if($ind ==1){
                	my @coord ;
                	$coord[0] = $i ;
                	$coord[1] = $j ;
                	my $name = "$i:$j" ;
                	push ( @{$DIAGONAL{$name}}, @coord );
                	push (@{$DIAGONAL{$name}}, $MAT[$i][$j] );
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
        my $x = shift (@_);
        my $y = shift (@_);
	my $gap_penalty = 1;
        $delta  = 2;
        ($ai,$bi,$c)=0;
        if ($MAT[$x-1][$y-1] > 0){
                $ai = $MAT[$x-1][$y-1] + $delta ;
        }elsif ($MAT[$x][$y-1] > 0){
                $bi = $MAT[$x][$y-1] - $gap_penalty ;
        }elsif ($MAT[$x-1][$y] > 0){
                $c = $MAT[$x-1][$y] - $gap_penalty ;
        }
        $MAT[$x][$y] = max($MAT[$x][$y],$ai,$bi,$c);
	return($MAT[$x][$y]);
}
	

foreach my $k (sort {$a<=>$b} keys %DIAGONAL){
        my $threashold = pop @{$DIAGONAL{$k}};
        if ($threashold >= 18){
                $ai = 0;
                $bi = 1;
                foreach (values @{$DIAGONAL{$k}}){
                        if ((${$DIAGONAL{$k}}[$ai] ne "") and ( ${$DIAGONAL{$k}}[$bi] ne "")){
				my ($vector1_ini,$vector2_ini)= split(":",$k);
				print "${$DIAGONAL{$k}}[$ai]\t${$DIAGONAL{$k}}[$bi]\t$vector2_ini:${$DIAGONAL{$k}}[1]\n";
                        }
			$ai ++; $ai ++; $bi ++; $bi ++;
                }
        }
}



