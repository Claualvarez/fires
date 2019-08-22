#!/usr/bin/perl
# … uso: ~/Desktop/Bioinformatics/Scripts/perl/get_FASTA_from_kegg.pl <lista> <salida>
# Autor: Claudia Alvarez-Carreño
# Fecha: 19 enero 2014
use strict ; 
use List::Util qw[min max] ;
#use SSDivide;

my $infile      ; #= $ARGV[0];
my $infile2     ; #= $ARGV[1];
my $infile3     ;
my $sampling   = 30;
my $redundancy = 150 ;
my $options    = "";



for (my $i = 0; $i < @ARGV ; $i++) {$options .= " $ARGV[$i] " }
if ($options =~ / -i\s+(\S+) /){$infile     = $1} 
if ($options =~ / -j\s+(\S+) /){$infile2    = $1}
if ($options =~ / -k\s+(\S+) /){$infile3    = $1}
if ($options =~ / -R\s+(\S+) /){$redundancy = $1}
if ($options =~ / -S\s+(\S+) /){$sampling   = $1}


my (@MAT,@X,@Y,@SIZEX,@SIZEY);
save($infile,0);
save($infile2,1);


my $size_x = max($SIZEX[0],$SIZEX[1]);
my $size_y = max($SIZEY[0],$SIZEY[1]);
my $total;
for (my $x = 0 ; $x < $size_x; $x ++){
	for (my $y = 0 ; $y < $size_y; $y ++){
#		my $value =  $MAT[1][$x][$y] - $MAT[0][$x][$y];
		if($MAT[0][$x][$y] == 5){
			my $value =  $MAT[0][$x][$y] - $MAT[1][$x][$y];
#			print "$x\t$y\t$value\n";
			$total = $value + $total ;
		}
		#if ($MAT[0][$x][$y] == 5){
		#	$total ++;
		#}
	}
	#print "\n";
}

open KMER, "$infile3" or die "ERROR APERTURA \n$!" ;
my (@SAMPLED,@SAMPLE);
while (my $line = <KMER>){
	chomp $line;
	my ($coord1,$coord2) = split (/\t/,$line);
	$SAMPLED[$coord1] ++;
	$SAMPLED[$coord1] ++;
}

foreach my $element (keys @SAMPLED){
	if ($SAMPLED[$element] =~ /^$/){}
	else{push @SAMPLE, $SAMPLED[$element]}
	#print "$element $SAMPLED[$element]\n";
}
my @sSAMPLE = sort {$a <=> $b} (@SAMPLE);
my $ult = scalar @SAMPLE;
#print "primero de la sample: $sSAMPLE[$0]\t ultimo $sSAMPLE[$ult-1]\n";
#<STDIN>;

if ($total < $sampling or $sSAMPLE[$ult-1] > $redundancy) {print "NO\n"}
else{ print "$total\n"}


sub save{
        my $in=@_[0];
        my $ref=@_[1];

        open IN, "$in" or die "ERROR APERTURA \n$!" ;
        while (my $lin = <IN>){
                if ($lin =~/^\s$/){
                }else{
                        chomp $lin;
                        my ($x,$y,undef) = split (" ",$lin);
        #               print "$ref $x $y\n";
                        push @X, $x ;
                        push @Y, $y ;
                        $MAT[$ref][$x][$y] = 5;
                        #print "$ref $x $y $MAT[$ref][$x][$y]\n";
                }
        }

        my $size_x = scalar @X;
        push @SIZEX, $size_x;
        my $size_y = scalar @Y;
        push @SIZEY, $size_y;
        close IN;
        
}

exit;
