#!/usr/bin/perl

use strict ;


my $infile  ;
my $inref   ;
my $pdbcode ;
my $options = "";


for (my $i = 0; $i < @ARGV ; $i++) {$options .= " $ARGV[$i] " }
if ($options =~ / -i\s+(\S+) /)    {$infile   = $1}
if ($options =~ / -dssp\s+(\S+) /) {$inref    = $1}
if ($options =~ / -pdb\s+(\S+) /)  {$pdbcode  = $1}


open IN, "$infile" or die "ERROR APERTURA \n$!" ;

my ($name,@PIECES);

while (my $line = <IN>){
	push @PIECES, $line;	
}

my%CODE=(
	"A"=>"ALA", "R"=>"ARG", "N"=>"ASN", "D"=>"ASP",
	"C"=>"CYS", "E"=>"GLU", "Q"=>"GLN", "G"=>"GLY",
	"H"=>"HIS", "I"=>"ILE", "L"=>"LEU", "K"=>"LYS",
	"M"=>"MET", "F"=>"PHE", "P"=>"PRO", "S"=>"SER",
	"T"=>"THR", "W"=>"TRP", "Y"=>"TYR", "V"=>"VAL"
);

my $xn = scalar(@PIECES);

for (my $x = 0 ; $x < $xn; $x++){
	my @INDEXES = split (" ",$PIECES[$x]);
	my $number  = 0;
	foreach my $index (@INDEXES){
		my ($ini,$fin) = split (":",$index);
		if ($number == 0){
			$name = $pdbcode;
			
			print "$name\t";
		}else{
			
			print "$name\t";
		}
		$name   = cut($ini,$fin);		
		$number ++;
	}
		
}
	

sub cut{
	my $ini = @_[0];
	my $fin = @_[1];

	open REF, "$inref" or die "ERROR APERTURA \n$!" ;
	my $initial;
	while (my $ref = <REF>){
		chomp $ref;
		if ($ref =~/\#/) {$initial ="INICIO"}
		if (($initial eq "INICIO") and ($ref =~/^\s+[0-9]+/)){
			$ref =~/.(....).(....).(.).(.)..(.)/;
			my $nb     = $1; 
			my $res_nb = $2; 
			my $chain  = $3; 
			my $aa     = $4; $aa=~s/[a-z]/C/g;
			my $str    = $5; $str=~s/ /-/;$str=~s/[GI]/H/;$str=~s/E/B/;$str=~s/S/T/; 
			$res_nb    =~s/ //g; 
			if ($ref =~/!/) {next}
			if ($res_nb eq $ini){
				
				print "$chain\t$CODE{$aa}\t$res_nb\t";
				$name = "t_$name\_$chain\_$CODE{$aa}$res_nb";
			}elsif($res_nb eq $fin){
			
				print "$chain\t$CODE{$aa}\t$res_nb\n";
				$name = "$name\_$CODE{$aa}$res_nb";
				close REF;

				return ($name);
				
			}
		}
	}
}


