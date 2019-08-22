#!/usr/bin/perl
# … 
# Autor: Claudia Alvarez-Carreño
# Fecha: 14 de noviembre de 2013
use strict ;

my $infile = $ARGV[0];


system ("sort -h $infile > $infile.sorted");
open INS, "$infile.sorted" or die "ERROR APERTURA \n$!" ;

my (@B,%BOXB,@PIECES,@BOX_keys,%BOX);
while (my $line = <INS>){
	chomp $line;
	push @PIECES, $line;	
}

my $token = 0;
my $xn    = scalar(@PIECES);

for (my $x = 0 ; $x < $xn; $x++){
	my $ref          = shift @PIECES;
	my $compared_ind = 0;
	foreach my $compared (@PIECES){		
		my $comp = $PIECES[$compared_ind];
		my $add  = class($ref, $PIECES[$compared_ind], $compared_ind);
		$compared_ind ++;
		$compared_ind = $compared_ind + $add ;	
	}
}


foreach my $element (sort {$b<=>$a} keys %BOX){
#foreach my $element (keys %BOX){
	%BOXB="";
	if ($BOX{$element} eq ""){
		if($element ne "") {print "$element\n"}
	}else{
		my @ELEMENTS = split ("-",$BOX{$element});
		my (@INISa,@ENDSa,@Bs);
 		foreach my $e (@ELEMENTS){
			my ($A, $B)           = split("\t",$e);
			my ($ini_ap, $fin_ap) = split (":",$A);
			my ($ini_bp, $fin_bp) = split (":",$B);
			push @INISa, $ini_ap;
			push @ENDSa, $fin_ap;
			push @Bs, $B;
		}
		
		my @INISORTEDa = sort {$a<=>$b} @INISa;
		my @ENDSORTEDa = sort {$a<=>$b} @ENDSa;
		my @Bsa        = sort {$a<=>$b} @Bs;
		
		print "$INISORTEDa[0]:$ENDSORTEDa[-1]";
        
	        arreglab(@Bsa);
		
		print "\n";
	}#if...else
}#foreach

sub arreglab{
	@B     = @_;
	my $xb = scalar(@B);
	my ($token, @INISb, @ENDSb);
	for (my $x = 0 ; $x < $xb  ; $x++){
		if ($B[$x+1] ne ""){
			my ($ini_b1,$fin_b1) = split (":",$B[$x]);
			my $length_b1        = $fin_b1 - $ini_b1;
			my ($ini_b2,$fin_b2) = split (":",$B[$x+1]);
			if ($ini_b2 < $fin_b1+$length_b1){
				push @INISb , ($ini_b1,$ini_b2);
				push @ENDSb , ($fin_b1,$fin_b2);
				$token = "used";
			}else{
				if($x == 0){
					$token = "unsed";
					print "\t$B[$x]";

					if ($x == $xb -2) {print "\t$B[$x+1]"}
					}

			}
			if (($token eq "unused" and  $INISb[0] ne "") or ($token eq "used" and $x == $xb-2) ){
				my @INISORTEDb = sort {$a<=>$b} @INISb;
				my @ENDSORTEDb = sort {$a<=>$b} @ENDSb;
				my (@INISb,@ENDSb);

				print "\t$INISORTEDb[0]:$ENDSORTEDb[-1]";
			
			}
		}
	}
}

sub class{
	my $A                 = @_[0];
	my $B                 = @_[1];
	my $ind               = @_[2];
	my $indicador         = "off";
	my ($Ap, $C)          = split(" ", $A);
	my ($Bp, $D)          = split(" ", $B);
	my ($ini_ap, $fin_ap) = split (":", $Ap);
	my $length_A          = $fin_ap - $ini_ap;
    	my ($ini_bp, $fin_bp) = split (":", $Bp);
	my $length_B          = $fin_bp - $ini_bp;
	
	if ($B eq "") {return(0)}

    	if ($ini_bp < $fin_ap+$length_A){
		$BOX{$A} = "$BOX{$A}$A-$B\t";		
		shift @PIECES,$ind,1;
    	
    		return(-1);

    	}
    
    	$BOX{$A} = "$BOX{$A}";	
    
    	return(0);
}
