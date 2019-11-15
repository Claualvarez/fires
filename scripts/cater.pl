#!/usr/bin/perl
# … 
# Autor: Claudia Alvarez-Carreño

use strict ;

my $options  ;
my $infile   ;
my $stored   ;
my $chain   ;
my $outname;
my $reference ;
my $chain_name ;
my $pdb ;

for (my $i = 0; $i < @ARGV ; $i++) {$options .= " $ARGV[$i] " }
if ($options =~ / -i\s+(\S+) /)    {$infile   = $1}
if ($options =~ / -c\s+(\S+) /)    {$chain    = $1}
if ($options =~ / -ch\s+(\S+) /)   {$chain_name    = $1}
if ($options =~ / -pdb\s+(\S+) /)  {$pdb    = $1}
if ($options =~ / -out\s+(\S+) /)  {$outname    = $1}
chomp $infile;


#print "$chain_name----\n";
#my ($pdb,undef,undef) = split (/\./,$infile);

open IN, "$infile" or die "ERROR APERTURA \n$!" ;
my ($number,$nb) = 0;
while (my $line = <IN> ){

	if ($line =~ /^([0-9]+)\..*Query: \[(.*)\].*:/){
		$nb = $1;
		$number = 0;

		my $range = $2;
		my ($ini_pos,$last_pos) = split ("-",$range);
		open OUT, ">$outname\_$chain_name\_$nb.pdb";
		open INPDB, "$pdb\_$chain_name.pdb" or die "\n$!";;
		
		while (my $l = <INPDB>){

		        if ($l =~/^ATOM\s+[0-9]+\s+\w+\s+(\w+)\s+(\w+)\s+([0-9]+)\s+/){
                        my $aa       = $1;
                        my $chain_   = $2;
                        my $pos      = $3;
                        $l =~ /^(....................)..(.*)/ ;
                        if ($chain_ eq $chain and $ini_pos == $pos) {print OUT "$1 $number$2\n"}
                        elsif ($chain_ eq $chain and $pos > $ini_pos and $pos <= $last_pos) {print OUT "$1 $number$2\n"}
                	}
		}

		print OUT "TERM\n";
	}
	if ($line =~ /^\s+[0-9]+-[0-9]+/){
		$number ++;
		my @LINE = split (" ",$line);
		my ($f1,$f2) = split ("-",$LINE[1]);
		my ($f3,$f4) = split ("-",$LINE[0]);

		system( "cat $pdb.*$f1-*$f2-$pdb.*$f3-*$f4*.pdb > temp");

		open TEMP, "temp";
		while (my $l = <TEMP>){
			if ($l =~ /^$/){}
			else{
				$l =~ /^(....................)..(.*)/ ;
				print OUT "$1 $number$2\n";
			}
		}
		print OUT "TERM\n";
		
	}
}

__END__
	if ($line =~ /([0-9]+)\.\s+/){
		my $nb = $1;
                chomp $line;
                my @LINE  = split (" ",$line); 
                my $query = $LINE[1];
                $query =~ s/-/./;
                my $target = $LINE[2];
                $target =~ s/-/./;
                my $A1 = "$pdb.$query.$pdb.$target.$chain.pdb";
		my $A2 = "$pdb.$target.$pdb.$query.$chain.pdb";
		my $outfile = "$pdb\_$chain\_$nb.pdb";
		open OUT, ">$outfile";
 		open A1, "$A1";
		while (my $la1 = <A1>){
			if ($la1 =~ /^$/){}
			else{
				
				$la1 =~ /^(....................)..(.*)/ ;
				print OUT "$1 1$2\n";
			}
		}
		close A1;
		print OUT "TERM\n";
		open A2, "$A2";
		while (my $la2 = <A2>){
			if ($la2 =~ /^$/){}
			else{
				$la2 =~ /^(....................)..(.*)/ ;
				print OUT "$1 2$2\n";
			}
		}
		close A2;
		print OUT "TERM\n";
#        	print "$pdb.$query.$pdb.$target.$chain.pdb\n";
        }
}
