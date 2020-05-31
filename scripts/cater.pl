#!/usr/bin/perl
# … 
# Autor: Claudia Alvarez-Carreño    ccarreno6(at)gatech.edu
use strict ;

my $options  ;
my $infile   ;
my $stored   ;
my $chain   ;
my $outname;
my $reference ;
my $chain_name ;
my $pdb ;
my $outfile_name ;
my $path_to_fires;
my $path_to_click ;
my $range;
my $pairs ;

for (my $i = 0; $i < @ARGV ; $i++) {$options .= " $ARGV[$i] " }
if ($options =~ / -i\s+(\S+) /)    {$infile   = $1}
if ($options =~ / -c\s+(\S+) /)    {$chain    = $1}
if ($options =~ / -ch\s+(\S+) /)   {$chain_name    = $1}
if ($options =~ / -pdb\s+(\S+) /)  {$pdb    = $1}
if ($options =~ / -out\s+(\S+) /)  {$outname    = $1}
if ($options =~ / -path\s+(\S+) /) {$path_to_fires    = $1}
if ($options =~ / -click\s+(\S+) /) {$path_to_click    = $1}

chomp $infile;
my $dssp_name = "$pdb\_$chain_name.dssp";

open IN, "$infile" or die "$infile \n$!" ;
my ($number,$nb) = 0;
while (my $line = <IN> ){
	if ($line =~ /^([0-9]+)\..*Query: \[(.*)\].*:(.*)/){
		if(-f "temp_ali" ){
			print "\n  Structure-derived pairwise sequence alignment of Query [ $range ] and Pairs$pairs\n";
			system("cat temp_ali "); system ("rm temp_ali")
		}
		print "\n$line";
		$nb = $1;
		$number = 0;
		$range = $2;
		$pairs = $3;
		my ($ini_pos,$last_pos) = split ("-",$range);
		$outfile_name = "$outname\_$chain_name\_$nb.pdb";
		open OUT, ">$outfile_name";
		open REF, ">query.pdb";
		open HEADER, ">head.temp";

		open INPDB, "$pdb\_$chain_name.pdb" or die "\n$!";;
		while (my $l = <INPDB>){
				if($l =~ /^REMARK RES/){print OUT $l; print REF $l; print HEADER $l}
		        if($l =~/^ATOM\s+[0-9]+\s+\w+\s+(\w+)\s+(\w+)\s+([0-9]+)\s+/){
                        my $aa       = $1;
                        my $chain_   = $2;
                        my $pos      = $3;
                        $l =~ /^(....................)..(.*)/ ;
                        if ($chain_ eq $chain and $ini_pos == $pos) {print OUT "$1 $number$2\n";print REF "$1 $number$2\n" }
                        elsif ($chain_ eq $chain and $pos > $ini_pos and $pos <= $last_pos) {print OUT "$1 $number$2\n"; print REF "$1 $number$2\n"}
                	}
		}

		print OUT "TERM\n";
	}
	if ($line =~ /^\s+Query/){print "$line"}
	if ($line =~ /^\s+[0-9]+-[0-9]+/){
		print "$line";
		$number ++;
		my @LINE = split (" ",$line);
		my ($f1,$f2) = split ("-",$LINE[1]);
		my ($f3,$f4) = split ("-",$LINE[0]);
		system( "cat $pdb.*$f1-*$f2-$pdb.*$f3-*$f4*.pdb > temp.pdb");
		open TEMP, "temp.pdb";
		system("cat head.temp > target.pdb");
		open TARG, ">>target.pdb";

		while (my $l = <TEMP>){
			if ($l =~ /^$/){}
			else{
				$l =~ /^(....................)..(.*)/ ;
				print OUT "$1 $number$2\n";
				print TARG "$1 $number$2\n";
			}
		}
		print OUT "TERM\n";
		system("perl $path_to_fires/str-base-msa_v-1.pl query.pdb target.pdb -path $path_to_click -dssp1 $dssp_name -dssp2 $dssp_name >> temp_ali");
	}

}
if(-f "temp_ali" ){
	print "\n  Structure-derived pairwise sequence alignment of Query [ $range ] and Pairs$pairs\n";
	system("cat temp_ali ")
}

__END__
