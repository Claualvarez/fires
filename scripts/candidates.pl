#!/usr/bin/perl
# Author: Claudia Alvarez-Carreno   ccarreno6(at)gatech.edu

use strict ;

my $infile  ;
my $inref   ;
my $pdbcode ;
my $chain   ;

my $options = "";

for (my $i = 0; $i < @ARGV ; $i++) {$options .= " $ARGV[$i] " }
if ($options =~ / -c\s+(\S+) /)    {$chain     = $1}
if ($options =~ / -i\s+(\S+) /)    {$infile    = $1}
if ($options =~ / -dssp\s+(\S+) /) {$inref     = $1}
if ($options =~ / -pdb\s+(\S+) /)  {$pdbcode   = $1}


chomp $infile;
open IN, "$infile" or die "\n$!" ;

my ($name,@PIECES);

while (my $line = <IN>){
	my $nonoverlapping = filter($line);
	push @PIECES, $nonoverlapping;	
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
	foreach my $index (keys @INDEXES){
		my ($ini,$fin) = split (":",$INDEXES[$index]);
		cut($ini,$fin,$index,$chain,$pdbcode);		
	}
}

my $name ;

sub filter{
	my $line = @_[0];
	my ($query,$pair) = split (/\t/,$line);
	my ($query_i,$query_f) = split (":",$query);
	my ($pair_i,$pair_f) = split (":",$pair);
	if($query_f > $pair_i){return()}
	else{return($line)}
}

sub cut{
	my $ini   = @_[0];
	my $fin   = @_[1];
	my $index = @_[2];
	my $chain = @_[3];
	my $pdbcode = @_[4];
	open REF, "$inref" or die "\n$!" ;
	my $initial;
	while (my $ref = <REF>){
		chomp $ref;
		if ($ref =~/\#/) {$initial ="INICIO"}
		if (($initial eq "INICIO") and ($ref =~/^\s+[0-9]+/)){
			$ref =~/.(....).(....)...(.)..(.)/;
			my $nb     = $1; #$nb     =~s/ //g;
			my $res_nb = $2; 
			my $aa     = $3; $aa=~s/[a-z]/C/g;
			my $str    = $4; $str=~s/ /-/;$str=~s/[GI]/H/;$str=~s/E/B/;$str=~s/S/T/; 
			$res_nb    =~s/ //g; 
			if ($ref =~/!/) {next}
			if   ($res_nb eq $ini) {$name  = "$index\t$pdbcode\t$chain\t$CODE{$aa}\t$res_nb\t";}
			elsif($res_nb eq $fin) {$name .= "$chain\t$CODE{$aa}\t$res_nb";}
		}
	}
	close REF;
	corta($name);
	
}

my $cortado ;
sub corta{
	my $param = @_[0];
	my ($index, $pdb, $chain1, $aa1, $res_nb1, $chain2, $aa2, $res_nb2) = split (/\t/,$param);

	open INPDB,"$pdb.pdb" or die "\n$!" ;
	if   ($index == 0){print "$pdb\.$aa1$res_nb1-$aa2$res_nb2\.pdb\t"}
	elsif($index == 1){print "$pdb\.$aa1$res_nb1-$aa2$res_nb2\.pdb\n"}

	open OUT, ">$pdb\.$aa1$res_nb1-$aa2$res_nb2\.pdb" ;
	while (my $line = <INPDB>){
		if($line =~ /^REMARK RES/){print OUT $line}
		if ($line =~/^ATOM\s+[0-9]+\s+\w+\s+(\w+)\s+(\w+)\s+([0-9]+)\s+/){
			my $aa = $1; my $chain =$2; my $pos= $3;
			$line =~/^ATOM\s+[0-9]+(\s+\w+\s+\w+\s+\w+\s+[0-9]+\s+.*)/;
			my $columns = $1;
			if(($chain1 eq $chain) and ($aa1 eq $aa) and ($res_nb1 == $pos)){
				$cortado ++;
                                print OUT "$line";
                        }elsif(($chain2 eq $chain) and ($pos > $res_nb1)and ($pos <= $res_nb2)){
				$cortado ++;
                                print OUT "$line";
                        }


                }

	}
	close INPDB;
}

__END__


