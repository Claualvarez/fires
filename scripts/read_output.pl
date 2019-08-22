#!/usr/bin/perl

use strict ;
use Math::Complex;
use List::MoreUtils qw(uniq);
use POSIX;

my $in    ;
my $c ;
my $options = "";
my $chain_name ;

for (my $i = 0; $i < @ARGV ; $i++) {$options    .= " $ARGV[$i] " }
if ($options =~ / -i\s+(\S+) /)    {$in          = $1}
if ($options =~ / -c\s+(\S+) /)    {$c           = $1}
if ($options =~ / -ch\s+(\S+) /)   {$chain_name      = $1}

my %TAG =(
	0 => "A", 1 => "B", 2 => "C", 3 => "D", 4 => "E", 5 => "F"
);

my %CODE=(
        "A"=>"ALA", "R"=>"ARG", "N"=>"ASN", "D"=>"ASP",
        "C"=>"CYS", "E"=>"GLU", "Q"=>"GLN", "G"=>"GLY",
        "H"=>"HIS", "I"=>"ILE", "L"=>"LEU", "K"=>"LYS",
        "M"=>"MET", "F"=>"PHE", "P"=>"PRO", "S"=>"SER",
        "T"=>"THR", "W"=>"TRP", "Y"=>"TYR", "V"=>"VAL"
);


my (@NET,@INDEX,@SPANN, %QUERYS, @QUERYS,%RESULTS);

chomp $in;

my ($pdb,undef,undef) = split (/\./,$in);


open IN, "$in" or die "\n$!" ;
while (my $line = <IN>){
	if ($line =~ /^$/){next}
	my ($Q,$T,undef) = split (" ",$line,3);
	my $q = $Q;
	$q =~ s/[A-Z]//g;  #remove letters
	my ($q_i, $q_f) = split ("-",$q);
	my $t = $T;    
	$t =~ s/[A-Z]//g;
	my ($t_i, $t_f) = split ("-",$t);
	push @INDEX , ($q_f , $t_f);
	$NET[$q_i][$q_f] .= "$Q:$T!";
	$SPANN[$q_i][$q_f] .= "$q_i:$q_f!";
	$SPANN[$q_i][$q_f] .= "$t_i:$t_f!";
}

my @S = sort {$a <=> $b} @INDEX;


my ($candidate);
open SCORES, ">scores.temp";
for (my $x = 0 ; $x < $S[$#S] ; $x ++){
	for (my $y = 0 ; $y < $S[$#S] ; $y ++){
		 
		if ($NET[$x][$y] =~ /^$/){}
		else {
			my @FRAGMENT = split ("!",$SPANN[$x][$y]);
			my (@SUMA, $spann, @LAST) ; 
			my ($flag,@FRAGMENTS) = decide(@FRAGMENT);
			
			if ($flag == 1){
				foreach my $element (@FRAGMENTS){
					my ($ini,$fin) = split (":",$element);
					push @LAST , $fin ;
		
					for (my $i = $ini ; $i < $fin ; $i ++){$SUMA[$i] = 1}

					open CAND, "candidates.txt";
					while (my $ca = <CAND>){
						chomp $ca;
						my @CA = split (/\t/,$ca);
						$CA[0] =~ s/[A-Z]//g;
						$CA[0] =~ s/-/:/g;
						if($element =~ /$CA[0]/){$candidate = $CA[0]}
					}	
				}

				foreach my $j (keys @SUMA){$spann = $spann + $SUMA[$j]}
				my @FRAGMENT= uniq sort {$a <=>$b} @FRAGMENT;
				print SCORES "$candidate @FRAGMENT $spann\n";
				my @LASTO = sort {$a <=> $b} @LAST ;
				my $key = "$x:$y";
 				push @QUERYS, $key;
				$QUERYS{$key} = "$spann\t$NET[$x][$y]" ;
				
			}

			if ($flag == 2){
				my ($ini,$fin);
				foreach my $element (@FRAGMENTS){
					($ini,$fin) = split (":",$element);
					push @LAST , $fin ;
					open CAND, "candidates.txt";
					while (my $ca = <CAND>){
						chomp $ca;
						my @CA = split (/\t/,$ca);
						$CA[0] =~ s/[A-Z]//g;
						$CA[0] =~ s/-/:/g;
						if($element =~ /$CA[0]/){$candidate = $CA[0] }
					}
				}
				
				for (my $i = $ini ; $i < $fin ; $i ++){$SUMA[$i] = 1}
				foreach my $j (keys @SUMA){$spann = $spann + $SUMA[$j]}

				print SCORES "$candidate $spann\n";
				my @LASTO = sort {$a <=> $b} @LAST ;
				my $key = "$x:$y";
				push @QUERYS, $key;
				$QUERYS{$key} = "$spann\t$NET[$x][$y]" ;
			}
		}
	}
}

sub decide{
	my @FRAGMENTS = @_ ;
	my $last = $#FRAGMENTS ;
	$last = ($last + 1) ;  
	my $j = 0;
	my @OVERLAP ;
	for (my $i = 0 ; $i < $last ; $i = $i + 2){
		my $j = $i + 1;
		if ($FRAGMENTS[1] ne $FRAGMENTS[$j]){
			push @OVERLAP , $FRAGMENTS[1] ;
			push @OVERLAP , $FRAGMENTS[$j] ;
		}

	}
	@OVERLAP = sort {$a <=> $b} @OVERLAP;
	@OVERLAP = uniq @OVERLAP;
	$last = $last / 2;
	my $flag = 1;
	my @NONOVER = @FRAGMENTS;
	
	for (my $i= 1 ; $i < $last ; $i ++){
		
		my $k = $i-1 ;
		my ($i1,$f1) = split (":",$OVERLAP[$k]);
		my ($i2,$f2) = split (":",$OVERLAP[$i]);
		
		if ($i2 < $f1){
			$flag = 2;
			last;
		}
		if($i1 == $i2 and $f1 == $i2){
			@NONOVER = uniq @OVERLAP;
			$flag = 1;
		}
	}
	
	return ($flag,@NONOVER);
}


my $tag = 0;

@QUERYS=uniq@QUERYS;
foreach my $llave (keys @QUERYS){
	$RESULTS{$TAG{$tag}} .= "$QUERYS{$QUERYS[$llave]}\n";
	if ($llave == 0 or $llave == $#QUERYS){}
	else{
		my ($i,$f) = split (":",$QUERYS[$llave]);
		my ($I,$F) = split (":",$QUERYS[$llave+1]);
		$tag ++;
	}
	
}

if (scalar (keys %RESULTS) == 0){print "\tNo repeats found!";}
else{ 
	print "\n";
}


my @DISP;

my $nb = 0;
my %LINE;
foreach my $element (sort {$a<=>$b} keys %RESULTS){
	my %POINTER;
	my @RES = split (/\n/, $RESULTS{$element});
	my @ORES = sort {$a <=> $b} @RES;

	my (undef,$chain) = split (" ",$ORES[$#ORES],2);
	my @FRAGMENTS = split ("!", $chain);
	my $allfrag;
	foreach my $l (keys @FRAGMENTS){
		$FRAGMENTS[$l] =~ s/:/\t/;
		my ($tm,$rmsd,$matched) = score($FRAGMENTS[$l]);
		push @DISP, "$tm,$FRAGMENTS[$l],$rmsd,$matched";
	}

}


my @DISPLAY = sort {$b<=>$a} @DISP ;
print "\n";

open INF, "tab.temp";
open OOU, ">candidates.txt";
while (my $inf=<INF>){
	chomp $inf;
	open REV, "scores.temp";
	my $token = "no";
	my @INF = split (/\t/,$inf);
	while (my $rev = <REV>){
		chomp $rev;
		$rev =~ s/:/-/g;
		my @REV = split (" ",$rev);
		if ($INF[0] =~ /$REV[0]/){
			print OOU "$INF[0]\t$INF[1]\t$INF[2]\t$INF[3]\t$INF[4]\t";
			foreach my $ou (@REV[1..$#REV-1]){print OOU "$ou;";}
			print OOU "\t$REV[$#REV]\n";
			$token = "yes";
			last;
		}
	}
	if ($token eq "no"){print OOU "$inf\n"}
}

sub score{
	my ($cumulative,$ltarget,$d0,$rmsd,$matched) ;
	my $input = @_[0];
	my ($query,$target) = split (/\t/,$input);
	$query =~ s/-/./;
	$target	=~ s/-/./;
	$in =~ /(\S+)\.(\S+)\.temp/ ;

	my $pdb = $1;
	my $clique = "$pdb.$query.$pdb.$target.pdb.1.clique\n";
	my $first_raw = "$pdb.$target.$pdb.$query.1.pdb\n";
	my $second_raw = "$pdb.$query.$pdb.$target.1.pdb\n";
	
	open CLI, "$clique" or die "\n$!";
	while (my $line = <CLI>){
		my @LINE = split (" ",$line);
		my ($fx,$fy,$fz,$sx,$sy,$sz);
		if ($line =~ /The number of matched atoms =\s+(\S+)\s/){
			$matched = $1;
		}
		if ($line =~ /RMSD =\s+(\S+)/){$rmsd = $1}
		if ($line =~ /^Number.*=\s+([0-9]+)\s+and\s+([0-9]+)/){
			#$ltarget = $1 ;
			my $optionA = $1;
			my $optionB = $2;
			my @OPTIONS;
			push @OPTIONS , $optionA ;
			push @OPTIONS , $optionB ;
			my @SOPTIONS = sort {$a <=> $b} @OPTIONS;
			$ltarget =  $OPTIONS[0];
		}
		$d0    = 1.24 * (($ltarget - 15) **(1/3)) - 1.8;
		if ($line =~ /^$c/){
			open FIRST, "$first_raw" or die "\n$!";
			while (my $f = <FIRST>){
				my @F = split (" ",$f);
				if ($f =~ /CA\s+$CODE{$LINE[2]}\s+$c\s+$LINE[1]/) {$fx=$F[6]; $fy=$F[7]; $fz=$F[8];
				last}
			}
			open SECOND, "$second_raw" or die "\n$!";
			while (my $s = <SECOND>){
				my @S =	split (" ",$s);	
				if ($s =~ /CA\s+$CODE{$LINE[6]}\s+$c\s+$LINE[5]/) {$sx=$S[6]; $sy=$S[7]; $sz=$S[8]; 
				last}
			}
			my $di = (sqrt(($fx-$sx)**2 + ($fy-$sy)**2 + ($fz-$sz)**2)) ;
			my $chunk = 1 / (1+( ($di /$d0) ** 2)) ;
			$cumulative = $cumulative + $chunk ;
		}
	
	}
	if ($ltarget <= 0){return("err")}
	else{
		my $unround =  $cumulative * (1/$ltarget);
		my $tm = sprintf("%.4f", $unround);
	
		return($tm,$rmsd,$matched)
	}
}
