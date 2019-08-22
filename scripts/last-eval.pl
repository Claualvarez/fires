#!/usr/bin/perl

use strict ;
use List::Util qw[min max] ;
use List::MoreUtils qw(uniq);
use POSIX;

my $chain    ;
my $overlap  = 75  ;
my $coverage = 75  ;
my $rmsd     = 3 ;
my $options  = ""  ;
my $min_aa   = 15   ;
my $mma      = 15;
my $infile   ;
my $chaini ;


for (my $i = 0; $i < @ARGV ; $i++)    {$options .= " $ARGV[$i] " }
if ($options =~ / -c\s+(\S+) /)       {$chain    = $1}
if ($options =~ / -i\s+(\S+) /)       {$infile   = $1}
if ($options =~ / -o\s+(\S+) /)       {$overlap  = $1}
if ($options =~ / -n\s+(\S+) /)       {$coverage = $1}
if ($options =~ / -rmsd\s+(\S+) /)    {$rmsd = $1}
if ($options =~ / -ma\s+(\S+) /)      {$min_aa = $1}
if ($options =~ / -mma\s+(\S+) /)     {$mma = $1}
if ($options =~ / -ch\s+(\S+) /)      {$chaini = $1}

open IN, "$infile" or die "\n$!" ;

my (%FILES, $token,$pdb,$r);

my %CODE=(
        "A"=>"ALA", "R"=>"ARG", "N"=>"ASN", "D"=>"ASP",
        "C"=>"CYS", "E"=>"GLU", "Q"=>"GLN", "G"=>"GLY",
        "H"=>"HIS", "I"=>"ILE", "L"=>"LEU", "K"=>"LYS",
        "M"=>"MET", "F"=>"PHE", "P"=>"PRO", "S"=>"SER",
        "T"=>"THR", "W"=>"TRP", "Y"=>"TYR", "V"=>"VAL"
);



while (my $lines = <IN>){
	chomp $lines;
	$lines =~ /(.*).pdb\s+(\S+).pdb/;
	my $arg1       = $1;
	my $arg2       = $2;
	my (undef,$arg1s) = split (/\./,$arg1);
	my (undef,$arg2s) = split (/\./,$arg2);
	$arg1s =~ s/[A-Z]//g;
	$arg2s =~ s/[A-Z]//g;
	my ($ini1,$fin1)       = split ("-",$arg1s);
	my ($ini2,$fin2)       = split ("-",$arg2s);
	my $resta1 = $fin1-$ini1;
	my $resta2	= $fin2-$ini2;
	my $resta3      = $fin1-$ini2;
	if ($resta3 > 0){next}
	if ($fin1 >= $fin2 or $resta1 < $min_aa or $resta2 < $min_aa){next}
	my $cliquefile = "$arg1-$arg2\.pdb\.1\.clique\n"; 
	$cliquefile    =~ /(\w+)\.(\w+-\w+)-\w+\.\w+-\w+\.pdb\.1\.clique/;
	$pdb           = $1;
	my $iniciof    = $2;
	$FILES{$iniciof} .= "$cliquefile\n" ;
}


open CAND, ">candidates.txt";
foreach my $files (sort {$a<=>$b} keys %FILES){
	my $ind = 0;
	my @EACH = split (/\n/,$FILES{$files});
	foreach my $file (@EACH){
		my (@XCHAIN,@YCHAIN,@XAA,@SXAA,$tag);
		
		if ($file =~ /^$/){next}
		
		else{
			open CLIQUE, "$file" or  next; 			
        		while (my $l = <CLIQUE>){
				chomp $l;
				if ($l =~ /RMSD /){
					my (undef, $rmsdi) = split ("= ",$l);
					my ($tmi,undef,undef)=score($file,$chain,$chaini);
					print CAND "$rmsdi\t$tmi";
					if ($rmsdi > $rmsd){last}
				}
				if($l =~ /The number of matched atoms =/){	
					my (undef,$mma1) = split ("= ",$l);
					if ($mma1 < $mma){last}
					my @info = split(/\./,$file) ;
					my @in = split (/-/,$info[1]);
					print CAND "\n$in[0]-$in[1]\t$info[2]\t$mma1\t";
				}
				elsif ($l =~ /Structure Overlap/){
					my (undef,$score1) = split ("= ",$l);
					if ($score1 < $overlap){print CAND "\twarning: coverage < 75%"; last}
				}elsif ($l =~ /Number of/){
					
					my (undef,$variables) = split ("= ",$l);
					$variables =~ s/ //g;
					my @VARIABLES = split ("and",$variables);
					my @V = sort {$a <=> $b} @VARIABLES;
					my $score2 = floor($V[0]*100/$V[1]);
					if($score2 < $coverage){print CAND "\twarning: coverage < 75%"; last}
				}else{
					if ($l=~ /^\S? /){
						my @DATA = split (" ",$l);
						$XCHAIN[$DATA[1]] = $DATA[2];
						$YCHAIN[$DATA[1]] = $DATA[6];	
						$tag=  "yes";	
					}
				}
        		}
			
		}
		if ($tag eq "yes"){
			$ind ++;
			$file    =~ /(\w+)\.(\w+-\w+)-\w+\.(\w+-\w+)\.pdb\.1\.clique/;
			my $inname1  = "$1\.$2-$1\.$3\.1\.pdb";
			my $inname2  = "$1\.$3-$1\.$2\.1\.pdb";
			my $pdbname  = $1;
			my $query    = $2;
			my $target   = $3;
			print "$query\t$target\t";
			$query =~ s/-/\./g ;
			$target =~ s/-/\./g ;
			my $outname1 = "$pdbname\.$query.$pdbname\.$target\.1\.pdb";
			my $outname2 = "$pdbname\.$target.$pdbname\.$query\.1\.pdb"; 
			my $first;
			system ("cp $inname1 $outname1");
			system ("cp $inname2 $outname2");
			my $newfile = $file;
			$newfile =~ s/-/\./g;
			system ("cp $file $newfile");
			foreach my $element(keys@XCHAIN){
				if ($YCHAIN[$element]=~ /^$/){}
				else{$first = $element;last}
			}
			for (my $element = $first ; $element < $#XCHAIN; $element ++){
				if ($YCHAIN[$element]=~ /^$/){print "-";}
				else{print"$YCHAIN[$element]"}
			}
		}
	}
	
}

sub score{
        my ($cumulative,$ltarget,$d0,$rmsd,$matched,$c,$in) ;
        my $input  = @_[0];
	my $c      = @_[1];
	my @INPUT  = split (/\./,$input);
	my $query  = "$INPUT[1]";
	my $target = "$INPUT[2]";
	my $pdb    = "$INPUT[0]";
	$query =~ s/-$pdb//;
	$target =~ s/-$pdb//; 
	
        my $chaini = @_[3];

        my $clique = "$pdb.$query-$pdb.$target.pdb.1.clique\n";
        my $first_raw = "$pdb.$target-$pdb.$query.1.pdb\n";
        my $second_raw = "$pdb.$query-$pdb.$target.1.pdb\n";

        open CLI, "$clique" or die "\n$!";
        while (my $line = <CLI>){
                my @LINE = split (" ",$line);
                my ($fx,$fy,$fz,$sx,$sy,$sz);
                if ($line =~ /The number of matched atoms =\s+(\S+)\s/){
                        $matched = $1;
                }
                if ($line =~ /RMSD =\s+(\S+)/){$rmsd = $1}
                if ($line =~ /^Number.*=\s+([0-9]+)\s+and\s+([0-9]+)/){
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
                                if ($f =~ /CA\s+$CODE{$LINE[2]}\s+$c\s+$LINE[1]/) {$fx=$F[6]; $fy=$F[7]; $fz=$F[8]; last}
                        }
                        open SECOND, "$second_raw"or die "\n$!";
                        while (my $s = <SECOND>){
                                my @S = split (" ",$s);
                                if ($s =~ /CA\s+$CODE{$LINE[6]}\s+$c\s+$LINE[5]/) {$sx=$S[6]; $sy=$S[7]; $sz=$S[8]; last}
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


