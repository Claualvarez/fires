#!/usr/bin/perl
# … 
# Autor: Claudia Alvarez-Carreño

use strict ;
use POSIX;
use List::MoreUtils qw(uniq);

my $options  ;
my $infile   ;
my $chain    ;

for (my $i = 0; $i < @ARGV ; $i++) {$options .= " $ARGV[$i] " }
if ($options =~ / -i\s+(\S+) /)    {$infile   = $1}
if ($options =~ / -c\s+(\S+) /)    {$chain    = $1}
chomp $infile ;

open IN, "$infile" or die "ERROR APERTURA \n$!" ;

my ($display,%KEEPER,@KEPT,%CLUSTER,%REFERENCE,%FINAL,$referred,%PAIRS,@SUMMARY,@RESULTS,%PAIRS_INFO,%DISP_INFO);

while (my $line = <IN> ){
	if ($line =~ /\S+/){
		chomp $line;
		#print "$line\n";
		my ($ref,$cplt,$info) = split (/\t/,$line,3);
		my @INF               = split (/\t/,$info);

		if ($INF[2] >= 0.5){
			my $score                  = $INF[0] * $INF[2];
			#my $score                 = $INF[0] / $INF[2] ;
			$PAIRS{"$ref\t$cplt"}      = $score;
			$PAIRS_INFO{"$ref\t$cplt"} = "$INF[0]\t$INF[1]\t$INF[2]" ;
			if ($info =~ /warning/ ){}
			else                    {$KEEPER{$ref} .= "$cplt\t";}
		}
	}
}
close IN;

for (my $i = 0 ; $i < 20000 ; $i ++){
	for (my $j = 0 ; $j < 20000 ; $j ++){
		my $check_1 = "$i-$j";
		if ($KEEPER{$check_1} =~ /\S+/){			
			my $flag             = pointer($check_1);
			my @ENDING           = split (/\t/, $KEEPER{$check_1});
			my @LINE             = ($check_1,@ENDING);
			$REFERENCE{$check_1} = $check_1 ;
		}
	}
}

foreach my $element (sort {$a <=> $b} keys %CLUSTER){
	my @ELEMENTS = split (/\t/, $CLUSTER{$element});
	@ELEMENTS    = uniq @ELEMENTS;
	@ELEMENTS    = sort {$a<=>$b} @ELEMENTS;
	my @FINAL    = split(/\t/,$CLUSTER{$element});
	@FINAL       = uniq @FINAL;
	
	foreach my $final (@FINAL){
		if ($element ne $final ){$FINAL{$element} .= "$final\t"  }
		else                    {$FINAL{$element} .= "$element\t"}
	}

}

my ($linef,@LINESF,%FIRST_POSITION);
foreach my $f (sort {$a<=>$b} keys %FINAL){

	my @INI = split (/\t/,$FINAL{$f});
        my (@FIN,$iINI,$iFIN) ;

        foreach my $i (@INI){ my @ENDS = split (/\t/,$KEEPER{$i}) ; push @FIN, @ENDS }

	unshift @INI, $f ;
	
	$referred = "";

	my (@iINI,@iFIN);
	foreach my $i (@INI){
		my ($ini,$fin) = split (/-/,$i);
		push @iINI, $ini ;
		push @iFIN, $fin ;
		#$iINI          = $iINI + $ini;
		#$iFIN          = $iFIN + $fin;
		$referred      = "$referred\t$i";
	}	
	
	@iINI = sort {$a<=>$b} @iINI;
	@iFIN = sort {$a<=>$b} @iFIN;
	#$iINI              = $iINI / ($#INI + 1);
	#$iINI              = floor $iINI;
	#$iFIN              = $iFIN / ($#INI + 1);
	#$iFIN              = ceil $iFIN;
	$linef             = "$iINI[0]-$iFIN[$#iFIN]";
	$REFERENCE{$linef} = $referred;
	@FIN               = uniq @FIN;
	@FIN               = sort {$a<=>$b} @FIN;
	
	foreach my $elem (@FIN){$linef="$linef\t$elem"}
	
	push @LINESF, $linef;
}

foreach my $line (sort {$a <=> $b} @LINESF){
	my @LINE = split (/\t/,$line);
	my $disp;

	if ($#LINE == 1){ $disp = findScore (@LINE) ; push @RESULTS, $disp }
	else{entscheiden (@LINE)}
}

#@RESULTS = sort {$a <=> $b} @RESULTS;

sub doublecheck {
	my @summary = "";
	my @ELEMS = @_;
	my $flag ;
	foreach my $elem (@ELEMS){
                my ($ini,$fin) = split ("-",$elem);
                for (my $i=$ini+1; $i<$fin+1; $i++) { $summary[$i] ++ }
        }
	foreach my $elem (keys @ELEMS){

                my ($ref1a,$ref2a) = split("-",$ELEMS[$elem]);
                my $mid_a = floor (($ref1a + $ref2a) / 2);
                if ($summary[$ref1a]>=2 or $summary[$ref2a] >= 2 or $summary[$mid_a]>=2){
                        $flag="internal_problem";
                }
        }
	return ($flag);
}

my ($rank,$overall_rank);
if ($#RESULTS >= 0){open OUTTAB, ">$infile.out"}

else{print "\tNo repeats found!\n";}
my $message ;
foreach my $results (sort {$b <=> $a} @RESULTS){
	my $flag = "pass";
	$results =~ /Pair\(s\): \[ (.*) \]/;
	my @ELEM = split ("; ",$1);
	foreach my $elem (@ELEM){
		my ($ini,$fin) = split ("-",$elem);
		
		$message = doublecheck (@ELEM);
	
		for (my $i=$ini+1; $i<$fin+1; $i++) {$SUMMARY[$i] ++}
	}
	if ($message eq "internal_problem"){@SUMMARY = ""; $flag ="no"}
	foreach my $elem (keys @ELEM){
		#print "$results\n";
		my ($ref1a,$ref2a) = split("-",$ELEM[$elem]);
		my $mid_a = floor (($ref1a + $ref2a) / 2);
		if ($SUMMARY[$ref1a]>=2 or $SUMMARY[$ref2a] >= 2 or $SUMMARY[$mid_a]>=2 ){
			$flag="no";
		}
	}

	if ($flag eq "pass"){
		$rank ++ ; 
		my @DISP_INFO = split (/\n/,$DISP_INFO{$results});
		#$DISP_INFO{$results} =~ /.*Pair\(s\): \[ (\S+) \]/
		#print "";
		@DISP_INFO = uniq @DISP_INFO ;
		@DISP_INFO = sort {$a <=> $b} @DISP_INFO ;
		print "\n$rank.\tGlobalScore= $results\n\n" ;
		printf '%20s', "Query";
		printf '%10s', "Target";
		printf '%10s', "AlnRes";
		printf '%10s', "RMSD";
		printf '%10s', "TM-score";		
		print "\n";
		foreach my $in (@DISP_INFO){
			my @IN = split (" ",$in);
				printf '%20s', "$IN[1]";
				printf '%10s', "$IN[0]";
				printf '%10s', "$IN[2]";
				printf '%10s', "$IN[3]";
				printf '%10s', "$IN[4]";
				print "\n";
			}
		}
        $overall_rank ++ ;
        my @DISP_INFO = split (/\n/,$DISP_INFO{$results});
        @DISP_INFO = uniq @DISP_INFO ;
        @DISP_INFO = sort {$a <=> $b} @DISP_INFO ;
        print  OUTTAB "\n$overall_rank.\tGlobalScore= $results\n\n" ;
        printf OUTTAB '%20s', "Query";
        printf OUTTAB '%10s', "Target";
        printf OUTTAB '%10s', "AlnRes";
        printf OUTTAB '%10s', "RMSD";
        printf OUTTAB '%10s', "TM-score";
        print  OUTTAB "\n";
        foreach my $in (@DISP_INFO){
                my @IN = split (" ",$in);
                printf OUTTAB '%20s', "$IN[1]";
                printf OUTTAB '%10s', "$IN[0]";
                printf OUTTAB '%10s', "$IN[2]";
                printf OUTTAB '%10s', "$IN[3]";
                printf OUTTAB '%10s', "$IN[4]";
                print  OUTTAB "\n"

	}
	$overall_rank ++ ;
	my @DISP_INFO = split (/\n/,$DISP_INFO{$results});
	@DISP_INFO = uniq @DISP_INFO ;
	@DISP_INFO = sort {$a <=> $b} @DISP_INFO ;
	print  OUTTAB "\n$overall_rank.\tGlobalScore= $results\n\n" ;
	printf OUTTAB '%20s', "Query";
        printf OUTTAB '%10s', "Target";
        printf OUTTAB '%10s', "AlnRes";
        printf OUTTAB '%10s', "RMSD";
        printf OUTTAB '%10s', "TM-score";
        print  OUTTAB "\n";
        foreach my $in (@DISP_INFO){
        	my @IN = split (" ",$in);
        	printf OUTTAB '%20s', "$IN[1]";
                printf OUTTAB '%10s', "$IN[0]";
                printf OUTTAB '%10s', "$IN[2]";
                printf OUTTAB '%10s', "$IN[3]";
                printf OUTTAB '%10s', "$IN[4]";
                print  OUTTAB "\n"

	}
}
print "\n";


sub findScore{
	my @LINE = @_;
	my ($line,%DISP,$disp,$megaline);
	$line = "Query: [ $LINE[0] ]\tPair(s): [ " ;
	for (my $l = 1 ; $l < $#LINE+1 ; $l ++ ){ $line = "$line$LINE[$l]; " }

	my $last = $#LINE +1;
	my @S    = split (/\t/,$LINE[0]);

	foreach my $s (@LINE){
		my @INIS = split (/\t/,$REFERENCE{$LINE[0]});
		foreach my $ini (@INIS) {
			my @FINS = split (/\t/,$KEEPER{$ini});
			foreach my $fin (@FINS){
				foreach my $ini (@INIS) {
					for (my $k = 1 ; $k < $last ; $k ++){
						my $p = "$ini\t$fin";

						if ($fin eq $LINE[$k] and $PAIRS{$p} =~ /\S+/){
							my $key = "$LINE[0]\t$fin" ;
							$DISP{$key} = "$PAIRS{$p}";
							$megaline .= "$fin\t$ini\t$PAIRS_INFO{$p}\n";
							#last;
							#print "$LINE[0]--->>>$p $PAIRS_INFO{$p}\n";
						}
					}
					
				}
			}
		}
	}

	my ($cumulate,$prom);
	for (my $k = 0  ; $k < $#LINE ; $k ++){
		my $key    = "$LINE[0]\t$LINE[$k + 1]";
		$cumulate += $DISP{$key};
		$prom      = sprintf("%.3f",$cumulate) ;
	}

	$disp = "$prom\t$line]";
	$disp =~ s/; ]/ ]/;
	$DISP_INFO{$disp} = $megaline;
	return ($disp);
}

sub entscheiden{
	my @LINE     = @_;
	my @FRAGMENTS;
	my @ELEMENTS = split(/\t/,$LINE[0]);
	my @E        = split (/\t/,$REFERENCE{$ELEMENTS[0]});

	foreach my $x (@E) {
		my @F = split (/\t/,$KEEPER{$x});
		foreach my $f (keys @F){
			my $key = "$x\t$F[$f]";
			push @FRAGMENTS, $key;
		}
	}

	decide ($ELEMENTS[0],@FRAGMENTS);
}

sub decide {
	my $firstElem = shift @_;
	my @FRAGMENTS = @_ ;
	my (@ELEMS1, @ELEMS2,@SCORES);
	foreach my $frag (@FRAGMENTS){
		my ($elem1,$elem2) = split (/\t/,$frag);
		push @ELEMS1, $elem1;
		push @ELEMS2, $elem2;
	}
	@ELEMS2 = sort {$a <=> $b} @ELEMS2;
	@ELEMS2 = uniq @ELEMS2 ;

	for (my $e = 0 ; $e < $#ELEMS2; $e ++){
		my ($ini,$fin)           = split (/-/,$ELEMS2[$e]);
		my ($ini_next,$fin_next) = split (/-/,$ELEMS2[$e+1]);
		if ($ini_next < $fin){
			my ($score, $scoree);
			foreach my $e1 (@ELEMS1){
				my $paire  = "$e1\t$ELEMS2[$e]";
				my $pairee = "$e1\t$ELEMS2[$e+1]";
				if ($PAIRS{$paire}  =~ /\S+/) {$score  = $PAIRS{$paire}}
				if ($PAIRS{$pairee} =~ /\S+/) {$scoree = $PAIRS{$pairee}}
			}
			if ($score > $scoree){
				my $pos = $e + 1;
				splice @ELEMS2, $pos, 1;
				$e -- ;
			}
			if ($score < $scoree){
				my $pos = $e ;
				splice @ELEMS2, $pos, 1;
				$e --;
			}
		}else{
			#
		}
	}
	unshift @ELEMS2, $firstElem;

	my $line = findScore (@ELEMS2);
	push @RESULTS, $line;
	#print "$line\n";
}


sub pointer{
	my $point      = @_[0];
	my (@POINTED);
	my $flag       = 1 ;
	my ($ini,$fin) = split ("-",$point);

	for (my $i = 1 ; $i < 8 ; $i ++){
		for (my $j = 1 ; $j < 8 ; $j ++){

			my $pace_fin   = $j+$fin;
			my $check_fin  = "$ini-$pace_fin";
			my $pace_ini   = $i+$ini;
			my $check_ini  = "$pace_ini-$fin";
			my $check_both = "$pace_ini-$pace_fin";

			if ($KEEPER{$check_ini}  =~ /\S+/) {push @POINTED , $check_ini  }
			if ($KEEPER{$check_fin}  =~ /\S+/) {push @POINTED , $check_fin  }
			if ($KEEPER{$check_both} =~ /\S+/) {push @POINTED , $check_both }

		}
	}

	@POINTED = uniq @POINTED;
	my ($mit);

	foreach my $clus (sort {$a <=> $b} keys %CLUSTER){
		
		my @GROUP = split (/\t/,$CLUSTER{$clus});
		
		if ($#GROUP == 0){ $flag = 0 }

		foreach my $element (@GROUP){ foreach my $po (@POINTED) { if ($element eq $po ){ $point = $clus } } }
	}

	$CLUSTER{$point} .= "$point\t";
	
	foreach my $element (@POINTED){$mit .= "$element\t";}
	$CLUSTER{$point} .= $mit;

	return ($flag);
}
