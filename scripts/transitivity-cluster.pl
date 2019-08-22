#!/usr/bin/perl

use strict ;
use List::Util qw[min max] ;
use List::MoreUtils qw(uniq);
use POSIX;

my $in    ;
my $dssp  ;
my $verb  ;
my $chain ;
my $separation = 10;
my $gap        = 15;
my $options    = "";
my $ss         = 3 ;
my $overlay    = 3;
my $force      = 0;

for (my $i = 0; $i < @ARGV ; $i++)  {$options    .= " $ARGV[$i] " }
if ($options =~ / -i\s+(\S+) /)     {$in          = $1}
if ($options =~ / -add\s+(\S+) /)   {$verb        = $1}
if ($options =~ / -dssp\s+(\S+) /)  {$dssp        = $1}
if ($options =~ / -gap\s+(\S+) /)   {$gap         = $1}
if ($options =~ / -d\s+(\S+) /)     {$separation  = $1}
if ($options =~ / -ss\s+(\S+) /)    {$ss          = $1}
if ($options =~ / -c\s+(\S+) /)     {$chain       = $1}
if ($options =~ / -force\s+(\S+)\s?/) {$force       = $1}

my (%Y, %PAIRS, @CHAIN,@CENTERS,@ALLE);
my (@MATCENT,%HASHCENT,%ALLCENTER,%REPEATS,@CENTERSV);
open VER, "$verb" or die "\n$!" ;
my (@DISJOINT,$disjoint,@INF,@SUP);
while (my $v = <VER>){
	if ($v =~ /\w+\t\w/){
		chomp $v;
		push @DISJOINT, $v;
		my ($Xs , $Ys) = split (/\t/,$v);
		my ($x_i,$x_s) = split (":",$Xs);
		my (undef,$y_s) = split (":",$Ys);
		push @ALLE, $x_s;
		push @ALLE, $y_s;
		push @INF, $x_i;
		push @SUP, $x_s;
	}
}

close VER;
sort {$a <=> $b} @INF;
my $first = $INF[0];

sort {$a <=> $b} @SUP;
my $last_sup  = int @SUP;
my $last      = $SUP[$last_sup-1];
my $intervalo =  floor(($last - $first)/2);

sort @DISJOINT;
foreach my $element (@DISJOINT){$disjoint .= "$element\n"}

my @UNIQDISJOINT = dobles ($disjoint) ;
my @VERBRACHT    = verbringe(@UNIQDISJOINT);
my @PRELIMINARY  = sort {$a <=> $b} @VERBRACHT;

my ($step,$tagged) = 0;
my ($tag,@DOUBLES);
foreach my $v (keys @PRELIMINARY){
	for (my $w = 1 ; $w < $separation ; $w ++){
		if ($PRELIMINARY[$v] eq $PRELIMINARY[$v+$w]){
		$tag = "double";

	}else{
		my ($xi,$xf,$xc,$yi,$yf,$yc) = split (":", $PRELIMINARY[$v]);
		my ($xii,$xff,$xcc,$yii,$yff,$ycc) = split (":", $PRELIMINARY[$v+$w]);
		if ($xi <= $xii + $separation and $xi >= $xii - $separation){
			if ($xf <= $xff + $separation and $xf >= $xff - $separation){
				if ($yi <= $yii + $separation and $yi >= $yii - $separation){
					if ($yf <= $yff + $separation and $yf >= $yff - $separation){$tag = "double"}
					else{$tag = "single"}
				}else{$tag = "single";}
			}else{$tag = "single"}
		}else{$tag = "single"}
		
	}

	if ($tag eq "double"){
		chomp $PRELIMINARY[$v];
		chomp $PRELIMINARY[$v+$w];
		$DOUBLES[$v+$w] = "$PRELIMINARY[$v+$w]-f";
		$DOUBLES[$v]   = "$PRELIMINARY[$v]-i"; 
	}else{chomp $PRELIMINARY[$v];}
		$tagged = $v;
	}
}

my (@PRE_X);

for (my $double = 0 ; $double < $tagged ; $double ++){
	if ($DOUBLES[$double] =~ /^$/){
		my ($xi,$xf,$xc,$yi,$yf,$yc) = split (":", $PRELIMINARY[$double]);
		$PRE_X[$xc][$yc] = "$xi:$xf\t$yi:$yf";
		print "$xi:$xf\t$yi:$yf\n";
		$HASHCENT{$xc} .= "$xi:$xf\t";
		$HASHCENT{$yc} .= "$yi:$yf\t";
	}else{
		my ($element,$key) = split ("-",$DOUBLES[$double]);
		if ($key eq "f"){
			my ($xi,$xf,$xc,$yi,$yf,$yc) = split (":",$element);
			$PRE_X[$xc][$yc] = "$xi:$xf\t$yi:$yf"
		}
	}
}

if ($tagged > 1600){
	print STDERR "\nPlease wait\nFires will test > 1600 query-target pairs!\n";
}

my @SALLE = sort {$a <=> $b} @ALLE ;
for (my $i = 0 ; $i < $SALLE[$#SALLE] ; $i ++){
	for (my $j = 0 ; $j < $SALLE[$#SALLE] ; $j ++){
		if ($PRE_X[$i][$j] =~ /^$/){}
		else {
			if($PRE_X[$j][$i] =~ /^$/){}
			else{print "$PRE_X[$i][$j]\n"; $PRE_X[$j][$i] = ""}
		}
	}
}

sub verbringe{
	my (@RES1,$RES1, $distance_x_abs);
	my @DISJOINT = sort {$a <=> $b} @_;
	foreach my $d (@DISJOINT){
		my $A_centerc = center($d);
		my ($A_domino, $A_center) = split (/\|/,$A_centerc);
		foreach my $dd (@DISJOINT){
			my $B_centerc =  center($dd);
			my ($B_domino, $B_center) = split (/\|/,$B_centerc);
			my ($A_domino_x,$A_domino_y) = split (/\t/, $A_domino);
			my ($B_domino_x,$B_domino_y) = split (/\t/, $B_domino);
			my ($A_domino_xi,$A_domino_xs) = split (":", $A_domino_x);
			my ($B_domino_xi, $B_domino_xs) = split (":", $B_domino_x);
			my ($A_domino_yi,undef) = split (":", $A_domino_y);
			my (undef, $B_domino_ys) = split (":", $B_domino_y);
			if ($A_center eq $B_center){}
			else{
				if ($A_domino_xs < $B_domino_xi){$distance_x_abs = $B_domino_xi - $A_domino_xs}
				elsif($B_domino_xs < $A_domino_xi){$distance_x_abs = $A_domino_xi - $B_domino_xs}
				else {next}
				my ($x_Ac, $y_Ac) = split (":", $A_center);
				my ($x_Bc, $y_Bc) = split (":", $B_center);
				my $distance_x = $x_Ac - $x_Bc;
				my $distance_y = $y_Ac - $y_Bc;
				if ($distance_x_abs > $separation and $distance_x_abs < $intervalo){
					my $linf = abs($distance_y) - $gap;
					my $lsup = abs($distance_y) + $gap;
					if (abs($distance_x) > $linf and $distance_x < $lsup ){
						my @X_region = ($A_domino_xi, $B_domino_xs);
						my @X_R = sort {$a <=> $b} @X_region;
						my @Y_region = ($A_domino_yi, $B_domino_ys);
						my @Y_R = sort {$a <=> $b} @Y_region;			
						my $new_center_x = abs(floor($distance_x/2)) + $X_R[0] ;
						my $new_center_y = abs(floor($distance_y/2)) + $Y_R[0] ;
						if($new_center_x < $new_center_y){
							push @RES1, "$X_R[0]:$X_R[1]:$new_center_x:$Y_R[0]:$Y_R[1]:$new_center_y\n";
						}elsif($new_center_y < $new_center_x){
							push @RES1, "$X_R[0]:$X_R[1]:$new_center_x:$Y_R[0]:$Y_R[1]:$new_center_y\n";
						}
					}
				}
			}
		}
	}
	return (@RES1);
}

sub center{
	my $element = @_[0];
	my ($x_tips, $y_tips) = split (/\t/,$element);
	my ($x_inf, $x_sup)   = split (":", $x_tips);
	my ($y_inf, $y_sup)   = split (":", $y_tips);
	my $center_x = $x_inf + floor (($x_sup - $x_inf) / 2);
	my $center_y = $y_inf + floor (($y_sup - $y_inf) / 2);
	my $centered = "$element\|$center_x:$center_y";
	
	return ($centered);
}

open IN, "$in" or die "ERROR APERTURA \n$!" ;
while (my $l = <IN>){
	chomp $l;
	my ($x, $y, $ref) = split (" ", $l);
 	$Y{$x}           .= "$y-$ref\t";
}
close IN;

foreach my $x  (sort {$b <=> $a} keys %Y){
	my @Ys = split (" ",$Y{$x});
	foreach my $element (@Ys){
		my ($y, $ref) = split ("-",$element);
		my @Xs        = split (" ",$Y{$y});
		foreach my $e (@Xs){
			my ($x_de_y, $ref_2) = split ("-",$e);
			if ($x eq $x_de_y) {my $key = "$ref---$ref_2" ; $PAIRS{$key} ++ }
		}

	}
}

foreach my $k  (sort {$a <=> $b} keys %PAIRS){
	if ($PAIRS{$k} > $overlay){
		my $check = sse ($k,$dssp,$chain);
		open IN, "$in" or die "\n$!" ;
		while (my $l = <IN>){
			my ($x,$y,$refb) = split (" ",$l);
			my ($ref,$what)  = split ("---",$k); 
	
			if ($ref eq $refb){
				my ($inix,$finx) = split (":",$what);
				my ($iniy,$finy) = split (":",$refb);
				my $restax       = $finx - $inix;
				my $bischenx     = int ($restax/2);
				my $centrox      = int ($restax/2) + $inix ;
				my $restay       = $finy - $iniy;
                		my $centroy      = int ($restay/2) + $iniy ;
				my $bischeny     = int ($restay/2) + $centroy;

				
				if ($check eq "SI") {

					$CENTERS[$centroy] .= "$centrox\t";
					$CENTERS[$centrox] .= "$centroy\t";
					$HASHCENT{$centrox} .= "$what\t";
					$HASHCENT{$centroy} .= "$refb\t"; 
				}
			}
		}
	}
}

domino(@CENTERS);
domino(@CENTERSV);
sub domino{
	my $cont    = 0;
	my @CENTERS = @_;
	foreach my $element (sort keys @CENTERS){
		$cont ++;
		my ($acumula,@NEWCHUNK,$thrd);
		for (my $x = 0 ; $x < $cont ; $x ++){
			if ($CENTERS[$x] =~ /^$/){}
			else{
			
				my @CENT      = mono ($CENTERS[$x]); 
				my $agregados = int @CENT;
				$acumula      = $acumula + $agregados;
				for ($thrd = 0 ; $thrd < $agregados ; $thrd ++){
					@NEWCHUNK = mono($CENTERS[$CENT[$thrd]]);
					sort @NEWCHUNK;
					my @ACUMULADO = @CENT;
		 			push @ACUMULADO, @NEWCHUNK;
					sort @ACUMULADO;
					my $next = $CENT[$thrd] + 0; #ANTES +1
					if ($CENTERS[$next] =~ /^$/){}
					else {push @ACUMULADO, $next;}
					my $prev = $CENT[$thrd] - 0; #ANTES -1
					if ($CENTERS[$prev] =~ /^$/){}
					else {push @ACUMULADO, $prev;}
					@CENT      = uniq (@ACUMULADO);
					$agregados = int  @CENT;
				}

				my @OCENT           = sort {$a <=> $b} @CENT;
				my $orederedCenter  = join " ", @OCENT;
				$REPEATS{$OCENT[0]} = $orederedCenter;
			}
		}
	}

}

my ($res);

foreach my $element (sort {$a <=> $b} keys %REPEATS){
	my @CENTROS = split (" ",$REPEATS{$element});
	my $uno     = edita ($HASHCENT{$element});
	foreach my $c (@CENTROS){
		my $pair = edita($HASHCENT{$c});
		my ($xp,$yp) = split (":",$pair);
		my ($xu,$yu) = split (":",$uno);
		if ($uno eq $pair){}
		elsif($xp eq $yp){}
		elsif(($xu < $xp and $yu > $yp) or ($xu > $xp and $yu < $yp)){}
		else {$res .= "$uno\t$pair\n"}
	}
}


my @FINAL = dobles($res);

foreach my $fin (@FINAL){print "$fin\n"}
if ($#FINAL < 20 and $force eq "on"){system ("cat $verb") }


sub dobles{
	my $raw = @_[0];
	my @RA = split (/\n/,$raw);
	my @RAW = sort {$a <=> $b} @RA;
	foreach my $dom_f (keys @RAW){
		my ($A, $B) = split (" ",$RAW[$dom_f]);
		my ($AA, $BB) = split (" ",$RAW[$dom_f+1]);
		my $spot = $dom_f+1;
		if ($RAW[$dom_f+1] =~ /^$/){}
		else{
			my ($ini_A,$fin_A)= split (":",$A);
			my ($ini_AA,$fin_AA)= split (":",$AA);
			my ($ini_B,$fin_B)= split (":",$B);
			my ($ini_BB,$fin_BB)= split (":",$BB);
			if ($ini_A == $ini_B or $ini_A == $ini_B-1){splice @RAW , $dom_f , 1}
			if ($A eq $AA or $ini_A == $ini_AA -1  or $ini_A == $ini_AA +1){
				if ($ini_B == $ini_BB or ($ini_B == $ini_BB -1) or ($ini_B == $ini_BB +1)){
					if ($fin_B == $fin_BB or ($fin_B == $fin_BB -1) or ($fin_B == $fin_BB +1)){splice @RAW , $spot , 1}
				}
			}
		}
	}
	return (@RAW);
}


sub edita{
	my $AB    = @_[0];
	my @TODOS = split (/\n/,$AB);
	sort {$a <=> $b} @TODOS;
	my (@A,@B);
	foreach my $dom (@TODOS){
		my ($a,$br)   = split (":",$dom);
		my ($b,undef) = split (" ",$br);
		push @A, $a ;
		push @B, $b ;
	}

	my @SORTEDA = sort {$a <=> $b} @A;
	my @SORTEDB = sort {$a <=> $b} @B;
	my $ultimo  = int @SORTEDB;
	my $domino  = "$SORTEDA[0]:$SORTEDB[$ultimo-1]";

	return ($domino);	

}

sub mono{
	my $saver;
	my $arreglo = @_[0];
	my @REP     = split (" ",$arreglo);
	sort {$a   <=> $b}@REP;
	my @NOREP   = uniq @REP ;
	return (@NOREP);
}


sub sse{
	my $k            = @_[0];
	my ($A,$B)       = split ("---",$k);
	my ($iniA,$finA) = split (":",$A);
	my ($iniB,$finB) = split (":",$B);
	$iniA            = $iniA -1 ;
	my $dssp         = @_[1];
	my ($initial,@SS);
	my $stored_str   = "FIRST";
	my $stored_chain = @_[2];
	my $count        = 0;
	my ($H, $B, $T);

	open IN, "$dssp" or die "\n$!" ;
	while (my $line = <IN>){
        	chomp $line;

        	if ($line =~/\#/) {$initial = "INICIO"}
        	if ($initial eq "INICIO" and $line =~/^\s+[0-9]+/){ 
                	$line      =~/.(....).(....).(.).(.)..(.)/;
                	my $nb     = $1; 
                	my $res_nb = $2 + 0;
			if ($res_nb >= $iniA and $res_nb <= $finA){
                		my $chain  = $3;
                		my $aa     = $4; $aa=~s/[a-z]/C/g;
                		my $str    = $5; 
				$str       =~ s/ /-/;
				$str       =~ s/[GI]/H/;
				$str       =~ s/E/B/;
				$str       =~ s/S/T/;
                		$res_nb    =~s/ //g;
                		if ($line =~/!/) {next}
                		
                		if ($stored_chain eq $chain and $stored_str ne $str){
                        		$stored_str = $str;					
					if ($str eq "H") {$H ++}
					if ($str eq "B") {$B ++}
					if ($str eq "T") {$T ++}
					if ($str eq "-") {$H  = 0; $B=0; $T = 0}
                		}
				
				elsif ($str eq $stored_str and $stored_chain eq $chain){	
					if    ($str eq "H" and $H == 1) {$count ++}
					elsif ($str eq "B" and $B == 1) {$count ++}
					elsif ($str eq "T" and $T == 1) {$count ++}
					$H  = 0; $B=0; $T = 0 ;
				}

                		if ($stored_chain ne $chain) {next}
			}
        	}
	}

	close IN;

	if($count >= $ss){ return("SI") }
	else{ return("NO") }

}
