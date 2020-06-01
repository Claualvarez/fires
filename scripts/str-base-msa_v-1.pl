use strict;
use List::MoreUtils qw(uniq);
use List::Util qw[min max] ;

my $options ;
my $inclick ;
my $lower_pdb ;
my $upper_pdb ;
my $pdb_1 ;
my $pdb_2 ;
my $out ;
my $method = "click" ;
my $dssp_file1;
my $dssp_file2;
my @dssp;
my $path_to_click;

my %SS = (
    "H" => "H",
    "G" => "H",
    "I" => "H",
    "B" => "E",
    "E" => "E",
    "T" => "c",
    "S" => "c",
    "X" => "c",
);

for (my $i = 0; $i < @ARGV ; $i++) {$options .= " $ARGV[$i] " }
if ($options =~ / -i\s+(\S+) /) {$inclick  = $1}
if ($options =~ / -o\s+(\S+) /) {$out      = $1}
if ($options =~ / -m\s+(\S+) /) {$method   = $1}
if ($options =~ / -dssp1\s+(\S+)/ ) {$dssp_file1  = $1}
if ($options =~ / -dssp2\s+(\S+)/) {$dssp_file2  = $1}
if ($options =~ / -path\s+(\S+)/) {$path_to_click  = $1}


my $upper_fasta_name = $ARGV[0]; 
my $lower_fasta_name = $ARGV[1]; 


my %CODE=(
        "A"=>"ALA", "R"=>"ARG", "N"=>"ASN", "D"=>"ASP",
        "C"=>"CYS", "E"=>"GLU", "Q"=>"GLN", "G"=>"GLY",
        "H"=>"HIS", "I"=>"ILE", "L"=>"LEU", "K"=>"LYS",
        "M"=>"MET", "F"=>"PHE", "P"=>"PRO", "S"=>"SER",
        "T"=>"THR", "W"=>"TRP", "Y"=>"TYR", "V"=>"VAL",
		"ALA"=>"A", "ARG"=>"R", "ASN"=>"N", "ASP"=>"D",
        "CYS"=>"C", "GLU"=>"E", "GLN"=>"Q", "GLY"=>"G",
        "HIS"=>"H", "ILE"=>"I", "LEU"=>"L", "LYS"=>"K",
        "MET"=>"M", "PHE"=>"F", "PRO"=>"P", "SER"=>"S", "MSE"=>"M",
        "THR"=>"T", "TRP"=>"W", "TYR"=>"Y", "VAL"=>"V"
);

my (%ALN,@ROUND);
my @aligned_sequences;
if ($inclick =~ /^$/){
	$pdb_1 = $ARGV[0] ; 
	$pdb_2 = $ARGV[1] ; 
}

if ($dssp_file2 =~ /\S+/){
	linearize ($dssp_file1,1);
	linearize ($dssp_file2,2);

}

my ($upper_index,$lower_index,$upper_fasta,$lower_fasta,@str_pairs) = superimpose($pdb_1,$pdb_2,$upper_fasta_name,$lower_fasta_name,$method);
my $length_ufasta = length($upper_fasta);
my $length_lfasta = length($lower_fasta);
prepare_to_leave();

sub linearize{
    my $dssp_file = $_[0];
	my $index = $_[1];
    my $switch = "off";
    open DSSP, "$dssp_file" or die "$! $dssp_file\n";
    while (my $line = <DSSP>){
        if ($line =~ /^\s+\#/){
            $switch = "on";
        }
        elsif ($switch eq "on"){
            $line =~/^(.....)(.....)(..)(..)..(.).*/;
			my $pos = $2 -1  ;
            my $chain = $3;
            my $aa = $4;
            my $ss = $5;
            $aa =~ s/\s+//g;
            $ss =~ s/\s+/X/;
			$dssp[$index][$pos] = $SS{$ss}; 
        }
    }
}

sub prepare_to_leave{
	system("rm -f background background1 background2 background3 dotplot plot lalign.temp.txt t.temp ");
	system("rm -f *.pdb.temp.fa ");
}

sub superimpose{
	my $pdb_1 = $_[0];
	my $pdb_2 = $_[1];
	my $upper_fasta_name = $_[2]; 
	my $lower_fasta_name = $_[3];
	my $method = $_[4];
	my $superimp  = "$pdb_1-$pdb_2"; 
	my $superimpr = "$pdb_2-$pdb_1";
	my $tm_score_s;
	my $tm_score_r;
	my @all_pairs;
	$superimp =~ s/\.pdb-/-/;
	$superimp =~ s/\.pdb$/.pdb.1.clique/;
	$superimpr =~ s/\.pdb-/-/;
	$superimpr =~ s/\.pdb$/.pdb.1.clique/;

	my ($fasta_1,$chain_1) = fasta_from_pdb($pdb_1); 
	my ($fasta_2,$chain_2) = fasta_from_pdb($pdb_2);

	if($method eq "click"){
		system("$path_to_click/click $pdb_1 $pdb_2 > t.temp ");
		$tm_score_s = quality_check($superimp,$chain_2,$chain_1);

		system("$path_to_click/click $pdb_2 $pdb_1 > t.temp ");
		$tm_score_r = quality_check($superimpr,$chain_2,$chain_1);
	}

	my ($pairsr,$pairs) ;
	if ($tm_score_s > $tm_score_r){
		if ($tm_score_r >= 0.3){pair($superimpr,"2","reverse")}
		if ($tm_score_s >= 0.3){$pairs  = pair($superimp,"2","sequential")}
		if ($tm_score_s < 0.3 and $tm_score_r < 0.3){print STDERR "Failed superimposition between $pdb_1 and $pdb_2\n"; prepare_to_leave(); exit}
		my @pairs  = split (/\n/,$pairs);
		push @all_pairs, @pairs;
	}elsif($tm_score_s <= $tm_score_r){
		if ($tm_score_s >= 0.3){pair($superimp,"2","sequential",$method)}
		if ($tm_score_r >= 0.3){$pairsr = pair($superimpr,"2","reverse",$method)}
		if ($tm_score_s < 0.3 and $tm_score_r < 0.3){print STDERR "Failed superimposition between $pdb_1 and $pdb_2\n"; prepare_to_leave(); exit}
		my @pairsr = split (/\n/,$pairsr);
		push @all_pairs, @pairsr;
	}

	@all_pairs = uniq @all_pairs;
	my @rlines  = gap_extension(@all_pairs);
	my @lines = dobles(@rlines);
	@lines = sort {$a <=> $b} @lines ;

	my ($Lower_line,$Upper_line);
	foreach  my $line (@lines){
		my @COL = split (/\t/,$line);
		$Upper_line .= "$COL[0]\n";
		$Lower_line .= "$COL[1]\n";
	}
	my ($tempname1,$tempname2,$upper_index,$lower_index) = subdivide($Lower_line,$Upper_line,$fasta_1,$fasta_2,$upper_fasta_name,$lower_fasta_name);
	return($upper_index,$lower_index,$fasta_1,$fasta_2,@all_pairs);
}

sub quality_check{
	my $superimp = $_[0];
	my $c1 = $_[1];
	my $c2 = $_[2];
	my $warning ;
    my ($cumulative,$ltarget,$d0,$rmsd,$matched,$in) ;
	$superimp =~ /(\S+)\-(\S+)\.pdb\.1\.clique/;
	my $query  = $1 ;
	my $target = $2 ;
	my $denominador = 0.1 ;
	my $first_raw  = "$target-$query.1.pdb\n";
    my $second_raw = "$query-$target.1.pdb\n";

    open CLI, "$superimp" or die "\n$!";
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
				if ($optionA > $optionB){$denominador = $optionB}
				else{$denominador = $optionA}
                my @OPTIONS;
                push @OPTIONS , $optionA ;
                push @OPTIONS , $optionB ;
                my @SOPTIONS = sort {$a <=> $b} @OPTIONS;
                $ltarget =  $OPTIONS[0];
            }
            $d0    = 1.24 * (($ltarget - 15) **(1/3)) - 1.8;
            if ($line =~ /^[A-Z0-9]\s+/){
                open FIRST, "$first_raw" or die "\n$!";
                while (my $f = <FIRST>){
                	my @F = split (" ",$f);
                	if ($f =~ /CA\s+$CODE{$LINE[2]}\s+$c2\s+$LINE[1]/) {$fx=$F[6]; $fy=$F[7]; $fz=$F[8]; last}
                }
                open SECOND, "$second_raw"or die "\n$!";
                while (my $s = <SECOND>){
                    my @S = split (" ",$s);
                    if ($s =~ /CA\s+$CODE{$LINE[6]}\s+$c1\s+$LINE[5]/) {$sx=$S[6]; $sy=$S[7]; $sz=$S[8]; last}
                }
                my $di = (sqrt(($fx-$sx)**2 + ($fy-$sy)**2 + ($fz-$sz)**2)) ;
                my $chunk = 1 / (1+( ($di /$d0) ** 2)) ;
                $cumulative = $cumulative + $chunk ;
            }

        }
        if ($ltarget <= 0){return("err")}
        else{
			my $unround =  $cumulative * (1/$denominador);
            my $tm = sprintf("%.4f", $unround);

            return($tm)
        }
}

sub dobles{
	my @RAW = @_;
	my @RAW = sort {$a <=> $b} @RAW;
	foreach my $dom_f (keys @RAW){
		my ($A, $B, $C) = split (" ",$RAW[$dom_f]);
		my ($AA, $BB) = split (" ",$RAW[$dom_f+1]);
		my $spot = $dom_f+1;
		if ($RAW[$dom_f+1] =~ /^$/){

		}
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

sub subdivide{
	my $lower_line = $_[0]; my @Lower_line = split (/\n/,$lower_line);
	my $upper_line = $_[1]; my @Upper_line = split (/\n/,$upper_line);
	my $upper_fasta = $_[2]; 
	my $lower_fasta = $_[3];
	my $upper_fasta_name = $_[4]; 
	my $lower_fasta_name = $_[5];
	my (@upper_fragment,@lower_fragment);
	my $tracker = 1;
	my $anterior ;
	my $upper_coords = "$upper_fasta_name\t";
	my $lower_coords = "$lower_fasta_name\t";
	for (my $i = 0 ; $i <= $#Lower_line  ; $i ++ ){
		my ($l_st,$l_end) = split (/:/,$Lower_line[$i]) ;
		my ($u_st,$u_end) = split (/:/,$Upper_line[$i]) ;
		my ($l_st_nxt,$l_end_nxt) = split (/:/,$Lower_line[$i+1]);
		my ($u_st_nxt,$u_end_nxt) = split (/:/,$Upper_line[$i+1]);
		my $gap_u = $l_st_nxt - $l_end -1;
		my $gap_l = $u_st_nxt - $u_end -1;
		my $d_total_covered = ($u_end_nxt - $u_st) * 1.2; ##### CHANGE TO 1.2??
		$upper_coords .= "$u_st-$u_end,";
		$lower_coords .= "$l_st-$l_end,";
		if ($gap_u >= 0){
			my @ltip = split (/:/,$lower_fragment[$tracker]);
			my $ltip = $ltip[$#ltip] ; $ltip  =~ s/\n//;
			my @utip = split (/:/,$upper_fragment[$tracker]);
			my $utip = $utip[$#utip] ; $utip  =~ s/\n//;
            ##print "uno\t";
			if ($tracker >= 1){
                ##print "uno.uno\t";
				my @other_ltip = split (/:/,$lower_fragment[$tracker-2]);
                my $other_ltip = $other_ltip[$#other_ltip] ; $other_ltip =~ s/\n//;
                my $distance_to_ltip = $l_st_nxt - $ltip ;
                my $distance_to_other_ltip = $l_st_nxt - $other_ltip ;

			}
			if ($ltip < $l_st_nxt and $utip < $u_st_nxt and $gap_u < $d_total_covered and $gap_l < $d_total_covered){
				$upper_fragment[$tracker] .= "$u_st:$u_end\n$u_st_nxt:$u_end_nxt\n";
				$lower_fragment[$tracker] .= "$l_st:$l_end\n$l_st_nxt:$l_end_nxt\n";
			}else{
                $upper_fragment[$tracker] .= "$u_st:$u_end\n";
				$lower_fragment[$tracker] .= "$l_st:$l_end\n";
				$tracker ++;
				$upper_fragment[$tracker] .= "$u_st_nxt:$u_end_nxt\n";
				$lower_fragment[$tracker] .= "$l_st_nxt:$l_end_nxt\n";
			}
			
		}elsif($gap_u < 0){
            my @ltip = split (/:/,$lower_fragment[$tracker]);
            my $ltip = $ltip[$#ltip] ; $ltip =~ s/\n//;
            my @utip = split (/:/,$upper_fragment[$tracker]);
            my $utip = $utip[$#utip] ; $utip =~ s/\n//;
			my ($l_st_skip,$l_end_skip) = split (/:/,$Lower_line[$i+2]);
			my ($u_st_skip,$u_end_skip) = split (/:/,$Upper_line[$i+2]);
			my $gap_u_skip = $l_st_nxt - $l_end_skip -1;
			$upper_fragment[$tracker] .= "$u_st:$u_end\n";
			$lower_fragment[$tracker] .= "$l_st:$l_end\n";
			$tracker ++;

			if($ltip < $l_st_nxt and $utip < $u_st_nxt){
				$upper_fragment[$tracker] .= "$u_st_nxt:$u_end_nxt\n";
				$lower_fragment[$tracker] .= "$l_st_nxt:$l_end_nxt\n";
				

				if ($gap_u_skip < $gap_u){
					$tracker -- ;
					$upper_fragment[$tracker] .= "$u_st_skip:$u_end_skip\n";
					$lower_fragment[$tracker] .= "$l_st_skip:$l_end_skip\n";
				}
			}
			
			else{
				$tracker++;
                $upper_fragment[$tracker] .= "$u_st_nxt:$u_end_nxt\n";
				$lower_fragment[$tracker] .= "$l_st_nxt:$l_end_nxt\n";
			}
		}
	}
	$upper_coords =~ s/,$//; 
	$upper_coords =~ /.*(.)\.pdb/; my $chain_1 = $1;$upper_coords =~ s/\.pdb/\t$chain_1/;
	$lower_coords =~ s/,$//;
	$lower_coords =~ /.*(.)\.pdb/; my $chain_2 = $1;$lower_coords =~ s/\.pdb/\t$chain_2/;
	open OCOORD, ">$upper_fasta_name-$lower_fasta_name.coord.txt";
	print OCOORD "$upper_coords\n$lower_coords\n";

	my $total_fragments = $#lower_fragment +1 ;
	my (@upper_bs,@lower_bs);
	my ($O1,$O2,$O3,$O4,$O5);
	for (my $j = 0 ; $j < $total_fragments ; $j ++){
		my @current_ufragment = split (/\n/,$upper_fragment[$j]);
		@current_ufragment = uniq @current_ufragment;
		my $current_u = join "\n", @current_ufragment;
		my @current_lfragment = split (/\n/,$lower_fragment[$j]);
		@current_lfragment = uniq @current_lfragment;
		my $current_l = join "\n", @current_lfragment;

		my ($upper_bi,$upper_bf,$lower_bi,$lower_bf,$o1,$o2,$o3,$o4,$o5) = align($current_l,$current_u,$upper_fasta,$lower_fasta,$upper_fasta_name,$lower_fasta_name);
		$O1 .= $o1;
		$O2 .= $o2;
		$O3 .= $o3;
		$O4 .= $o4;
		$O5 .= $o5;
		if ($upper_bi =~ /\S+/){push @upper_bs, $upper_bi}
		if ($upper_bf =~ /\S+/){push @upper_bs, $upper_bf}
		if ($lower_bi =~ /\S+/){push @lower_bs, $lower_bi}
		if ($lower_bf =~ /\S+/){push @lower_bs, $lower_bf}
	}
	print "  query SS  |$O1\n  query     |$O2\n            |$O3\n  target    |$O4\n  target SS |$O5\n\n";
	my @upper_fasta = split (//,$upper_fasta);
	my @lower_fasta = split (//,$lower_fasta);
	@upper_bs = sort {$a <=> $b} @upper_bs;
	@lower_bs = sort {$a <=> $b} @lower_bs;

	my $uf = join('',@upper_fasta[$upper_bs[0]-1..$upper_bs[$#upper_bs]-1]);
	my $lf = join('',@lower_fasta[$lower_bs[0]-1..$lower_bs[$#lower_bs]-1]);

	open O1, ">$upper_fasta_name.temp.fa";
	open O2, ">$lower_fasta_name.temp.fa";
	print O1 ">upper\n$upper_fasta\n";
	print O2 ">lower\n$lower_fasta\n";

	return("$upper_fasta_name.temp.fa","$lower_fasta_name.temp.fa",$upper_bs[0],$lower_bs[0]);

}

sub align{
	my $lower_line = $_[0]; my @Lower_line = split (/\n/,$lower_line);
    my $upper_line = $_[1]; my @Upper_line = split (/\n/,$upper_line);
	my $lower_fasta = $_[3]; my @Lower_fa  = split (//,$lower_fasta);
	my $upper_fasta = $_[2]; my @Upper_fa  = split (//,$upper_fasta);
	my $upper_fasta_name = $_[4]; 
	my $lower_fasta_name = $_[5];

    my (@upper_fragment,@lower_fragment);
    my (@LOWER,@UPPER,$u1,$l1,$ulast,$llast,@upper_bs,@lower_bs);
	my $limit ;
	my $warning ;
    my $tracker = 0;

    for (my $i = 0 ; $i < $#Upper_fa +1 ;$i ++){$UPPER[0][$i] = $Upper_fa[$i]}         
    for (my $i = 0 ; $i < $#Lower_fa +1 ;$i ++){$LOWER[0][$i] = $Lower_fa[$i]}

	if ($#Lower_line == 0){$limit =1;$warning="on"}
	else{$limit = $#Upper_line}
        for (my $i = 0 ; $i < $limit   ; $i ++ ){
            my ($l_st,$l_end) = split (/:/,$Lower_line[$i]) ;
            my ($u_st,$u_end) = split (/:/,$Upper_line[$i]) ;

            my ($l_st_nxt,$l_end_nxt) = split (/:/,$Lower_line[$i+1]);
            my ($u_st_nxt,$u_end_nxt) = split (/:/,$Upper_line[$i+1]);
			if ($tracker == 0){$u1 = $u_st ; $l1 = $l_st}
			$tracker ++;
            my $gap_u = $l_st_nxt - $l_end -1;
            my $gap_l = $u_st_nxt - $u_end -1;
            if ($gap_l >=1 and $gap_u >=1){
                for (my $g = 0 ; $g < $gap_l ; $g ++){
                    $LOWER[1][$l_end] .= "-";
                }
                for (my $g = 0 ; $g < $gap_u ; $g ++){
                    $UPPER[1][$u_end+$gap_l] .= "-";
                }
            }
            elsif ($gap_l >= 1){
                for (my $g = 0 ; $g < $gap_l ; $g ++){
                    $LOWER[1][$l_end] .= "-";
                }
            }
        	elsif ($gap_u >= 1){ 
                for (my $g = 0 ; $g < $gap_u ; $g ++){
                    $UPPER[1][$u_end] .= "-";
                }
            }
			$ulast = $u_end_nxt;
			$llast = $l_end_nxt;
		if ($warning eq "on"){
			$ulast = $u_end;$llast = $l_end;$u1=$u_st;$l1=$l_st
		}
	}
	push @upper_bs, $u1;
	push @upper_bs, $ulast;
	push @lower_bs, $l1;
	push @lower_bs, $llast;
	
	my $line_SS_query;
	my $line_upper_fasta;
	my $line_tics;
	my $line_lower_fasta;
	my $line_SS_target;

	if ($u1 =~ /\S+/){
		if (($ulast - $u1) < 1){}
		else{
			my $upper_dssp ;
			my $lower_dssp;
			$upper_fasta_name =~ s/\.pdb$//;
			$lower_fasta_name =~ s/\.pdb$//;
			#$line_SS_query  .= sprintf '%10s',"SS_query";
			$line_SS_query .= sprintf '%6s';
			for (my $i = $u1-1 ; $i < $ulast ; $i ++){
				$line_SS_query .= "$UPPER[1][$i]$dssp[1][$i]";
				$upper_dssp .= "$UPPER[1][$i]$dssp[1][$i]";
			}

			$line_SS_query .= sprintf '%+9s', "|";
			#$line_upper_fasta .= sprintf '%10s', "$upper_fasta_name"  ;
    		$line_upper_fasta .= sprintf '%6s', " $u1 ";
	        for (my $i = $u1-1 ; $i < $ulast ; $i ++){
				$line_upper_fasta .= "$UPPER[1][$i]$UPPER[0][$i]";
			}
			$line_upper_fasta .= sprintf '%+6s', " $ulast ";

			$line_upper_fasta .= sprintf '%+3s' , "|";
			#$line_tics .= sprintf '%10s',"";
			$line_tics .= sprintf '%6s';
			for (my $i = $l1-1 ; $i < $llast ; $i ++){
				$lower_dssp .= "$LOWER[1][ $i]$dssp[2][$i]";
			}
			my @upper_dssp = split (//,$upper_dssp);
			my @lower_dssp = split (//,$lower_dssp);
    	   	for (my $i = 0 ; $i < $#upper_dssp + 1; $i ++){
				if($lower_dssp[$i] =~ /-/ or $upper_dssp[$i] =~ /-/){$line_tics .= " "}
				else{
					if ($lower_dssp[$i] eq $upper_dssp[$i]){
						if($lower_dssp[$i] =~ /[A-Z]/){$line_tics .= ":"}
						else{$line_tics .= ":"}
					}
					else{$line_tics .= "."}
				}
			}

			$line_tics .= sprintf '%+9s', "|";
			#$line_lower_fasta .= sprintf '%10s', $lower_fasta_name  ;
	       	$line_lower_fasta .= sprintf '%6s'," $l1 ";
    	   	for (my $i = $l1-1 ; $i < $llast ; $i ++){
				$line_lower_fasta .= "$LOWER[1][ $i]$LOWER[0][$i]";
			}
			$line_lower_fasta .= sprintf '%+6s', " $llast ";	
			$line_lower_fasta .= sprintf '%+3s', "|";
			#$line_SS_target .= sprintf '%10s', "SS_target";
			$line_SS_target .= sprintf '%6s';
			$line_SS_target .= "$lower_dssp";			
			$line_SS_target .= sprintf '%+9s', "|";
		}
	}

	#print "$line_SS_query\n";
	#print "$line_upper_fasta\n";
	#print "$line_tics\n";
	#print "$line_lower_fasta\n";
	#print "$line_SS_target\n\n";


	@upper_bs = sort {$a <=> $b} @upper_bs;
	@lower_bs = sort {$a <=> $b} @lower_bs;
	my $upper_bi = $upper_bs[0];
	my $upper_bf = $upper_bs[$#upper_bs];
	my $lower_bi = $lower_bs[0];
	my $lower_bf = $lower_bs[$#lower_bs];
	return($upper_bi,$upper_bf,$lower_bi,$lower_bf,$line_SS_query,$line_upper_fasta,$line_tics,$line_lower_fasta,$line_SS_target);
}


sub fasta_from_pdb{
	my $inpdb    = $_[0];
	my @sequence ;
	my $sequence ;
	my $chain    ;

	open IN, $inpdb or die "$! $inpdb\n";
	while (my $line = <IN>){
		if ($line =~ /REMARK RES ([A-Z]{3}) (\S) (.*)$/){
			my $aa     = $CODE{$1};
			my $c      = $2;
			my $resnum = $3;
			$sequence[$resnum]= $aa;
			$chain = $c;
		}
		
	}
	
	for(my $i = 1 ; $i <= $#sequence; $i ++){
		if ($sequence[$i] =~ /^$/){$sequence[$i] = "X"}
		$sequence .= $sequence[$i];
	}
	
	return($sequence,$chain);
}

sub pair{
    my $clique              = @_[0];
    my $eval                = @_[1];
	my $order               = @_[2];
	my $method              = @_[3];
    my ($resnum1,$resnum2,$last2);
    my $last1               = "INI";
    my ($dir,$rep_p,$rep_n) = 0;
    my $pos                 = -1;
    my $round;
    my $pairs ;
	my $type ;
    open CLICK, "$clique" or print "$! $clique\n" ;
    while (my $line = <CLICK> ){
        chomp $line;
		if ($line =~ /^MATCHED ATOM PAIRS BETWEEN/){$type = "click"}
        if ($line =~ /^\S\S? / or $line =~ /^\s+[0-9]+\s+/){
			my ($chain1,$resnum1,$chain2,$resnum2);

            if ($order eq "sequential" and $type eq "click"){
				($chain1,$resnum1,undef,undef,$chain2,$resnum2,undef) = split (" ",$line,7);
			}elsif($order eq "sequential"){
				($resnum1,undef,$chain1,$resnum2,undef,$chain2,undef) = split (" ",$line,7);
			}
			if ($order eq "reverse" and $type eq "click"){
				($chain2,$resnum2,undef,undef,$chain1,$resnum1,undef) = split (" ",$line,7);
			}elsif($order eq "reverse"){
				($resnum1,undef,$chain1,$resnum2,undef,$chain2,undef) = split (" ",$line,7);
			}
            if ($last1 eq "INI"){
                $dir   = 0;
				$last1 = $resnum1;
				$last2 = $resnum2;
                $rep_p = 1;
                $rep_n = 1;
                next;
            }

            if ($resnum1 == $last1 + 1){
                $dir = $resnum2 - $last2;
                if ($dir == 1){
                    $last1 = $resnum1;
                    $last2 = $resnum2;
                    $rep_p ++;
                }elsif ($dir == -1){
                    $last1 = $resnum1;
                    $last2 = $resnum2;
                    $rep_n ++; ###
                }else{
                    $last1 = $resnum1;
                    $last2 = $resnum2;
                    $rep_p = 1;
                    $rep_n = 1;
                }
            if ($rep_p == $eval or $rep_n == $eval){
                my ($res1,$res2);
                for (my $i = $eval - 1; $i >= 0; $i-- ) {
                    if ($rep_p == $eval){
                        $res1 = $resnum1 - $i ;
                        $res2 = $resnum2 - $i ;
			            $pairs .= "$res1\t$res2\n";
                    }
                    if ($rep_n == $eval){
                        $res1 = $resnum1 - $i ;
                        $res2 = $resnum2 + $i ;
                    }
                }
            }elsif ($rep_p > $eval) {$pairs .= "$resnum1\t$resnum2\n"}

            }else{
                $last1 = $resnum1;
                $last2 = $resnum2;
                $rep_p = 1;
                $rep_n = 1;
            }
			print BACK "$resnum1\t$resnum2\n";
        }
    }
    return($pairs);
}

sub gap_extension{
	my @inpairs     = @_ ;
	my $reward      = 2 ;  
	my $gap_penalty = 1 ;
	my $threshold   = 1 ; 
	my $options     = "";
	my @matrix      ;

	my (%DIAGONAL,@MAT);

	foreach my $lin (@inpairs){
        my ($i, $j)   = split (" ", $lin);
		$MAT[$i][$j]  = 1;
        my $ind       = maximal($i, $j, $gap_penalty, $reward, $threshold);

		store($i,$j);

		if($ind ==1){
            my @coord ;
            $coord[0] = $i;
            $coord[1] = $j;
            my $name  = "$i:$j";
            push (@{$DIAGONAL{$name}}, @coord);
            push (@{$DIAGONAL{$name}}, $MAT[$i][$j]);
		}
	}

	sub store{
    	my $x = shift (@_);
    	my $y = shift (@_);
		foreach my $k (sort {$a<=>$b} keys %DIAGONAL){
    	    my $anteriorx = ${$DIAGONAL{$k}}[0];
    	    my $anteriory = ${$DIAGONAL{$k}}[1];
    	    if ($x-1 == $anteriorx and $y-1 == $anteriory){
    	        my @coord ;
    	        $coord[0] = $x ;
    	        $coord[1] = $y ;
    	        unshift( @{$DIAGONAL{$k}}, @coord );
    	        pop @{$DIAGONAL{$k}};
    	        push (@{$DIAGONAL{$k}}, $MAT[$x][$y] );
    	    }
    	}
	}

	sub maximal{
    	my $x            = shift (@_);
    	my $y            = shift (@_);
		my $gap_penalty  = shift (@_);
    	my $reward       = shift (@_);
		my $threshold    = shift (@_);
    	my ($ai, $bi, $c)  = 0;

    	if ($MAT[$x-1][$y-1] > 0)  {$ai = $MAT[$x-1][$y-1] + $reward}
		$MAT[$x][$y] = max($MAT[$x][$y], $ai, $bi, $c);
        
		return($MAT[$x][$y]);
	}
	

	foreach my $k (sort {$a<=>$b} keys %DIAGONAL){
    	my $baseline = pop @{$DIAGONAL{$k}};
    	if ($baseline >= $threshold){
    	    my $ai = 0;
    	    my $bi = 1;
    	    if (${$DIAGONAL{$k}}[$ai] ne "" and ${$DIAGONAL{$k}}[$bi] ne ""){
    	        my ($vector1_ini,$vector2_ini)= split(":",$k);
        	    push @matrix,  "$vector1_ini:${$DIAGONAL{$k}}[$ai]\t$vector2_ini:${$DIAGONAL{$k}}[$bi]\t$baseline\n";
        	    $ai ++; $ai ++; $bi ++; $bi ++;
        	}
    	}
	}

	return(@matrix);

}


__END__
