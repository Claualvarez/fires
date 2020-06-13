#########################################################################
#									#
#		FiRES: Find Repeated Elements in Structure		#
#   A computational method for the de novo identification of repeats	#
#	             Author: Claudia Alvarez-Carreno			#
#			calvarez(at)ifc.unam.mx				#
#########################################################################


my $path_to_click="/usr/bin/Click_X86_64";
my $path_to_dssp="/usr/bin";


#########################################################################

use strict ;
#use integer;
use List::Util qw[min max] ;
use List::MoreUtils qw(uniq);

my $path_to_fires ;
my $outname ;
my $dssp;                # -dssp  name of the dssp file
my $pdb;                 # -pdb   name of the pdb file
my $chaini      = "A"  ; # -c     chain of interest, default all
my $pdbcode            ; # -f     pdb 3-letter code
my $min         = 1    ; # -min   minimum number of secondary structural elements in the initial fragmentation
my $max         = 5    ; # -max   maximum number of secondary structural elements in the initial fragmentation (antes 7)
my $divider     = 4    ; #        number of times the total number of secondary structural elements are going to be divided to determine the initial fragmentation
my $rounds      = 15   ; # -i     maximum number of iterations
my $eval        = 3    ; #        minimum number of sequential elements to be incorporated as a diagonal
my $reward      = 2    ; # -ext   rewarding points for a diagonal extension
my $overlap     = 75   ; # -o     minimal Structure Overlap Score on the last click evaluation
my $coverage    = 75   ; # -n     minimal percentage of overlapped atoms of the shortest chain over the longest chain compared 
my $separation  = 8    ; # -d     minimal distance between adjacent  diagonals to be considered for diagonal-joining
my $gap         = 55   ; # -gap   gap spaces allowed between adjacent diagonals
my $ss          = 4    ; # -ss    minimal number of secondary structure elements in last evaluation
my $threshold   = 10   ; # -b     minimum score for a diagonal (antes 18)
my $gap_penalty = 1    ; # -p
my $sampling    = 15   ; # -S     number of new coordinates in next round
my $redundancy  = 150  ; # -R     redundancy of a single coordinate
my $rmsd        = 3    ; # -rmsd  antes 2.5
my $min_aa      = 20   ; # -ma
my $mma         = 20   ; # -mma 
my $force       = "yes"; # -force antes "off"
my @CHAINS;
my $chain;
my $warning ;
my $input ;
my $options     = ""   ;
my $chainname ;
my $directory ;

for (my $i = 0; $i < @ARGV ; $i++)       {$options    .= " $ARGV[$i] " }

if ($options =~ / -pdb\s+(\S+) /)        {$pdb           = $1}
if ($options =~ / -dssp\s+(\S+) /)       {$dssp          = $1}
if ($options =~ / -c\s+(\S+) /)          {$chaini        = $1}
if ($options =~ / -b\s+(\S+) /)          {$threshold     = $1}
if ($options =~ / -c\s+(\S+) /)          {$chaini        = $1}
if ($options =~ / -cu\s+(\S+) /)         {$chaini        = uc $1}
if ($options =~ / -d\s+(\S+) /)          {$separation    = $1}
if ($options =~ / -dssp\s+(\S+) /)       {$dssp	         = $1}
if ($options =~ / -ext\s+(\S+) /)        {$reward        = $1}
if ($options =~ / -i\s+(\S+) /)          {$rounds        = $1}
if ($options =~ / -pdb\s+(\S+) /)        {$pdb           = $1}
if ($options =~ / -max\s+(\S+) /)        {$max           = $1}
if ($options =~ / -min\s+(\S+) /)        {$min           = $1}
if ($options =~ / -ma\s+(\S+) /)         {$min_aa        = $1}
if ($options =~ / -mma\s+(\S+) /)        {$mma           = $1}
if ($options =~ / -o\s+(\S+) /)          {$overlap       = $1}
if ($options =~ / -n\s+(\S+) /)          {$coverage      = $1}
if ($options =~ / -p\s+(\S+) /)          {$gap_penalty   = $1}
if ($options =~ / -R\s+(\S+) /)          {$redundancy    = $1}
if ($options =~ / -S\s+(\S+) /)          {$sampling      = $1}
if ($options =~ / -gap\s+(\S+) /)        {$gap           = $1}
if ($options =~ / -ss\s+(\S+) /)         {$ss            = $1}
if ($options =~ / -f\s+(\S+)\s?/)        {$pdbcode       = $1}
if ($options =~ / -rmsd\s+(\S+)\s?/)     {$rmsd          = $1}
if ($options =~ / -input\s+(\S+)\s?/)    {$input         = $1}
if ($options =~ / -force\s+(\S+)\s?/)    {$force         = $1}
if ($options =~ / -path\s+(\S+)\s?/)     {$path_to_fires = $1}
if ($options =~ / -pdssp\s+(\S+)\s?/)    {$path_to_dssp  = $1}
if ($options =~ / -pclick\s+(\S+)\s?/)   {$path_to_click = $1}


######################################################################

print STDOUT "\tInput file : $pdb\n";

######################################################################
system("cp $path_to_fires/Parameters.inp . ");

if ($pdb =~/\w.*/ and $chaini =~ /\S+/ ){
	$outname   = $pdb ; 
	$pdbcode   = int(rand(10000)) ;

	if ($chaini =~ /\-path/){print STDERR "Invalid input\nUsage: perl fires.pl -pdb <pdbid> -c <chain> [options]\n";exit}

	$chainname = $chaini;
	my $code   = $pdbcode;
	$dssp      = "$code\_$chainname.dssp";
	prepare_pdb("$outname",$chainname,"$code\_$chainname.pdb",$path_to_dssp,$dssp);
	$pdb       = $pdbcode	;
	system ("cp $code\_$chainname.pdb $pdbcode.pdb ; cp $code\_$chainname.dssp $pdbcode.dssp");
	$dssp      = "$pdbcode.dssp" ;
	$pdb       = "$pdbcode.pdb" ;
	$chaini    = "A";
}
else{
	print "Error: invalid input\nUsage: perl fires.pl -pdb <pdbid> -c <chain> [options]\n";exit
}


push @CHAINS, $chaini ;

my $onemore = $rounds +1 ;
my $oneless = $rounds -1 ;

foreach $chain (@CHAINS){
#	my    ($reference,$warning) = SSDivide::SS($dssp,$chain,$pdbcode,$min,$max,$divider);
	my  ($reference,$warning) = SS($dssp,$chain,$pdbcode,$min,$max,$divider);
	for (my $i = 0 ; $i < $onemore ; $i++){

		open  OUT, ">>$pdbcode\_$chain.reference";
        	print OUT  "Ref\n";
        	close OUT;

		serialClick($reference,$i,$pdbcode,$chain,$path_to_click);
		$reference = "";
		limpia($i,$pdbcode,$eval,$chain);

		system("perl $path_to_fires/gap-extension.pl -i $pdbcode\_$chain.kmer.$i.sorted -ext $reward -p $gap_penalty -b $threshold > $pdbcode\_$chain.diag.$i");
		system("perl $path_to_fires/new_query-target.pl $pdbcode\_$chain.diag.$i> $pdbcode\_$chain.domino.$i");
		system("perl $path_to_fires/query-target.pl -i $pdbcode\_$chain.domino.$i -dssp $dssp -pdb $pdbcode > $pdbcode\_$chain.reference.$i");

		open IN, "$pdbcode\_$chain.reference.$i" or die "$pdbcode\_$chain.reference.$i \n$!";
		while (my $line = <IN>) {
			my $token = 0;
			open COMP, "$pdbcode\_$chain.reference" or die "$pdbcode\_$chain.reference \n$!";
			while (my $r = <COMP>) {
				chomp $line;
				chomp $r;
				if ($r eq $line) {$token ++;last}
			}
			if ($token == 0){$reference .= "$line\n"}
		
		}
		close IN;
	

		if ($i >= 1){		
			my $j = $i - 1;
			system("cat $pdbcode\_$chain.kmer.* >> $pdbcode\_$chain.kmer ");
			system("perl $path_to_fires/new_info.pl -R $redundancy -S $sampling -i $pdbcode\_$chain.kmer.$j.sorted -j $pdbcode\_$chain.kmer.$i.sorted -k $pdbcode\_$chain.kmer > score ");

			open IN, "score" or die "ERROR APERTURA score \n$!";
			while (my $line = <IN>) {
				chomp $line;
				if ($line eq "NO" or $i == $oneless) {
					system ("cat $pdbcode\_$chain.kmer*sorted > kmersorted ");
					my $elema = "kmersorted";
					my $elemb = "$pdbcode\_$chain.dotplot";
					csuc($elema,$elemb);

					system ("perl $path_to_fires/square_mtx.pl     $pdbcode\_$chain.dotplot     > $pdbcode\_$chain.dotplot.int");
					system ("perl $path_to_fires/smith-waterman.pl $pdbcode\_$chain.dotplot.int > $pdbcode\_$chain.dotplot.out");
					my $elema = "dotplotin";
					my $elemb = "$pdbcode\_$chain.dotplot.in";					
					system ("cat $pdbcode\_$chain.dotplot.int > dotplotin");
					csuc($elema,$elemb);

					system ("perl $path_to_fires/gap-extension.pl -i $pdbcode\_$chain.dotplot.in > $pdbcode\_$chain.verbinden");
					system ("cut -f 1,2 *$pdbcode\_$chain*domino* | awk '(NF==2)' > $pdbcode\_$chain.2verbinden");
					system ("cut -f 1,3 *$pdbcode\_$chain*domino* | awk '(NF==2)' >> $pdbcode\_$chain.2verbinden");
					system ("cut -f 2,3 *$pdbcode\_$chain*domino* | awk '(NF==2)' >> $pdbcode\_$chain.2verbinden");
					system ("perl $path_to_fires/transitivity-cluster.pl -i $pdbcode\_$chain.dotplot.out -c $chainname -dssp $dssp -add $pdbcode\_$chain.2verbinden -gap $gap -d $separation -ss $ss -force $force > $pdbcode\_$chain.domino.2last");
					system ("perl $path_to_fires/transitivity-cluster.pl -i $pdbcode\_$chain.dotplot.out -c $chainname -dssp $dssp -add $pdbcode\_$chain.verbinden -gap $gap -d $separation -ss $ss -force $force > $pdbcode\_$chain.domino.last");
					system ("perl $path_to_fires/candidates.pl -i $pdbcode\_$chain.domino.last -dssp $dssp -pdb $pdbcode -c $chain > $pdbcode\_$chain.repeats.simple.list");
					system ("perl $path_to_fires/candidates.pl -i $pdbcode\_$chain.domino.2last -dssp $dssp -pdb $pdbcode -c $chain > $pdbcode\_$chain.repeats.simple.2list");
					my $elema = "simple_list";
					my $elemb = "$pdbcode\_$chain.repeats.simple_list";
					system ("cat *repeats.simple* > simple_list ");
					csuc($elema,$elemb);
					system ("perl $path_to_fires/se_click.pl -i $pdbcode\_$chain.repeats.simple_list -c $chain -click $path_to_click ");
					system ("perl $path_to_fires/last-eval.pl -i $pdbcode\_$chain.repeats.simple_list -c $chain -ch $chaini -o $overlap -n $coverage -rmsd $rmsd -ma $min_aa -mma $mma  > $pdbcode.$chainname.temp ");
					system ("cat candidates.txt > tab.temp");					
					system ("perl $path_to_fires/read_output.pl -i $pdbcode.$chainname.temp -c $chaini -ch $chainname > $pdbcode.$chainname.out ");
					$i = $rounds;
					last;
				}
				
			}
				
		}
	}
	

}

my $suboptimal = "$outname.tab";

suboptimal("candidates.txt",$suboptimal);
system ("perl $path_to_fires/clean_table.pl -i $outname.tab -c $chaini > $outname\_$chainname.best ");
system(" cp $pdbcode\_$chainname.pdb $outname.$chainname.pdb");
#system("cat $outname\_$chainname.best");
#system ("perl $path_to_fires/cater.pl -path $path_to_fires -out $outname -pdb $pdbcode -i $outname\_$chainname.best -ch $chainname -c $chaini ");
system ("perl $path_to_fires/cater.pl -click $path_to_click -path $path_to_fires -out $outname -pdb $pdbcode -i $outname\_$chainname.best -ch $chainname -c $chaini ");

system("mv $outname.tab.out $outname.out");
system("rm $outname.$chainname.pdb $outname\_$chainname.best");

my %CODE=(
        "A"=>"ALA", "R"=>"ARG", "N"=>"ASN", "D"=>"ASP",
        "C"=>"CYS", "E"=>"GLU", "Q"=>"GLN", "G"=>"GLY",
        "H"=>"HIS", "I"=>"ILE", "L"=>"LEU", "K"=>"LYS",
        "M"=>"MET", "F"=>"PHE", "P"=>"PRO", "S"=>"SER",
        "T"=>"THR", "W"=>"TRP", "Y"=>"TYR", "V"=>"VAL"
);

sub prepare_pdb{
	my $infile       = @_[0];
	my $chain        = @_[1];
	my $outfile      = @_[2];
	my $path_to_dssp = @_[3];
	my $dssp         = @_[4];
	my $number ;
	my $flag = 0;
	
	open IN, "$infile" or die "\n$!";
	open OUT, ">tempoutfile";
	open TMPR, ">tempremarks";
	my $remark_res_token = 0;
	my $new_aa_pos = 0;

	while (my $line = <IN>){
		if($line =~ /^REMARK RES/){$remark_res_token ++;print OUT $line}
		if($line =~ /^ATOM\s+/){
			$line =~ /(^ATOM.)(......)(.....)(....)(..)(....)(.*)/;
			my $head = $1 ;
			my $nb   = $2 ;
			my $atm  = $3 ;
			my $aa   = $4 ; 
			my $c    = $5 ;
			my $aa_pos = $6;
			my $inf  = $7 ;
			$c =~ s/ //g;
			if ($c eq $chain and $atm =~ /\s+N\s+/){
				$new_aa_pos ++;
				if ($remark_res_token < 1) {open TMPR, ">>tempremarks";} 
				print  TMPR "REMARK RES$aa A";
				printf TMPR "%5s", $new_aa_pos;
				print  TMPR "\n";
			}
			if ($c eq $chain and $flag == 0){
				$flag ++;	
				$number = $nb ;
			}
			if ($c eq $chain){
				if ($remark_res_token < 1) {$aa_pos = $new_aa_pos}
				if ($aa =~ /^A/){
					$aa =~ s/^A//;
					print  OUT  "$head";
					printf OUT '%6s',$number;
					print  OUT "$atm $aa";
					printf OUT '%2s', "A";
					printf  OUT "%4s", $aa_pos;
					print  OUT "$inf\n";
					$number ++;
				}elsif($aa =~ /^\s/){
					print  OUT "$head";
					printf OUT '%6s',$number;
					print  OUT "$atm$aa";
					printf OUT '%2s', "A";
					printf  OUT "%4s", $aa_pos;
					print  OUT "$inf\n";
					$number ++
				}
			}

		}elsif($line =~ /^HETATM/){
			$line =~ /(^HETATM)(.....)(.....)(....)(..)(....)(.*)/;
            my $head = $1 ;
    		my $nb   = $2 ;
    		my $atm  = $3 ;
            my $aa   = $4 ;
            my $c    = $5 ;
			my $aa_pos = $6;
            my $inf  = $7 ;
            $c =~ s/ //g;
			if ($c eq $chain and $atm =~ /\s+N\s+/){
				$new_aa_pos ++;
				close TMPR;
				if ($remark_res_token < 1) {open TMPR, ">>tempremarks";} 
				print  TMPR "REMARK RES$aa A";
				printf TMPR "%5s", $new_aa_pos;
				print  TMPR "\n";	
			}

            if ($c eq $chain and $flag == 0){
                $flag ++;
                $number = $nb ;
            }
            if ($c eq $chain){
				if ($remark_res_token < 1) {$aa_pos = $new_aa_pos}
                if ($aa =~ /^A/){	
                    $aa =~ s/^A//;
                    print  OUT  "ATOM  ";
                    printf OUT '%5s',$number;
                    print  OUT "$atm $aa";
                    printf OUT '%2s', "A";
					printf  OUT '%4s', $aa_pos;
                    print  OUT "$inf\n";
                    $number ++;
            	}elsif($aa =~ /^\s/){
                    print  OUT "ATOM  ";
                    printf OUT '%5s',$number;
                    print  OUT "$atm$aa";
                    printf OUT '%2s', "A";
					printf  OUT "%4s", $aa_pos;
                    print  OUT "$inf\n";
                    $number ++;
				}
            }
		}elsif($flag > 0 and $line =~ /^TER/){last}
	}
	if ($remark_res_token < 1){system ("cat tempremarks  >  $outfile")}
	system ("cat  tempoutfile >  $outfile");
	close OUT ;
	system("$path_to_dssp/dssp $outfile > $dssp");
}

######
sub SS{
	my $dssp         = @_[0];
	my $stored_str   = "FIRST";
	my $stored_chain = @_[1];
	my $pdbcode      = @_[2];
	my $min          = @_[3];
	my $max          = @_[4];
	my $divider      = @_[5];  
	my $count        = 0;
	my $warning;
	my $res_nbf;
	my $initial;
	my @SS;
	my $result;
	
	open IN, "$dssp" or die "dssp file:$dssp \n$!" ;
	while (my $line = <IN>){
		if ($line =~/\#/) {$initial ="INICIO"}
		if (($initial eq "INICIO") and ($line =~/^\s+[0-9]+/)){
			chomp $line;
			my ($nb, $res_nb, $chain, $aa, $str) = standard($line);
			if ($line =~/!/) {next}
			if ($str eq $stored_str and $stored_chain eq $chain) { $SS[$count] = "$SS[$count]$chain\t$CODE{$aa}\t$res_nb\n" }
			if ($aa eq "X") { next }
			if ($str ne $stored_str and $stored_chain eq $chain){
				$stored_str = $str;
				if ($str eq "H" or $str eq "B" or $str eq "S") {$count++}
				$SS[$count] = "$SS[$count]$chain\t$CODE{$aa}\t$res_nb\n";
			}
			$res_nbf = $res_nb;
		}
	}
	close IN;

	if    ($res_nbf < 130) {$warning = "small"}
	elsif ($res_nbf < 450) {$warning = "large"}
	else                   {$warning = "off"}

	my $totalelem = $#SS;
	my $i_size    = int($totalelem / $divider);
	my $size;

	if($i_size <= 2){ $size = 2}	
	else{
		if    ($i_size < $min)                    {$size = $min}
		elsif ($i_size > $min and $i_size < $max) {$size = $i_size}
		else                                      {$size = $max}
	}	

	foreach my $key (keys @SS){	
		my @FIRST          = split ("\n",$SS[$key]);
		my ($chain1,undef) = split (" ",$SS[$key],2);
		my @LAST           = split ("\n",$SS[$key+$size]);
		my ($chain2,undef) = split (" ",$SS[$key+$size],2);
			
		if ($chain1 eq $chain2) {$result .= "$pdbcode\t$FIRST[0]\t$LAST[$#LAST]\n"}
	}
	return($result, $warning);
}

############
sub standard{
	my $line = @_[0];
	$line =~/.(....).(....).(.).(.)..(.)/;
        my $nb     = $1;
        my $res_nb = $2;
        my $chain  = $3;
        my $aa     = $4;
        my $str    = $5;
        $aa        =~ s/[a-z]/C/g;
        $str       =~ s/ /-/;
        $str       =~ s/[GI]/H/;
        $str       =~ s/E/B/;
        $str       =~ s/S/T/;
        $res_nb    =~ s/ //g;

	return($nb, $res_nb, $chain, $aa, $str);

}

###############
sub serialClick{
	my $reference = @_[0];
	my $i         = @_[1];
	my $pdbcode   = @_[2];
	my $chain     = @_[3];
	my $path_to_click = @_[4];
	my @REFS      = split ("\n", $reference);
	
	open OUT, ">$pdbcode\_$chain.list.$i";
	foreach my $args (@REFS) {
		open  REF, ">>$pdbcode\_$chain.reference";
		print REF  "$args\n";
		my $click = corta($args,$path_to_click); 
		print OUT  "$click\n";
	}
	close OUT;	
}


#########
sub corta{
	my $params        = $_[0];
	my $path_to_click = $_[1];
	my ($pdb,$chain_o,$ini_aa,$ini_pos,$chain_check,$last_aa,$last_pos) = split("\t",$params);
	my $last       = $last_pos+1;
	my $stored_aa  = "X";
	my $stored_pos = "X";
	my $pace       = 0 ;
	my $cuenta     = 0 ;
	my $stored_chain;	
	my $query      = "q_$pdb\_$chain_o\_$ini_aa$ini_pos\_$last_aa$last_pos.pdb";
	my $target     = "t_$pdb\_$chain_o\_$ini_aa$ini_pos\_$last_aa$last_pos.pdb";

	open OUTPDB,">$query";
	open TRUNCATED, ">$target";
	open INPDB,"$pdb.pdb"  ;

	while (my $line = <INPDB> ){
		if ($line =~/^ATOM\s+[0-9]+\s+\w+\s+(\w+)\s+(\w+)\s+([0-9]+)\s+/){
			my $aa      = $1; 
			my $chain   = $2; 
			my $pos     = $3; 
			$line       =~ /^ATOM\s+[0-9]+(\s+\w+\s+\w+\s+\w+\s+[0-9]+\s+.*)/;
			my $columns = $1;
			
			if ($chain_o eq $chain and $ini_aa eq $aa and $ini_pos == $pos) {print OUTPDB "$line"}
			elsif ($chain_o eq $chain and $pos > $ini_pos and $pos <= $last_pos) {print OUTPDB "$line"}
			elsif ($chain_o eq $chain) {print TRUNCATED "$line"}
		}
	}

	close INPDB;
	print TRUNCATED "END\n";
	print OUTPDB "END\n";
	close OUTPDB;
	close TRUNCATED;

	click($query,$target,$path_to_click);
}

#########
sub click{
	my ($query,$target,$path_to_click) = @_; 
	my $query_name      = $query;
	$query_name         =~ s/\.pdb//; 
	system ("$path_to_click/click $query $target -a 1 >> click-file ");
	my $outfile         = "$query_name-$target\.1\.clique";
	return ($outfile);
}


##########
sub limpia{
	my $i                   = @_[0];
	my $pdbcode             = @_[1]; 
	my ($resnum1,$resnum2,$last2);
	my $eval                = @_[2];
	my $chain               = @_[3];
	my $last1               = "INI";
	my ($dir,$rep_p,$rep_n) = 0;
	my $pos                 = -1;
	my $round;
	my $output_lines;
	my $tm;
	my $tm_threshold = 0.30;

	open OUT, ">$pdbcode\_$chain.kmer.$i.s";
	open LIST, "$pdbcode\_$chain.list.$i" or die "$pdbcode\_$chain.list.$i \n$!";

	while (my $clique = <LIST> ){
        open CLICK, "$clique" or print "$! $clique\n" ;
		my $query  ;
		my $target ;
		my $first_raw;
		my $second_raw;
		my $rmsd;
		my $ltarget;
		my $c = $chain;
		my ($cumulative,$ltarget,$d0,$rmsd,$matched,$in) ;

		while (my $line = <CLICK> ){
			chomp $line; 
			my @LINE = split (" ",$line);
        	my ($fx,$fy,$fz,$sx,$sy,$sz);

			if ($line =~ /^Number of.*of\s+(\S+)\s+and\s+(\S+)\s+=\s+([0-9]+)\s+and\s+([0-9]+)/){
				$query      = $1;
				$target     = $2;
				$first_raw  = "$target-$query.1.pdb\n";
    			$second_raw = "$query-$target.1.pdb\n";
				my $optionA = 0 + $3;
           		my $optionB = 0 + $4;
            	my @OPTIONS;
            	push @OPTIONS , $optionA ;
            	push @OPTIONS , $optionB ;
            	my @SOPTIONS = sort {$a <=> $b} @OPTIONS;
            	$ltarget =  $OPTIONS[0];
			}
        	if ($line =~ /RMSD =\s+(\S+)/){$rmsd = $1}
        	$d0    = 1.24 * (($ltarget - 15) **(1/3)) - 1.8; 
        	if ($line =~ /^$c/){
            	open FIRST, "$first_raw" or die "\n$!";
            	while (my $f = <FIRST>){
                	my @F = split (" ",$f);
                	if ($f =~ /CA\s+$CODE{$LINE[2]}\s+$c\s+$LINE[1]/) {$fx=$F[6]; $fy=$F[7]; $fz=$F[8]; 
					last}
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
		if ($ltarget <= 0){$tm = 10}
    	else{
        	my $unround =  $cumulative * (1/($ltarget+1));
        	$tm = sprintf("%.4f", $unround);
    	}

        	if ($line =~ /^\w\w? /){
                my ($chain1,$resnum1,undef,undef,$chain2,$resnum2,undef) = split (" ",$line,7);
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
                        $rep_n ++;
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
								$output_lines .= "$res1\t$res2\n";
                            }
                            if ($rep_n == $eval){
                                $res1 = $resnum1 - $i ;
                                $res2 = $resnum2 + $i ;	
                            }
						}
                    }elsif ($rep_p > $eval) {
						$output_lines .= "$resnum1\t$resnum2\n" 
					}
					
                }else{
                    $last1 = $resnum1;
                    $last2 = $resnum2;
                    $rep_p = 1;
                    $rep_n = 1;
                }
			}
        }
		if ($tm >= $tm_threshold){print OUT "$output_lines"}
		
	}
	close OUT;
	my $input  = "$pdbcode\_$chain.kmer.$i.s";
	my $output = "$pdbcode\_$chain.kmer.$i.sorted";
	csuc($input,$output);	
}


sub csuc{
	my $input  = @_[0];
	my $output = @_[1];
	open IN, "$input" or die "not found $input\n" ;
	my @LINE ;
	while (my $line = <IN>){
		push @LINE , $line ;
	}
	@LINE = uniq @LINE ;
	@LINE = sort {$a <=> $b} @LINE;
	open OUT, ">$output";
	foreach my $lines (@LINE){print OUT "$lines"}	
}


sub suboptimal{
	my $input  = @_[0];
	my $output = @_[1];
	my @INFO ;
	open IN, "$input" or die "\n$!";
	open OUT, ">$output";
	while (my $line = <IN>){
		if ($line =~ /\S+/){
			chomp $line;
			my $warning;
			my @LINE = split (/\t/,$line);
			$LINE[0] =~ s/[A-Z]//g;
			$LINE[1] =~ s/[A-Z]//g;
			if ($LINE[5]=~ /warning/){$warning = "!"}
			my $info    = "$LINE[4]\t$LINE[0]\t$LINE[1]\t$LINE[2]\t$LINE[3]\t$LINE[4]\t$warning\n";
			push @INFO, $info;
		}
	}
	@INFO = sort {$a <=> $b} @INFO ;
	foreach my $element (@INFO){
		my @LINE = split (/\t/,$element);
		print OUT join "\t", @LINE[1..$#LINE]; 
	}
}
