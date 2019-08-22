#!/usr/bin/perl

use strict ;
use POSIX;
use List::MoreUtils qw(uniq);

my $infile ;
my $chain ;
my $path_to_click;
my $options     = "";

for (my $i = 0; $i < @ARGV ; $i++)   {$options    .= " $ARGV[$i] " }
if ($options =~ / -i\s+(\S+) /)      {$infile   = $1}
if ($options =~ / -c\s+(\S+) /)      {$chain   = $1}
if ($options =~ / -click\s+(\S+) /)  {$path_to_click   = $1}

open IN, "$infile" or die "\n$!" ;
while (my $line =<IN>){
       my ($query,$target)=split (" ",$line);
       my $queryn = $query ;
       $queryn =~ s/\.pdb//;
       my $targetn = $target ; 
       $targetn    =~ s/.pdb//;
       my $oldname = "$targetn-$queryn.pdb.1.clique";
       my $newname = "$queryn-$targetn.pdb.1.clique";
       open CHECK, "$oldname" ;
       while (my $ch = <CHECK>){
               if ( $ch =~ /\S+/ ){ ; next}
       }
       system ("$path_to_click/click $target $query -a 1 >> click-file");
       system ("cp $oldname $newname");
}


__END__
