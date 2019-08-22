#               FiRES: Find Repeated Elements in Structure              #
#   A computational method for the de novo identification of repeats    #


            INSTRUCTIONS

1) Prerequisites

   1.1 DSSP

       $ sudo apt-get update
       $ sudo apt-get install dssp

   1.2 Click

       You can find tha standalone version of click here:
       http://cospi.iiserpune.ac.in/click/Download/Download.jsp


2) Installation

   2.1 Extract the compressed package
      
       $ unzip fires-master.zip

   2.2 Modify paths to the fires_v0.pl script in fires

       path_to_fires='/path/to/fires-master/scripts'

       And then:

       $ chmod 555 fires

   2.2 Define the paths to your DSSP and Click directories at 
   the top of the ./scripts/fires_v0.pl script to include the path to DSSP 
   and the path to Click:
   
   	my $path_to_click="/path/to/click";
   	my $path_to_dssp="/path/to/dssp";

   (make sure to type only the path but not the name of the program)
   
   Then:
   
       $ chmod 555 *pl

   
3) Test

       $ perl ./fires 1wm5.pdb A
       $ perl ./fires 3e3x.pdb A

4) FiRES results

   4.1 FiRES generates three output files:
       - A list of similar structural elements 
       - A structural alignment of the similar structural elements
       - A table of the best candidates, including thos that were 
	 excluded durinf the final evaluation step (marked with !)

Note: FiRES takes as input a standard DPB file. It is strongly recommended to use
      PDB files generated by the get_pdb script from Rosetta modeling software.

______________________________________________________
Please cite

DSSP:
Kabsch, W., and Sander, C. (1983) Dictionary of protein secondary structure: Pattern
recognition of hydrogen‐bonded and geometrical features. Biopolymers. 22, 2577–2637

Click:
Nguyen, M. N., and Madhusudhan, M. S. (2011) Biological insights from topology
independent comparison of protein 3D structures. Nucleic Acids Res. 10.1093/nar/gkr348
