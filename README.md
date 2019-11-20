#               FiRES: A computational method for *de novo* identification of proteins with repeated elements    	#

FiRES (Find Repeated Elements in Structure) relies on a topology-independent structure alignment method to identify repeating elements in protein structure.
FiRES can be used for the discovery of proteins containing similar structural elements with very low
sequence identity (less than 20%)

## Benchmark:
https://github.com/Claualvarez/Internal_structure_similarity_benchmark

## Web server: 
http://fires.ifc.unam.mx/

## Standalone program:
INSTRUCTIONS 
(Linux)

1) Prerequisites

   1.1 List::MoreUtils module for perl
   
       $ cpanm List::MoreUtils
   
   1.2 DSSP
   
     The instructions are available here: https://swift.cmbi.umcn.nl/gv/dssp/
     
     or try a package manager:
 
       $ sudo apt-get update
       $ sudo apt-get install dssp

     or
     
       $ conda install -c salilab dssp

   1.3 Click

      The standalone version of click is available here:
       http://cospi.iiserpune.ac.in/click/Download/Download.jsp


2) Setup

   2.1 Extract the compressed package
      
       $ unzip fires-master.zip
       $ cd fires-master

   2.2 Modify paths to the scripts' directory in the fires controller script:

       path_to_fires='/path/to/fires-master/scripts'

   2.3 Define the paths to your DSSP and Click directories at 
   the top of the fires controller script:
   
        $path_to_click="/path/to/click";
        $path_to_dssp="/path/to/dssp";

   (make sure to type only the path but not the name of the program)

   
3) Test

       $ bash ./fires 1wm5.pdb A
       $ bash ./fires 3e3x.pdb A
       
   To save the top non-overlapping results run:
   
       $ bash ./fires [pdb_file] [chain_name] > [output_name]
       
   e.g.
   
       $ bash ./fires 1wm5.pdb A > 1wm5.A.fires_output.txt

4) FiRES results

   4.1 FiRES generates three output files:
       a) A list of similar structural elements
       b) A structural alignment of the similar structural elements
       c) A table of the best candidates, including those that were 
	 excluded during the final evaluation step (marked with !)
	 

Note: FiRES takes as input a standard DPB file. It is strongly recommended to use
      PDB files generated by the get_pdb script from Rosetta modeling software.

______________________________________________________
## Please cite

DSSP:
Kabsch, W., and Sander, C. (1983) Dictionary of protein secondary structure: Pattern
recognition of hydrogen‐bonded and geometrical features. Biopolymers. 22, 2577–2637

Click:
Nguyen, M. N., and Madhusudhan, M. S. (2011) Biological insights from topology
independent comparison of protein 3D structures. Nucleic Acids Res. 10.1093/nar/gkr348
