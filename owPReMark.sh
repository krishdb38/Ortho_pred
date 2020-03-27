#!/bin/sh 
python3 owPReMark.py -M 1 3  -g AAE CAC ECO ECU HIN   -m BLOSUM45 -c 8

echo "all tasks Completed"

# M ==> Mode 1.Blastp 2 BLASTP Using precalculated data and 3 is Clustering 
# g ==> genome. Write the name of the genome to Use
# m ==> Which Matrix to use default is BLOSUM62. BLOSUM45 and BLOSUM80 can be used
# c ==> Number of Cpu to Use not greater than 50
# i ==> inflation Factor ( default is 1.4 )
# s ==> path for species default is ok ("./species./")
# F ==> score_File(default is "./score_file"
# B ==> Path for Precalculated blastp data
# l ==> CLUSTER_OUT

# Availiable list of Species
# AAE SCE   SPO  SYN   YPE CAC   ECO  HIN   LLA SPY  TMA 

# chicken.faa           dog.faa          fruitfly.faa       human.faa 
#pufferfish.faa         seasquirt.faa    zebrafish.faa      mouse.faa
#chimpanzee.faa         elegans.faa     opossum.faa         rat.faa


# BlastP done 
#AAE CAC HIN SCE 
#