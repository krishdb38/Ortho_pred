import subprocess
from multiprocessing import Pool, Process
import os, sys


Species = "./species/"
Blastp = "blastp"
blastp_matrix = "BLOSUM62"
def RunBlast(subject, parallel_num):
    subject = Species+subject       
    run_blastp =subprocess.Popen([Blastp, "-query", "./query"+"_"+str(parallel_num), "-subject", subject,
                                "-matrix", blastp_matrix, "-outfmt", "10 qseqid sseqid score length"],
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    run_blastp_stream = run_blastp.communicate()
    run_blastp_output_stream = run_blastp_stream[0]
    run_blastp_error_stream = run_blastp_stream[1]          
    return run_blastp_output_stream


def Read_Species_List(pr=0):
    species = "species/"
    """ If pr is 1, it will print "Species_List"
    Other wise only return the Value in dic Format """ 
    read_species =os.listdir(species)  #os.listdir() can be used 
    selected_species_dic = {}  #list 
    backward_selected_species_dic = {} #list
    #in Python2 selected_species_dic , and  backward_selected_species_dic is global variable 
    #number = 0
    if sys.platform == "win32":
        print("Program Running in Windows")
        for i, species in enumerate(sorted(read_species), start=1): #
            selected_species_dic[i] = species.split('\\')[-1] #in window linux /
            backward_selected_species_dic[species.split('\\')[-1]] = i 
            if pr == 1 :
                print (str(i)+".", species.split('\\')[-1])
            number = i
    else:
        print("Non - Window Platform")
        for i, species in enumerate(sorted(read_species), start=1): #
            selected_species_dic[i] = species.split('//')[-1] #in window linux /
            backward_selected_species_dic[species.split('//')[-1]] = i 
            if pr == 1 :
                print (str(i)+".", species.split('//')[-1])
            number = i
    return selected_species_dic, backward_selected_species_dic, number
