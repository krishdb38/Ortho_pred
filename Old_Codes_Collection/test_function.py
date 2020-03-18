import glob
import sys
import subprocess
import pdb
import time
from progress.bar import Bar
import multiprocessing
import queue
import numpy.matlib
import numpy
import copy
import operator
import pprint
from itertools import *
import argparse
import datetime

Species = "./species/"
Blastp = "blastp"
blastp_matrix = "BLOSUM62"

def RunBlast(subject, parallel_num):
    import subprocess
    "By this Function it will create a Pipe line to run Blastp in Computer by input Parameter"
    subject = Species+subject 
    cmd = [Blastp, "-query", "./query"+"_"+str(parallel_num), "-subject", subject,
                                "-matrix", blastp_matrix, "-outfmt", "10 qseqid sseqid score length"]  

    #run_blastp =subprocess.Popen(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    run_blastp = subprocess.run(cmd,shell=True ,capture_output = True,text = True)
    #output Format 10 qseqid query (e.g. gene sequence id  ,  sseqid subject (e.g. reference genome) genome id)
    run_blastp_stream = run_blastp.stdout()

    #run_blastp_stream = run_blastp.communicate()
    #run_blastp_output_stream = run_blastp_stream[0]
    #run_blastp_error_stream = run_blastp_stream[1] 
    #          
    return run_blastp_stream
a= RunBlast("AAE","AAE")
print(a)