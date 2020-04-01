#!/usr/bin/python3
import sys
import os                  # For system operation
import glob                #list the file unix system
import subprocess          #to run blastp from OS
import time                # to display time
from alive_progress import alive_bar #For Progress Bar
import multiprocessing  #For Parallel Processing
import queue            # For Multiprocessing
import numpy.matlib     # Matrix operation 
import numpy as np      # For Mathmatical (Algebra) Operation
import copy             # To copy file
import operator   #to concate the list we will remove this later 
import argparse
import datetime         #to display current time & calculate difference
from Bio import SeqIO  # For Bio Python Sequence Object

__author__ = " SS Kim & Adhikari Krish"
__teammates__ = [" pk", "okg"]
__copyright__ = "  Copyright 2020"
__credits__ = "  THER................."
__license__ = "  GNU"
__version__ = "  1.00"
__maintainer__ = "  Adhikari Krish"
__status__ = "  Ortholog Detection"
__email__ = "  *********@gmail.com"

print("\n\n\n")
print("="*82)
print("||\tAuthor     : %s\t\t\t\t\t||" % __author__)
print("||\tTeammates  : %s\t\t\t\t\t\t\t||" % ",".join(__teammates__))
print("||\tCopyright  :%s \t\t\t\t\t\t||" % __copyright__)
print("||\tCredits    :%s \t\t\t\t\t||" % __credits__)
print("||\tLicense    :%s \t\t\t\t\t\t\t||" % __license__)
print("||\tVersion    :%s \t\t\t\t\t\t\t||" % __version__)
print("||\tMaintainer :%s \t\t\t\t\t\t||" % __maintainer__)
print("||\tEmail      :%s \t\t\t\t\t||" % __email__)
print("||\tStatus     :%s \t\t\t\t\t||" % __status__)
print("="*82)


parser = argparse.ArgumentParser(description='This program help you to detect ortholog between the protein sequences from different genomes\
     and to cluster orthologs to ortholog groups.',
                                 add_help=True, prefix_chars='-+')
essential_group = parser.add_argument_group('Essential arguments')
blastp_group = parser.add_argument_group('blastp optional arguments')
mcl_group = parser.add_argument_group('clustering optional arguments')
essential_group.add_argument('-M', '--mode',
                             action='store',
                             nargs='+',
                             dest='mode',
                             choices=['1', '2', '3'],
                             help="1 : BLASTP. 2 : BLASTP using precalculated data. 3 : clustering. ex) 1 or 2 or 1 3 or 2 3")
essential_group.add_argument('-g', '--genomes',
                             nargs='+',
                             dest='genomes',
                             help="Select a genomes to detect orthologs.")
parser.add_argument('-c', '--cpu',
                    action='store',
                    default='1',
                    dest='cpu_count',
                    type=int,
                    help='Set the number of CPU to use in program. The default is "1". This system has '+str(multiprocessing.cpu_count())+" CPUs")
blastp_group.add_argument('+r',
                          action='store_true',
                          default=True, #Default was false
                          dest='save_raw_blastp_score',
                          help='Save the raw score of blastp to the score_file dirctory.')
blastp_group.add_argument('-t', '--threshold',
                          action='store',
                          default=5,
                          type=int,
                          dest='threshold_score',
                          help='Set the threshold score. The threshold score is an allowable range to be an ortholog.' +
                          'If the score difference between backward best hit score and forward best hit score is less than threshold score, it is considered the forward best hit pair to be the ortholog.' +
                          'That is  called the "one-way threshold best hit" by us. The default value is "5". "backward best hit score - forward best hit score <= threshold score"')
blastp_group.add_argument('-m', '--matrix',
                          action='store',
                          default='BLOSUM62',
                          dest='blastp_matrix',
                          choices=['BLOSUM45', 'BLOSUM62', 'BLOSUM82'],
                          help='select the matrix to use in blastp. The default value is "BLOSUM62".')
blastp_group.add_argument('-s', '--species',
                          action='store',
                          default='./species/',
                          dest='Species',
                          help='set the path of species directory. The default is  "./species/".')
blastp_group.add_argument('-b', '--blastp',
                          action='store',
                          default="blastp",  # For Windows
                          dest='Blastp',
                          help='set the path of blastp file to run the blastp program. The default is "blastp".')
blastp_group.add_argument('-F', '--scorefile',
                          action='store',
                          default='./score_file/',
                          dest='Score_file',
                          help='The default path is "./score_file/".')
blastp_group.add_argument('-B', '--blastp_data',
                          action='store',
                          default='./blastp_data/',
                          dest='Blastp_data',
                          help='set the path of precalculated blastp data. The default path is "./blastp_data/".')
mcl_group.add_argument('-i', '--IF',
                       action='store',
                       default=1.4,
                       dest='inflation_factor',
                       type=float,
                       help='The default value is "1.4".')
mcl_group.add_argument('-l', '--loop',
                             action='store',
                             default=60,
                             type=int,
                             dest='infinite_loop',
                             help='prevent the infinite loop of mcl algorithm. The default value is "60".')
mcl_group.add_argument('-o', '--out',
                             action='store',
                             default='./cluster_out',
                             dest='Cluster_out',
                       help='set the path and name of ortholog cluster file, log and error_warning file. The default path and name is "./cluster_out".')
mcl_group.add_argument('+v',
                       action='store_true',
                       default=False,
                       dest='verbose',
                       help='verbosely show  information of a big matrix computaion.')
parser.add_argument('--version',
                    action='version',
                    version='%(prog)s Ver. 11.1')
command_options = parser.parse_args()


def matrix_name():
    """This Function will Return the Matrix name choosed by User.
    BLOSUM45 ,  BLOSUM62 , BLOSUM82 """
    print("BLOcks SUbstitution Matrix (BLOSUM) is a Substitution matrix used for sequence alignment of Proteins")
    print("""\n1. BLOSUM45 :-For more distantly related Proteins alignment DataBase\n2. BLOSUM62 :- MidRange Seq with more than 62%similarity\
         \n3. BLOSUM82 :- More related Proteins\nOther Keys will exit the Program""")
    metrix_num = input("\nEnter a matrix number: ")
    if metrix_num not in ("1", "2", '3'):
        print("Wrong input *%s*Sorry not in list\n" %
              metrix_num, "*"*20, "Good Bye", "*"*20, "\n")
        sys.exit(1)
    if metrix_num == "1":
        return "BLOSUM45"
    if metrix_num == "2":
        return "BLOSUM62"  # 62 is not availiable currently
    if metrix_num == "3":
        return "BLOSUM82"


def query_sequence(genome):
    """This Function Read Fastaq Files and return as a list Format with gene Position as a index
    This Function create a list of Seprate Query Sequence"""
    # This Function is created by using Bio Python we will test later
    gene_seq = ""
    gene_seq_list = []
    try:
        with open(Species+genome) as gene:
            for each_line in gene:
                if ">" in each_line:
                    if gene_seq != "":
                        gene_seq_list.append(gene_seq)
                        gene_seq = ""
                gene_seq = gene_seq+each_line
            gene_seq_list.append(gene_seq)
            #print("query_sequence() Run Successfully")
            return gene_seq_list

    except IOError as err:
        print("IOError occurred in query_sequence function : " + str(err))


def query_sequence_(genome):
    """This Function Read Fasta (Genome) file and return as a list Format with gene Position and Sequence Developed by Krish"""
    global Species
    try:
        return [str((seq_record.id+"\n"+seq_record.seq)) for seq_record in SeqIO.parse(Species+genome, "fasta")]
        # To understand this Function Bio Python library needs to be studied
    except IOError as err:
        print(str(err))


def write_query(query, parallel_num):
    "This Function Write Query with file Name query+ parallel_num in same directory and raise IO error if Error rises"
    # print("write_query(query,parallel_num)")
    try:
        with open("./query/query_"+str(parallel_num), "w") as write_query:
            write_query.write(query)
    except IOError as err:
        print("IOError occurred in write_query function : " + str(err))


def run_blast(subject, parallel_num):
    """By this Function it will create a Pipe line to run Blastp in Computer by input Parameter Subject is whole Genome
         and parallel_number is query file Created by early step.to run a big file it is time Consuming. So reduce a File size and run
         The Output Format 10 --> Comma Separated Values 
         qseqid --> Query Seq ID
         ssequid --> Subject of Seq. id ID
         Score -->Raw Score
         length-->Alignment length"""

    subject = Species+subject  # later remove this
    cmd = ["blastp", "-query", "./query/query_"+str(parallel_num), "-subject", subject,
           "-matrix", blastp_matrix, "-outfmt", "10 qseqid sseqid score length"]

    # query subject is inside query Folder
    #run_blastp =subprocess.Popen(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #run_blastp = subprocess.run(cmd,shell=True ,capture_output = True,text = True)
    # output Format 10 qseqid query (e.g. gene sequence id  ,  sseqid subject (e.g. reference genome) genome id)
    run_blastp = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #run_blastp_stream = run_blastp.stdout()
    run_blastp_stream = run_blastp.communicate()
    run_blastp_output_stream = run_blastp_stream[0]
    #run_blastp_error_stream = run_blastp_stream[1]

    # Later We can return by this
    # return run_blastp.communicate[0]
    return run_blastp_output_stream


def same_species_forward_best_hit(blastp_score):
    """Search the forward best hit among the blastp scores of same species.
    Because there are an duplicated genes in a same genome.
    blstp_score file is return value of run_blast() Function
    When the blastp score compare with blastp score of duplicate gene, if score and length are same, blasp score of duplicated gene is added to a second best score."""
    #print(blastp_score)
    blastp_score_split_list = []
    temp_best_score = ['-1', '-1', '-1']
    second_temp_best_score = []
    best_score = []
    second_best_score = []
    blastp_score_split = blastp_score.split("\n")
    # delete of [''] in the last index  ex) ['gi,gi,1,1','gi,gi,2,2','']
    del blastp_score_split[-1]
    for i in blastp_score_split:
        blastp_score_element = i.split(',')
        blastp_score_split_list.append(blastp_score_element)
    # ex) k is ['gi|15605613|ref|NP_212986.1|', 'gi|15605613|ref|NP_212986.1|', '3702', '699']
    for k in blastp_score_split_list:
        if k[0] == k[1]:  # if the Position (gene) is Same
            best_score.append(k)
        elif k[0] != k[1]:
            if int(k[2]) > int(temp_best_score[2]):  # Compare score
                temp_best_score = k
            elif int(k[2]) == int(temp_best_score[2]):  # Is Equal
                if int(k[3]) > int(temp_best_score[3]):  # compare length
                    temp_best_score = k
                elif int(k[3]) == int(temp_best_score[3]):
                    second_temp_best_score.append(k)
 #    print ("############ temp best score ############", temp_best_score)
    second_best_score.append(temp_best_score)
    for j in second_temp_best_score:
        if j[2] == temp_best_score[2] and j[3] == temp_best_score[3]:
            second_best_score.append(j)
    for m in second_best_score:
        if (best_score[0][2] == m[2] and int(best_score[0][3]) <= int(m[3])) or int(best_score[0][2]) < int(m[2]):
             # '104' < '23' is True because of string. So the int function is used.
            best_score.append(m)
    #print("best Score is", best_score)
    return best_score


def forward_best_hit(blastp_score):
    """Search the forward best hit among the blastp scores of same species."""
    #print("Rnunning GetForward BestHit")
    blastp_score_split_list = []
    temp_best_score = ['-1', '-1', '-1']
    second_temp_best_score = []
    best_score = []
    blastp_score_split = blastp_score.split("\n")
    # delete of ['']   ex) ['gi,gi,1,1','gi,gi,2,2','']
    del blastp_score_split[-1]
    for i in blastp_score_split:
        blastp_score_element = i.split(',')
        blastp_score_split_list.append(blastp_score_element)

    for k in blastp_score_split_list:
        # ex) k is ['gi|15605613|ref|NP_212986.1|', 'gi|15605613|ref|NP_212986.1|', '3702', '699']
        # print ">>>>>>>>>>>>>>>forward_best_hit   k", k
        if int(k[2]) > int(temp_best_score[2]):  # Compare the Score of the Blast
            temp_best_score = k
        elif int(k[2]) == int(temp_best_score[2]):
            if int(k[3]) > int(temp_best_score[3]):  # Compare the Lenth of Sequence
                temp_best_score = k
            elif int(k[3]) == int(temp_best_score[3]):
                second_temp_best_score.append(k)
  #print "############ temp best score ############", temp_best_score
    best_score.append(temp_best_score)
    for j in second_temp_best_score:
        if j[2] == temp_best_score[2] and j[3] == temp_best_score[3]:
            best_score.append(j)
    return best_score, blastp_score_split_list


def division_parallel_query(query_v, query_division_value, cpu_count, query_v_len):
    "query_v is list of queries, query_division_value is , cpu_count is the Number of CPU User Entered and query_v_len is the length of the query"
    parallel_query = []
    parallel_query_start = 0
    if query_v_len % cpu_count == 0:  # perfect division
        for i in range(cpu_count):
            i += 1
            if parallel_query_start == 0:
                parallel_query.append(
                    query_v[int(parallel_query_start):i*int(query_division_value)])
                parallel_query_start += 1
            else:
                parallel_query.append(query_v[int(
                    parallel_query_start)*int(query_division_value):i*int(query_division_value)])
                parallel_query_start += 1
    else:  # imperfect division
        for i in range(cpu_count):
            i += 1
            if parallel_query_start == 0:
                parallel_query.append(
                    query_v[int(parallel_query_start):i*int(query_division_value)])
                parallel_query_start += 1
            elif i < cpu_count:
                parallel_query.append(query_v[int(
                    parallel_query_start)*int(query_division_value):i*int(query_division_value)])
                parallel_query_start += 1
            elif i == cpu_count:
                parallel_query.append(
                    query_v[int(parallel_query_start)*int(query_division_value):])
                parallel_query_start += 1

    return parallel_query


def run_parallel_query(species_of_query, species_of_subject, query_v, parallel_num):
    """ run_parallel_query(i ,k , query_v , cpu_count) i and k are user selected number in a list Format
    Run the following functions. write_query, run_blast, same_species_forward_best_hit, forward_best_hit
    Save the files which are oneway_threshold_best_hit, second_oneway_threshold_best_hit, 
    blastp_score_split_list and raw_blastp_score (optional) by each species. 
    parallel_num is the Number of CPU selected by the user if CPu 1 selected then 1 """

    #print("run_parallel_query Running")
    global selected_number, selected_species_dic
    # bar = Bar('Processing '+str(parallel_num), max = len(query_v)) #progressing bar setting , Creating a Object
    # bar is not Supported in Python 3    We use alive_bar instead
    with alive_bar(len(query_v)) as bar:  # declare your set of items for loop
        for j in query_v:
            # bar.next() #progressing bar not supported
            bar()# Call after Consuming One Item
            write_query(j, parallel_num)  # if 1 Added Here Add also to Run
        # This Function Only Write a file with j name and parallel_num i.e CPU Count
            blastp_score = run_blast(
                selected_species_dic[species_of_subject], parallel_num)
        # " Return Byte File type, and this is only for one gene position query and Whole Genome of species"
            if blastp_score != '':  # if blastp run successfully
                best_score, blastp_score_split_list = forward_best_hit(
                    blastp_score.decode())  # .decode() Convert in to str format
                # ex) AAE == AAE. It will save best_score without reversing run_blast.
                if species_of_query == species_of_subject:
                    same_species_forward_best_score = same_species_forward_best_hit(
                        blastp_score.decode())
                    for best_score_element in same_species_forward_best_score:
                        # ex) [A1 of AAE, A1 of AAE, 30]
                        if best_score_element[0] == best_score_element[1]:
                            with open(Score_file+selected_species_dic[species_of_query] + "_" +
                                      selected_species_dic[species_of_subject]
                                      + "_oneway_threshold_best_hit_Score"
                                      + str(threshold_score), "a") as oneway_threshold_best_hit:
                                save_best_score = selected_species_dic[species_of_query]+"_"+best_score_element[0].split("\s")[0]\
                                    + " "+selected_species_dic[species_of_subject]+"_"+best_score_element[1].split(
                                        "\s")[0]+" "+best_score_element[2]+"\n"
                                # best_score_element[0].split("|") ==> ['gi', '15642790', 'ref', 'NP_227831.1', '']
                                oneway_threshold_best_hit.write(
                                    save_best_score)
                        else:  # ex) [A1 of AAE, A2 of AAE, 30]
                            with open(Score_file+selected_species_dic[species_of_query]+"_"
                                      + selected_species_dic[species_of_subject]
                                      + "_second_oneway_threshold_best_hit_Score"
                                      + str(threshold_score), "a") as second_oneway_threshold_best_hit:
                                second_save_best_score = selected_species_dic[species_of_query]\
                                    + "_"+best_score_element[0].split("\s")[0]\
                                    + " "+selected_species_dic[species_of_subject]+"_"\
                                    + best_score_element[1].split("\s")[0]+" "+best_score_element[2]+"\n"
                                second_oneway_threshold_best_hit.write(
                                    second_save_best_score)
                else:  # If species_of_query not equal with species_of_subject, run reversing run_blast
                    for best_score_element in best_score:
                        if not '-1' in best_score_element:
                            with open(Score_file+selected_species_dic[species_of_query]
                                      + "_" +
                                      selected_species_dic[species_of_subject] +
                                      "_"+"best_score_S"
                                      + str(threshold_score)+"_"+str(parallel_num), "a") as save_best_hit:
                                best_score_save = selected_species_dic[species_of_query]+"_"\
                                    + best_score_element[0].split("\s")[0]+" "+selected_species_dic[species_of_subject]\
                                    + "_"+best_score_element[1].split(
                                        "\s")[0]+" "+best_score_element[2]+" "+best_score_element[3]+"\n"
                                save_best_hit.write(best_score_save)

                    for blastp_score_split_list_element in blastp_score_split_list:
                        with open(Score_file+selected_species_dic[species_of_query]+"_"
                                  + selected_species_dic[species_of_subject]+"_"+"blastp_score_split_list_S"
                                  + str(threshold_score)+"_"+str(parallel_num), "a") as save_blastp_score_split_list:
                            blastp_score_split_list_save = selected_species_dic[species_of_query]+"_"\
                                + blastp_score_split_list_element[0].split("\s")[0]\
                                + " "+selected_species_dic[species_of_subject]+"_"+blastp_score_split_list_element[1].split("\s")[0]\
                                + " " + \
                                blastp_score_split_list_element[2]+" " \
                                + blastp_score_split_list_element[3]+"\n"
                            save_blastp_score_split_list.write(
                                blastp_score_split_list_save)

            else:
                print("No value in blastp_score")
            if save_raw_blastp_score:
                with open(Score_file+selected_species_dic[species_of_query]
                          + "_"+selected_species_dic[species_of_subject]+"_S"
                          + str(threshold_score)+"_"+str(parallel_num), "a") as save_blastp:
                    save_blastp.write(blastp_score.decode()) # Only string is format is supported
    # bar.finish() # progressing bar finish
    return  # None


def oneway_threshold_best_hit(mode):
    """ This Function accept the mode and Run program and return backward_best_hit_work_list mode 1 and 2 is supported """
    #print("oneway_threshold_best_hit running")
    global user_selected_number, cpu_count,precalculated_data_list,new_calculated_data_list
    b_info = "Running the blastp & forward best hit searches "
    process_list = []
    backward_best_hit_work_list = []
    if "1" in mode:
        """We have 3 Mode 1 is for blastp, Mode 2 is for BLASTP using precalcualted data and Mode 3 is for clustering"""
        for i in user_selected_number:  # Select species to write query 
            query_v = query_sequence(selected_species_dic[i])
            # query_v is a list Format  with a position gene id , Seq
            # User Selected  [1, 3, 5] is list of User input
            for k in user_selected_number:
                if k < i:  # gene ====> query 1->1 1->2 1->3 2->2 2->3  Forward checking will skip the same or less
                    continue
                else:
                    print(b_info+" between %s genome and %s genome"
                          % (selected_species_dic[i], selected_species_dic[k]))
                    # this will Create a loop though 2 times
                    # length of Genome gene ID in Sequence File
                    query_v_len = len(query_v)
                    if cpu_count == 1:
                        #"No Parallel Computing while cpu count == 1"
                        blastp_time_start = time.time()
                        run_parallel_query(i, k, query_v, cpu_count)
                        #"i is first species k is second species number queryv is list file of i Position genome"
                        blastp_time_end = time.time()
                        print(b_info+ "took %.2f minutes" % (
                            (blastp_time_end-blastp_time_start)/60))
                    else:
                        # If the number of query_v_len is less than cpu_count, Remark will select the number of query_v_len.
                        if query_v_len < cpu_count:
                            #"because the cpu_count will seprate list items from query_v file.In general this will not happened because length of query is always long"
                            blastp_time_start = time.time()
                            # 1 is query_division_value. Because query_v_len / query_v_len(=cpu_count) is 1.
                            parallel_query = division_parallel_query(
                                query_v, 1, query_v_len, query_v_len)
                            for m in range(query_v_len):
                                process = multiprocessing.Process(
                                    target=run_parallel_query, args=(i, k, parallel_query[m], m+1))
                                # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                process_list.append(process)
                                process.start()
                            for n in process_list:
                                n.join()
                            blastp_time_end = time.time()
                            print(b_info+"took %.2f minutes " % (
                                (blastp_time_end-blastp_time_start)/60))
                        else:
                            #"this wiil run if cpu_count is more than 1"
                            blastp_time_start = time.time()
                            query_division_value = query_v_len / cpu_count
                            parallel_query = division_parallel_query(
                                query_v, query_division_value, cpu_count, query_v_len)
                            for m in range(cpu_count):
                                process = multiprocessing.Process(
                                    target=run_parallel_query, args=(i, k, parallel_query[m], m+1))
                                # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                process_list.append(process)
                                process.start()
                            for n in process_list:
                                n.join()
                            blastp_time_end = time.time()
                            print(b_info+"took %.2f minutes " % (
                                (blastp_time_end-blastp_time_start)/60))
                if not i == k:
                    backward_best_hit_work_list.append((i, k, query_v_len))
                    print("Backward")
    elif "2" in mode:
        
        """If blast Files exist then this will pass the blast time other wise blast will run.
        This will reduce the time for blast and increase the speed of the system.
        """
        print("\nMode 2, Running from precalculated data sets which will save the Time")
        print(len(precalculated_data_list),"oneway_threshold_best_hit_ files are located") # !!! Delete This row later
        for i in user_selected_number:  # Select species to write query
            query_v = query_sequence(selected_species_dic[i])
            for k in user_selected_number:  # Select of subject
                if Blastp_data+selected_species_dic[i]+"_"+selected_species_dic[k]\
                    + "_oneway_threshold_best_hit_Score"\
                        + str(threshold_score) in precalculated_data_list: 
                    print("Skipped")
                    used_precalculated_data_list.append(selected_species_dic[i]+"_"+selected_species_dic[k])
                    #used_precalcualted _data_list is created before calling this function
                    #if files exist then passed
                    continue #File will skkipped
                else:
                    if k < i:  # gene ====> query 1->1 1->2 1->3 2->2 2->3
                        continue
                    else:
                        "same as Mode 1, if Precalculated file doesnot exist"
                        print(b_info+ " between %s genome and %s genome" % (
                            selected_species_dic[i], selected_species_dic[k]))
                        query_v_len = len(query_v)
                        if cpu_count == 1:
                            blastp_time_start = time.time()
                            run_parallel_query(i, k, query_v, cpu_count)
                            blastp_time_end = time.time()
                            print(b_info+ "took %.2f minutes" % (
                                (blastp_time_end-blastp_time_start)/60))
                        else:
                            # If the number of query_v_len is less than cpu_count, Remark will select the number of query_v_len.
                            if query_v_len < cpu_count:
                                blastp_time_start = time.time()
                                # 1 is query_division_value. Because query_v_len / query_v_len(=cpu_count) is 1.
                                parallel_query = division_parallel_query(
                                    query_v, 1, query_v_len, query_v_len)
                                for m in range(query_v_len):
                                    process = multiprocessing.Process(
                                        target=run_parallel_query,
                                        args=(i, k, parallel_query[m], m+1))
                                    # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                    process_list.append(process)
                                    process.start()
                                for n in process_list:
                                    n.join()
                                blastp_time_end = time.time()
                                print(b_info+" took %.2f minutes" % (
                                    (blastp_time_end-blastp_time_start)/60))
                            else:
                                blastp_time_start = time.time()
                                query_division_value = query_v_len / cpu_count
                                parallel_query = division_parallel_query(
                                    query_v, query_division_value, cpu_count, query_v_len)
                                for m in range(cpu_count):
                                    process = multiprocessing.Process(
                                        target=run_parallel_query,
                                        args=(i, k, parallel_query[m], m+1))
                                    # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                    process_list.append(process)
                                    process.start()
                                for n in process_list:
                                    n.join()
                                blastp_time_end = time.time()
                                print(b_info+ " took %.2f minutes" % (
                                    (blastp_time_end-blastp_time_start)/60))
                        new_calculated_data_list.append(
                            selected_species_dic[i]+"_"+selected_species_dic[k])
                        if not i == k:
                            backward_best_hit_work_list.append(
                                (i, k, query_v_len))
    return backward_best_hit_work_list


def backward_best_hit(args):
    bck_info = "Running the backward_best_hit"
    test =" "
    print(bck_info)
    species_of_query, species_of_subject, query_v_len = args
    start_time_bbh = time.time()
    forward_best_hit_score_list = []
    blastp_score_split_list = []
    print(bck_info+ " between %s genome %s genome" % (
        selected_species_dic[species_of_query], selected_species_dic[species_of_subject]))
    # If the number of query_v_len is less than cpu_count, the cpu_count is changed to query_v_len.
    if query_v_len < cpu_count:
        for parallel_num in range(query_v_len):
            parallel_num += 1
            with open(Score_file+selected_species_dic[species_of_query]+"_"
                      + selected_species_dic[species_of_subject] + "_"+"best_score_S"
                      + str(threshold_score)+"_"+str(parallel_num), "r") as best_hit_score:
                for each_line in best_hit_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
                    forward_best_hit_score_list.append(split_each_line)
            with open(Score_file+selected_species_dic[species_of_query]+"_"
                      + selected_species_dic[species_of_subject] +
                      "_"+"blastp_score_split_list_S"
                      + str(threshold_score)+"_"+str(parallel_num), 'r') as blastp_score:
                for each_line in blastp_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
                    blastp_score_split_list.append(split_each_line)
    else:
        for parallel_num in range(cpu_count):
            parallel_num += 1
            with open(Score_file+selected_species_dic[species_of_query]+"_"
                      + selected_species_dic[species_of_subject]+"_"+"best_score_S"
                      + str(threshold_score)+"_"+str(parallel_num), "r") as best_hit_score:
                for each_line in best_hit_score:
                    split_each_line = each_line.split(" ")

                    split_each_line[3] = int(split_each_line[3])
                    forward_best_hit_score_list.append(split_each_line)
            with open(Score_file+selected_species_dic[species_of_query] + "_"
                      + selected_species_dic[species_of_subject]
                      + "_"+"blastp_score_split_list_S"
                      + str(threshold_score)+"_"+str(parallel_num), 'r') as blastp_score:
                for each_line in blastp_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
                    blastp_score_split_list.append(split_each_line)

    # "bar = Bar("Searching : "+selected_species_dic[species_of_query]+"-"+selected_species_dic[species_of_subject], max = len(forward_best_hit_score_list))" #remove ""

    for forward_best_hit_score_element in forward_best_hit_score_list:
        matching_list = []
        backward_best_score = ['-1', '-1', '-1']
        # bar.next()
        for element in blastp_score_split_list:
            if element[1] == forward_best_hit_score_element[1]:
                matching_list.append(element)

        for element in matching_list:
            if int(element[2]) > int(backward_best_score[2]):
                backward_best_score = element
    #        with open('./'+selected_species_dic[species_of_query]\
    # +"_"+selected_species_dic[species_of_subject]+'_subtraction'+"_"+str(threshold_score), 'a') as subtraction :
    #            save_data = int(backward_best_score[2]) - int(forward_best_hit_score_element[2])
    #            subtraction.write(str(save_data)+"\n")

        if int(backward_best_score[2]) - int(forward_best_hit_score_element[2]) <= threshold_score:
            with open(Score_file+selected_species_dic[species_of_query]
                      + "_"+selected_species_dic[species_of_subject]
                      + "_oneway_threshold_best_hit_Score"+str(threshold_score), "a")\
                    as other_oneway_threshold_best_hit:
                save_data = forward_best_hit_score_element[0]\
                    + " "+forward_best_hit_score_element[1]\
                    + " " + str(int(forward_best_hit_score_element[2]))+"\n"
                other_oneway_threshold_best_hit.write(save_data)

    # bar.finish()
    finish_time_bbh = time.time()
    rbh_time = float((finish_time_bbh - start_time_bbh)/60)
    print(bck_info+ " %s-%s took %.2f minutes" %
          (selected_species_dic[species_of_query], selected_species_dic[species_of_subject], rbh_time))
    return rbh_time


def search_equal_bbh_data(target_A):
    """Search the equal backward best hit data. ex) AAE_AAE_backward_best_hit
    The if with pass statement will be later updated """
    # print("search_equal_bbh_data")
    put_data = equal_BBH_data_dic[target_A]
    if put_data[1] == 0:
        pass
    else:
        copy_put_data = copy.copy(put_data)
        copy_put_data.insert(0, target_A)
        results.put(copy_put_data)
        equal_BBH_data_dic[target_A][1] = 0
        for i in second_equal_BBH_data:
            if i[2] == 0:
                pass
            else:
                if i[0] == target_A or i[1] == target_A:
                    copy_second_put_data = copy.copy(i)
                    # Don't put results as queue. Because the tasks will put copy_second_put_data to results as queue.
                    tasks.put(copy_second_put_data)
                    i[2] = 0


def search_unequal_bbh_data(target_b):
    """Search the unequal backward best hit data. ex) AAE_CAC_backward_best_hit
    No Return !! vlaues means Null return"""
    global unequal_BBH_data
    print(unequal_BBH_data)
    for i in unequal_BBH_data:
        if i[2] == 0:
            " "
            pass
        else:
            if target_b[0] == i[0] or target_b[0] == i[1] or target_b[1] == i[0] or target_b[1] == i[1]:
                #sample of B  - ["AAE_gi|156|ref....","CAC_gi|1560..",85]
                copy_i = copy.copy(i)
                tasks.put(copy_i)
                unequal_BBH_data[unequal_BBH_data.index(i)][2] = 0


def matching_bbh(target):
    """ Match the backward best hit """
    #print("running matching_bbh")
    if target[2] == 0:
        return

    else:
        copy_target = copy.copy(target)
        search_equal_bbh_data(copy_target[0])
        search_equal_bbh_data(copy_target[1])
        results.put(copy_target)
        unequal_BBH_data[unequal_BBH_data.index(target)][2] = 0

    for j in unequal_BBH_data:
        if j[2] == 0:
            pass

        else:
            if copy_target[0] == j[0] or copy_target[0] == j[1] or copy_target[1] == j[0] or copy_target[1] == j[1]:
                copy_j = copy.copy(j)
                unequal_BBH_data[unequal_BBH_data.index(j)][2] = 0
                tasks.put(copy_j)

    while not tasks.empty():
        get_task = tasks.get()
        search_equal_bbh_data(get_task[0])
        search_equal_bbh_data(get_task[1])
        results.put(get_task)
        search_unequal_bbh_data(get_task)


def matrix_clustering_ortholog(element_set):
    # matrix_clustering_ortholog(element_set, bar): # this is Old Method old is generating is  removed
    """ Generate the matrix of clustering ortholog. """
    print("matrix_clustering_ortholog Function is running element_set is Passed is ")
    print(element_set)
    row_data = []
    col_data = []
    temp_results = queue.Queue()
    "Create a queue object with a given maximum size. If maxsize is <=0, teh queue size is infinite"
    # bar.next()
    for element in element_set:
        # if element[0] exist, returning the index in the row_data.
        if row_data.count(element[0]) > 0:
            # element[0] is data of row. ['gi|15606057|ref|NP_213434.1|', 'gi|15606057|ref|NP_213434.1|', '3823\n']
            row = row_data.index(element[0])

        else:
            row = len(row_data)
            row_data.append(element[0])
            # if element[0] doesn't exist, appending the element[0] to the col_data.
            if col_data.count(element[0]) < 1:
                col_data.append(element[0])

        if col_data.count(element[1]) > 0:
            col = col_data.index(element[1])  # element[1] is data of col.

        else:
            col = len(col_data)
            col_data.append(element[1])
            if row_data.count(element[1]) < 1:
                col_data.append(element[1])

        temp_results.put([row, col, element[2]])

    # create a new matrix of given shape(the size_resuls) and type, filled with zeros.
    score_matrix = numpy.matlib.zeros(
        (len(row_data), len(col_data)), dtype=np.float)
    # np.zeros() can be used, will test later,
    while not temp_results.empty():
        get_temp_results = temp_results.get()
        row = get_temp_results[0]
        col = get_temp_results[1]
        score_matrix[row, col] = get_temp_results[2]
        score_matrix[col, row] = get_temp_results[2]

    # If the elements of row and col is less than 2, it is excluded.
    if len(row_data)*len(col_data) > 4:
        # The big size of matrix(bigger than 1000 X 1000) will be computed by parallel_matrix_multiplication function.
        if len(row_data) > 1000 and cpu_count > 1:
            score_matrix = parallel_mcl(score_matrix)
        else:
            score_matrix = mcl(score_matrix)     ## !!! Error 
        clustering(row_data, col_data, score_matrix)


def parallel_mcl(score_matrix):
    """Input score_matrix is the user selected matrix llike BLOASUM45, BLOSUM62 or BLOSUM80"""
    count = 0
    infinitesimal_value = 10**-10
    idempotent_matrix = numpy.matlib.ones((2, 2))

    while idempotent_matrix.sum() > infinitesimal_value:  # > infinitesimal_value
        mcl_time_start = time.time()
        pool = multiprocessing.Pool(cpu_count)  # create a expansion_matrix
        multiplication_results = pool.map(parallel_matrix_multiplication,
                                          zip(score_matrix, repeat(score_matrix)))
        pool.close()
        pool.join()

        # create a inflation_matrix(part 1)
        pool = multiprocessing.Pool(cpu_count)
        power_results = pool.map(
            parallel_matrix_power, multiplication_results)
        pool.close()
        pool.join()

        sum_matrix = 0
        for i in power_results:
            sum_matrix = i + sum_matrix

        # create a inflation_matrix(part 2)
        pool = multiprocessing.Pool(cpu_count)
        divide_results = pool.map(parallel_matrix_division,
                                  zip(power_results, repeat(sum_matrix)))
        pool.close()
        pool.join()

        # Make a Combined matrix for results of parallel_matrix_multiplication function.
        for i in range(len(divide_results)):
            if i == 0:
                score_matrix = divide_results[i]
            else:
                score_matrix = np.concatenate(
                    (score_matrix, divide_results[i]), axis=0)

        sum_results = 0
        for i in multiplication_results:
            sum_results += i
        # identify whether inflation_matrix is idempotent matrix or not.
        idempotent_matrix = abs(np.sum(score_matrix) - sum_results)

        count += 1
        if count > infinite_loop:  # It will prevent the infinite loop of mcl algorithm.
            break
        mcl_time_finish = time.time()
        if verbose:
            print(" mcl time : %f, count : %d, matrix size : %d * %d"
                  % ((mcl_time_finish - mcl_time_start)/60, count, score_matrix[0].size, score_matrix[0].size))
    return score_matrix


def mcl(score_matrix):
    """Input score_matrix is the user selected matrix llike BLOASUM45, BLOSUM62 or BLOSUM80"""
    count = 0
    infinitesimal_value = 10**-10
    idempotent_matrix = numpy.matlib.ones((2, 2))
    #"idempotent_matrix = np.ones((2,2)) # i will later test with this"
    while idempotent_matrix.sum() > infinitesimal_value:  # > infinitesimal_value
        mcl_time_start = time.time()
        expansion_matrix = score_matrix ** 2
        
        print("shape of input score_matrix",score_matrix.shape)
        print("shape of expansion_matrix is ",expansion_matrix.shape )
        print("Shape of inflation matrix is ",expansion_matrix.shape)

        score_matrix = np.power(expansion_matrix, inflation_factor)  ## !!! Eroor
        score_matrix_sum = score_matrix.sum(axis=0)
        # create a inflation_matrix
        score_matrix = np.divide(score_matrix, score_matrix_sum)
        # identify whether inflation_matrix is idempotent matrix or not.
        idempotent_matrix = abs(score_matrix - expansion_matrix)
        count += 1
        if count > infinite_loop:  # It will prevent the infinite loop of mcl algorithm.
            break
        mcl_time_finish = time.time()
        if verbose:
            print(" mcl time : %f, count : %d, matrix size : %d * %d"
                  % ((mcl_time_finish - mcl_time_start)/60, count, score_matrix[0].size, score_matrix[0].size))
    return score_matrix

def clustering(row_data, col_data, score_matrix):
    global cluster_count
    global ortholog_count
    ortholog_temp_list = []
    for i in range(len(row_data)):
        ortholog_list = []
        ortholog = queue.Queue()  # It is Queue which is put the ortholog.
        # It is Queue which is put the ortholog having changed gene ID
        gene_id_queue = queue.Queue()
        for j in range(len(col_data)):
            if 0.1 <= score_matrix[i, j]:
                ortholog.put(col_data[j])
                ortholog_list.append(col_data[j])

        # If the ortholog queue  has the element of ortholog more than 2, it will be printed.
        if ortholog.qsize() >= 3:
            for element in ortholog_list:
                try:
                    ortholog_sum += ortholog_temp_list.index(element)+1

                except ValueError:
                    with open(Cluster_out+"_geneID_S"+str(threshold_score) +
                              "_"+str(inflation_factor), "a") as ortholog_list_save:
                        ortholog_print = "cluster "+str(cluster_count)+" :"
                        ortholog_list_save.write(ortholog_print)

                        while not ortholog.empty():
                            get_ortholog = ortholog.get()
                            ortholog_list_save.write("\t"+get_ortholog)

                            try:
                                # get_ortholog --> ECO_170082288,  get_ortholog.split('_') --> ['ECO', '170082288']
                                get_ortholog_split = get_ortholog.split('_')
                                gene_id_queue.put(
                                    get_ortholog_split[0]+"_"+gene_id_dic[get_ortholog_split[1]])
                            except KeyError:
                                # If the gene_id_dic don't have get_ortholog, it will print the original ID(get_ortholog).
                                gene_id_queue.put(get_ortholog)
                        ortholog_list_save.write("\n")
                    with open(Cluster_out+"_KO_ID_S"
                              + str(threshold_score)
                              + "_"+str(inflation_factor), "a")\
                            as ortholog_list_geneID_save:
                        ortholog_list_geneid_print = "cluster " + \
                            str(cluster_count)+" :"
                        ortholog_list_geneID_save.write(
                            ortholog_list_geneid_print)
                        while not gene_id_queue.empty():
                            ortholog_list_geneID_save.write(
                                "\t"+gene_id_queue.get())
                            ortholog_count += 1
                        cluster_count += 1
                        ortholog_list_geneID_save.write("\n")
                    break
            ortholog_temp_list = operator.concat(
                ortholog_temp_list, ortholog_list)


def parallel_matrix_multiplication(data):
    "Doing parallel_matrix_multiplication Using Numpy"
    matrix_element, matrix = data
    return  matrix_element * matrix


def parallel_matrix_power(matrix_element):
    "This Function Compute Parallel_Matrix Power by using Numpy np.power() Function"
    global inflation_factor
    return np.power(matrix_element, inflation_factor)


def parallel_matrix_division(data):
    "This Function Will Perform Matrix Division by Using Numpy library"
    matrix_element, sum_data = data
    return np.divide(matrix_element, sum_data)


def read_species(pr=1):  # default 1 which will always shows the name of Species
    """ If pr is 1, it will print "Species_List"  Other wise only return the Value in dic Format """
    read_species = os.listdir(command_options.Species)
    selected_species_dic = {}  # list
    backward_selected_species_dic = {}  # list
    for i, species in enumerate(sorted(read_species), start=1):
        selected_species_dic[i] = species
        backward_selected_species_dic[species] = i
        if pr == 1:
            print(str(i)+".", species)
        number = i
    return selected_species_dic, backward_selected_species_dic, number


def del_file(path, file):
    "This Function delete the file passed with path and file"
 #   "Delete the Unnecessary file according to path and file name passed"
 #  del_file =subprocess.Popen(["rm "+path+file], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
 #   del_file_stream = del_file.communicate()
 #   if not del_file_stream[1]:
 #       print ("Done to del "+path+file)
 #   elif del_file_stream[1]:
 #       print (del_file_stream[1])
    try:
        os.remove(path+file)
        print("File Successfully Removed")
        "to remove all files inside the folder we need to use loop over the folder"
    except:
        print("Check the File or Path")


def check_file(file):
    """When mode 3 is Selected.This function will run. file is user input value for clustering output.
    Check the file weather exists or not .
    if exist this will warn to use other name
    Cluster_out is file.Cluster_out(results) + _geneID_S+ str(threshold_score)+"_"\
        +str(inflation_factor) or Cluster_out+"_KO_ID_S"+str(threshold_score)+"_"+str(inflation_factor))
        are files created by mcl algorithm so if these file exist system will exit
    """
    # Declare all variable as a globally Added
    global Cluster_out, threshold_score, infinite_loop
    file_list = glob.glob(file+'*')  # list all related files in same directory local variable
    if (Cluster_out+"_geneID_S"+str(threshold_score)
        + "_"+str(inflation_factor) or Cluster_out+"_KO_ID_S"
            + str(threshold_score)+"_"+str(inflation_factor)) in file_list:
        print("Please, set other name of output.line Number 966")
        sys.exit(2)
    else:
        print(file,"is ")


def read_equal_bbh(path):
    """Read_Equal BBH by user path Blast Best Hit
    the score fiie location is ./score_file/ + species + species like ./score_file/A_C
    """
    global threshold_score

    with open(path+"_oneway_threshold_best_hit_Score"
              + str(threshold_score), 'r') as equal_RBH:
        for j in equal_RBH:
            split_data = j.split()
            # Split Data sample
            # ['A_gi|15605613|ref|NP_212986.1|', 'C_gi|15004707|ref|NP_149167.1|', '51']
            split_data[2] = int(split_data[2]) #Only Converting to int format
            equal_BBH_data.append(split_data)
            equal_BBH_data_dic[split_data[0]] = split_data[1:]
    try:
        with open(path+"_second_oneway_threshold_best_hit_Score"
                  + str(threshold_score), 'r') as second_equal_RBH:
            for j in second_equal_RBH:
                split_data = j.split()
                #"['AAE_gi|15606128|ref|NP_213505.1|', 'AAE_gi|15606877|ref|NP_214257.1|', 898] split_data sample"
                split_data[2] = int(split_data[2])
                second_equal_BBH_data.append(split_data)
    except:
        pass


def read_unequal_bbh(path):
    """ Read unequal BBH path passed by User.
    Since this Function will append value to unequal_BBH_data which created before running this function
    Later i will changed to return type of list from this function.create blank list and append and return"""

    with open(path+"_oneway_threshold_best_hit_Score"
              + str(threshold_score), 'r') as unequal_RBH:
        for j in unequal_RBH:
            split_data = j.split()
            split_data[2] = int(split_data[2])
            unequal_BBH_data.append(split_data)

print("="*82)
print("||\t****Default Variables(Values)****\t\t\t\t\t||")
print("||\tBlastp                 = %s\t\t\t\t\t\t||" % command_options.Blastp)
print("||\tBlastp_data            = %s\t\t\t\t\t||" %
      command_options.Blastp_data)
print("||\tblstp_matrix           = %s\t\t\t\t\t||" %
      command_options.blastp_matrix)
print("||\tcpu_count              = %s\t\t\t\t\t\t||" %
      command_options.cpu_count)
print("||\tgenomes                = %s\t\t\t\t\t\t||" %
      command_options.genomes)
print("||\tinfinite_loop          = %s\t\t\t\t\t\t||" %
      command_options.infinite_loop)
print("||\tinflation_factor       = %s\t\t\t\t\t\t||" %
      command_options.inflation_factor)
print("||\tmode                   = %s\t\t\t\t\t\t||" % command_options.mode)
print("||\tCluster_out            = %s\t\t\t\t\t||" %
      command_options.Cluster_out)
print("||\tthreshold_score        = %s\t\t\t\t\t\t||" %
      command_options.threshold_score)
print("||\tsave_raw_blastp_score  = %s\t\t\t\t\t\t||" %
      command_options.save_raw_blastp_score)
print("||\tScore_file             = %s\t\t\t\t\t||" %
      command_options.Score_file)
print("||\tSpecies                = %s\t\t\t\t\t||" % command_options.Species)
print("||\tverbose                = %s\t\t\t\t\t\t||" %
      command_options.verbose)
print("="*82)
print("\n")
print("Ortholog Detection Program Starts Now\n")

if not sys.argv[1:]:
    """ If not Parameter Passed the Manual Process will Start,
    in case of mode the input parameter can be different so we will not apply the below code function"""
    print("1. BLASTP. \n2. BLASTP using precalculated data. \n3. clustering.\n")
    mode = input(">> Select a mode or modes (1 2 *OR* 1 3 *OR* 2 3): ")
    selected_species_dic, backward_selected_species_dic, number_i = read_species(
        1)
    selected_number = input(
        ">> Select Genomes to detect Orthologs(e.g. 1 2 3 4 5 or 1-5) : ")

    if selected_number.find('-') > 0:
        # find() return the index position of first occurance
        SN = selected_number.split("-")
        if int(SN[-1]) > number_i:  # number_i is length of Genome file inside folder exit the process
            print("\nWrongInput\nInput must be less than", number_i)
            sys.exit(2)
        else:
            user_selected_number = range(int(SN[0]), int(SN[-1])+1)
            for j in user_selected_number:
                print(selected_species_dic[j], end=" ")  # loop in Dic
            print("Selected!!")

    else:
        user_selected_number = sorted(
            set([int(read_species) for read_species in selected_number.split()]))
        # Create a set (remove repeating)
        if int(user_selected_number[-1]) > number_i:
            print("\nWrongInput\nInput must be less than", number_i)
            sys.exit(2)
            # Greater than Genome list will system error
        else:
            for j in user_selected_number:
                print(selected_species_dic[j], end=" ")
            print("Selected!!")
    blastp_matrix = matrix_name()
    cpu_count = int(input("You can use %s Core.\nIf you input >= 2,\
        The Program will run a parallel computation for the blastp.\n" % multiprocessing.cpu_count()
                          + "Enter the number of Core to use in this program (1 ~ %s): " % multiprocessing.cpu_count()))

elif sys.argv[1:]:
    genomes = command_options.genomes
    mode = command_options.mode
    cpu_count = command_options.cpu_count
    blastp_matrix = command_options.blastp_matrix
    inflation_factor = command_options.inflation_factor
    selected_species_dic, backward_selected_species_dic, number_i = read_species()
    user_selected_number = [backward_selected_species_dic[ele]
                            for ele in genomes]   #select by value of dictionary 
    Cluster_out = command_options.Cluster_out     # File Name to save 

Species = command_options.Species
Blastp = command_options.Blastp
Score_file = command_options.Score_file
Blastp_data = command_options.Blastp_data
save_raw_blastp_score = command_options.save_raw_blastp_score
threshold_score = command_options.threshold_score
verbose = command_options.verbose
infinite_loop = command_options.infinite_loop
log_time =datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")

#print("This is for checking ")
#print(command_options.Species ,command_options.Blastp ,command_options.Score_file)
#  This is For test
if "3" in mode:
    "Mode 3 is for clustering the matrix data."
    print("3 mode clustering is selected")
    inflation_factor = input("Enter the inflation factor to cluster: ") #
    Cluster_out = input("Set the name of clustering output Folder : ")
    """Cluster_out is user input value --> the name of clustering output <--results
    """
    check_file(Cluster_out) #if same file for user input species and cluester out and inflation same no need to run cluster again

#del_file(Score_file, "*") #this Function will delete all files inside currently folder

if "3" in mode:
    Log_file_name = "./cluster_out/"+Cluster_out+"_S" + str(threshold_score)+"_"+str(inflation_factor)+log_time+"_log.txt"

elif not "3"  in mode:
    Log_file_name ='log_files/Log' +log_time+"_log.txt"

with open(Log_file_name, 'w') as log:
        log.write(str(datetime.datetime.now()))
        log.write("\nmode :")
        for i in mode :            
            log.write(" "+i)
        log.write("\ngenomes : ")
        for i in user_selected_number :
            log.write(selected_species_dic[i]+" ")    
        log.write("\ncpu_count : "+str(cpu_count))
        log.write("\nblastp matrix : "+blastp_matrix)
        if "3" in mode :
            log.write("\ninflation_factor : "+str(inflation_factor))
            log.write("\nCluster out : "+Cluster_out)
        log.write("\nSpecies : "+Species)
        log.write("\nBlastp : "+Blastp)
        log.write("\nScore file : "+Score_file)  
        log.write("\nBlastp_data : "+Blastp_data)
        log.write("\nsave rawblastp score : "+str(save_raw_blastp_score))
        log.write("\n")   

start_time_OBH = time.time()
if "1" in mode:    # 1 is for Blastp Run
    backward_best_hit_work_list = oneway_threshold_best_hit(mode)
    #print("Back Ward Best Hit result", backward_best_hit_work_list)
    pool = multiprocessing.Pool(cpu_count)
    results = pool.map(backward_best_hit, backward_best_hit_work_list)
    #"function_name is backward_best_hit and the parameter is backward_best_hit_work_list"
    pool.close()
    pool.join()

elif "2" in mode:
    used_precalculated_data_list = []
    new_calculated_data_list = []
    #precalculated_data_list = os.listdir(
    print("Threshold_Score passed is",threshold_score)
    precalculated_data_list = glob.glob(
        Blastp_data+"*oneway_threshold_best_hit_Score"+str(threshold_score)) 
        #_S changed to Score file name and copy score file from score_file to blastp_file folder manually to speed up the system
        # glob.glob function may not work properly in windows so  need to check the  result
    backward_best_hit_work_list = oneway_threshold_best_hit(mode)
    #" This function will also append"
    # If backward_best_hit_work_list is an empty list, pool instance can't finsh the work.
    if not backward_best_hit_work_list == []: 
        # If backward_best_hit_work_list is an empty list, pool instance can't finsh the work.
        "The pool only run if backward best_hit_work_list has some value."
        pool = multiprocessing.Pool(cpu_count)
        # multiprocessing.pool() for Parallel
        results = pool.map(backward_best_hit, backward_best_hit_work_list)
        #backward_best_hit is the function and backward_best_hit_work_list is list parameter
        pool.close()
        pool.join()
    else:
        results = [0, 0]

#del_file("./", "query*") #Delete all files start with query
"Since the  query files will replaced by new files"

finish_time_OBH = time.time()
blastp_time_log = float(((finish_time_OBH - start_time_OBH)/60))
print("BLASTP searches + forward best Hit + backwardbest hit took %f minutes" %
      blastp_time_log)
#Del_file("./","query*") #delete all files of query
with open(Log_file_name, 'a') as log:
    log.write("backward_best_hit took "+str(max(results))+"minutes\n")
    #results will not be empty because of pool not show any results then results = [0,0]
    log.write("BLASTP + Best_Hit + backward_best_hit searches took " +
              str(blastp_time_log)+" minutes\n")

if "3" in mode:
    start_time_clustering = time.time()
    ##########################################################################################
    # generate matrix and calculate the matrix using mcl algorithm and cluster the ortholog."""
    """queue is a linear data structure that stores items in First In First Out(FIFO) manner.
    Create a queue object with a given maximum size.If maxsize is <=0, the queue size is infinite.
    the default
    queue.Queue(maxsize = 0) 
    """
    print("\n>>>> Start mcl algorithm and clustering ortholog <<<<")
    "the values of these blank list and dic is added inside the function which is very bad"
    equal_BBH_data = []
    unequal_BBH_data = []
    equal_BBH_data_dic = {}
    second_equal_BBH_data = []
    results = queue.Queue()
    tasks = queue.Queue()
    cluster_count = 1
    ortholog_count = 0
    gene_id_dic = {}

    with open("./db/myva=gb", "r") as id_read:
        # myvba=gb is a database
        for i in id_read:
            gene_name, gene_id = i.split()
            gene_id_dic[gene_id.replace("\n", "")] = gene_name  # remove "\n"

    if "1" in mode:
        for i in user_selected_number:
            for k in user_selected_number:
                if k < i:
                    #print(k,"is less than ", i)
                    pass
                elif i == k:
                    "If the both species is same, then read_equal_bbh() function will run"
                    read_equal_bbh(
                        Score_file+selected_species_dic[i]+"_"+selected_species_dic[k])
                elif i != k:
                    "If both the species are different then read_unequal_bbh() function will run"
                    read_unequal_bbh(
                        Score_file+selected_species_dic[i]+"_"+selected_species_dic[k])

    elif "2" in mode:
        for used_data in used_precalculated_data_list:
            "used_precalculated_data_list is return from function one_way_threshold_best_hit(mode)  "
            first, second = used_data.split("_")
            if first == second:
                #if Both genes are same
                read_equal_bbh(Blastp_data+used_data)
            elif first != second:
                read_unequal_bbh(Blastp_data+used_data)

        for new_data in new_calculated_data_list:
            #"new_calculated_data_list is lilst appended in Function Oneway_threshold_best_hit(mode =2) "
            first, second = new_data.split("_")
            if first == second:
                read_equal_bbh(Score_file+new_data)
            elif first != second:
                read_unequal_bbh(Score_file+new_data)

    matched_BBH_data = []
    matched_BBH_element_data_set = []

    for unequal_RBH_element in unequal_BBH_data:
        matching_bbh(unequal_RBH_element)
        temp_results_list = []
        if results._qsize() != 0:  # return the number of results as Queue.
            while not results.empty():
                get_results = results.get()
                temp_results_list.append(get_results)
            matched_BBH_data.append(temp_results_list)

    for data in matched_BBH_data:
        matrix_clustering_ortholog(data)
        #matrix_clustering_ortholog(data, bar)
    finish_time_clustering = time.time()
    mcl_time_log = float((finish_time_clustering - start_time_clustering)/60)
    remark_time_log = float((finish_time_clustering - start_time_OBH)/60)
    print("mcl algorithm and Ortholog clustering took %.2f minutes" % mcl_time_log)
    print("owPRemark program took %.2f minutes" % remark_time_log)

    if "3" in mode:
        with open(Log_file_name, 'a') as log:
            log.write("Ortholog count : "+str(ortholog_count)+"," +
                      " Cluster count : "+str(cluster_count-1)+"\n")
            log.write("mcl algorithm and Ortholog clustering took " +
                      str(mcl_time_log)+" minutes\n")
            log.write("XXX program took "+str(remark_time_log) + "minutes\n")