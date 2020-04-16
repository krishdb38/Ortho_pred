#!/usr/bin/python

"""
This_is_a_program_to_detect_ortholog_between_species.

"""
import sys
import os                  # For system operation
import glob                # list the file unix system
import subprocess          # to run blastp from OS
import time                # to display time
import multiprocessing  # For Parallel Processing
import queue            # For Multiprocessing
import copy             # To copy file
import operator         # to concate the list we will remove this later
import argparse
import datetime         # to display current time & calculate difference
from itertools import repeat
import numpy as np
import numpy.matlib      # For Mathmatical (Algebra) Operation
from alive_progress import alive_bar  # For Progress Bar
from Bio import SeqIO

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


PARSER = argparse.ArgumentParser(description="""
This program help you to detect ortholog between the protein sequences
from different genomes and cluster the orthologs to ortholog groups.""",
                                 add_help=True, prefix_chars='-+')
ESSENTIAL_GROUP = PARSER.add_argument_group('Essential arguments')
BLASTP_GROUP = PARSER.add_argument_group('blastp optional arguments')
MCL_GROUP = PARSER.add_argument_group('clustering optional arguments')
ESSENTIAL_GROUP.add_argument('-M', '--MODE',
                             action='store',
                             nargs='+',
                             dest='MODE',
                             choices=['1', '2', '3'],
                             help="""
                                1 : BLASTP
                                2 : BLASTP using precalculated data.
                                3 : clustering. ex) 1 or 2 or 1 3 or 2 3
                                """)
ESSENTIAL_GROUP.add_argument('-g', '--genomes',
                             nargs='+',
                             dest='genomes',
                             help="Select a genomes to detect orthologs.")
PARSER.add_argument('-c', '--cpu',
                    action='store',
                    default='1',
                    dest='cpu_count',
                    type=int,
                    help='Set the number of CPU to use in program.\
                        The default is "1".\
                        This system has '
                    + str(multiprocessing.cpu_count()) + " CPUs")
BLASTP_GROUP.add_argument('+r',
                          action='store_true',
                          default=True,  # Default was false
                          dest='save_raw_blastp_score',
                          help='Save the raw score of blastp to\
                            the score_file dirctory.')
BLASTP_GROUP.add_argument('-t', '--threshold',
                          action='store',
                          default=5,
                          type=int,
                          dest='threshold_score',
                          help="""Set the threshold score.
                          Threshold score is an allowable range for ortholog
                          If the score difference between backward best hit score and
                          forward best hit score is less than threshold score,
                          it is considered the forward best hit pair to be the ortholog.'
                          'That is  called the "one-way threshold best hit" by us.
                          The default value is "5".
                          "backward best hit score - forward best hit score <= threshold score""")
BLASTP_GROUP.add_argument('-m', '--matrix',
                          action='store',
                          default='BLOSUM62',
                          dest='blastp_matrix',\
                          # 82 is not availiable
                          choices=['BLOSUM45', 'BLOSUM62', 'BLOSUM80'],\
                          help='select the matrix to use in blastp.\
                              The default value is "BLOSUM62".')
BLASTP_GROUP.add_argument('-s', '--species',
                          action='store',
                          default='./species/',
                          dest='SPECIES',
                          help='set the path of species directory.\
                              The default is  "./species/".')
BLASTP_GROUP.add_argument('-b', '--blastp',
                          action='store',
                          default="blastp",  # blastp must be in your system path
                          dest='BLASTP',
                          help='set the path of blastp file to run the blastp program.\
                              The default is "blastp".')
BLASTP_GROUP.add_argument('-F', '--scorefile',
                          action='store',
                          default='./score_file/',
                          dest='Score_file',
                          help='The default path is "./score_file/".')
BLASTP_GROUP.add_argument('-B', '--blastp_data',
                          action='store',
                          default='./blastp_data/',
                          dest='BLASTP_DATA',
                          help='set the path of precalculated blastp data.\
                              The default path is "./blastp_data/".')
MCL_GROUP.add_argument('-i', '--IF',
                       action='store',
                       default=1.4,
                       dest='inflation_factor',
                       type=float,
                       help='The default value is "1.4".')
MCL_GROUP.add_argument('-l', '--loop',
                       action='store',
                       default=60,
                       type=int,
                       dest='infinite_loop',
                       help='prevent the infinite loop of mcl algorithm.\
                                 The default value is "60".')
MCL_GROUP.add_argument('-o', '--out',
                       action='store',
                       default='./cluster_out',
                       dest='CLUSTER_OUT',
                       help=""""
                       set the path and name of ortholog cluster file,
                       log and error_warning file.
                       The default path and name is "./cluster_out""")
MCL_GROUP.add_argument('+v',
                       action='store_true',
                       default=False,
                       dest='verbose',
                       help='verbosely show information of a big matrix computaion.')
PARSER.add_argument('--version',
                    action='version',
                    version='%(prog)s Ver. 11.1')
COMMAND_OPTIONS = PARSER.parse_args()


def matrix_name():
    """This Function will Return the Matrix name choosed by User.
    BLOSUM45 ,  BLOSUM62 , BLOSUM80 """
    print("BLOcks SUbstitution Matrix (BLOSUM) is a Substitution matrix\
        used for sequence alignment of Proteins")
    print("""
    1. BLOSUM45 :-For more distantly related Proteins alignment DataBase
    2. BLOSUM62 :- MidRange Seq with more than 62%similarity
    3. BLOSUM80 :- More related Proteins\n
    Other Keys will exit the Program
    """)
    metrix_num = int(input("\nEnter a matrix number:\t"))
    matrix_dic = {1:"BLOSUM45", 2:"BLOSUM62", 3:"BLOSUM80"}
    if metrix_num  in (1, 2, 3):
        print(matrix_dic[metrix_num], "Matrix selected")
        return matrix_dic[metrix_num]
    else:
        print("Wrong input *%s*Sorry not in list\n" %
              metrix_num, "*"*20, "Good Bye", "*"*20, "\n")
        sys.exit(2)

def check_mode(mode_):
    check = any(x in mode_ for x in ["1","2"])
    if not check:
        print("You Must Enter either 1 or 2")
        sys.exit(2)


def query_sequence(genome):
    """
    This Function Read Fastaq Files and return as a list Format with gene Position as a index
    This Function create a list of Seprate Query Sequence
    # !! This Function can also be calculated busing Bio Python but some error rise !!
    """
    gene_seq = ""
    gene_seq_list = []  # container to  store gene id with seq
    try:
        with open(SPECIES+genome) as gene:
            # open genome file like AAE HIN
            for each_line in gene:
                if ">" in each_line and gene_seq != "":
                    # starting with ">" and not overlap
                    gene_seq_list.append(gene_seq)
                    # added gene name to the first and gene_seq will blank
                    gene_seq = ""
                gene_seq = gene_seq+each_line
                # add each line of seq till ">" not found
            gene_seq_list.append(gene_seq)  # last position gene will added
            return gene_seq_list
    except IOError as err:
        print("IOError occurred in query_sequence function : " + str(err))


def query_seq_(genome):
    """This Function Read Fasta (Genome) file and return as a list Format
    with gene Position and Sequence Developed by Krish.
    The result of this function is not running while blast """
    genome = SPECIES+genome
    try:
        return [(seq_record.id+"\n"+str(seq_record.seq)) for seq_record in SeqIO.parse(genome, "fasta")]
        # To understand this Function Bio Python library needs to be studied
    except IOError as err:
        print(str(err))


def write_query(query, parallel_num):
    "This Function Write Query with file Name query + parallel_num in same directory\
        and raise IO error (Folder?Files not Found) if Error occured"
    # print("write_query(query,parallel_num)")
    try:
        with open("./query/query_"+str(parallel_num), "w") as write_:
            write_.write(query)
    except IOError as err:
        print(
            "IOError occurred in write_query function Folder or File Not found: " + str(err))


def run_blast(subject, parallel_num):
    """
    By this Function it will create a Pipe line to run blastp in Computer by input Parameter.
    Subject is whole Genome and parallel_number is query file Created in early step.
    Running a big file is time Consuming.The parameters passed in blastp are
    The Output Format 10 --> Comma Separated Values
        >> qseqid --> Query Seq ID
        >> ssequid --> Subject of Seq. id ID
        >> Score -->Raw Score
        >> length-->Alignment length
    # query subject is inside query Folder
    # run_blastp =subprocess.Popen(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # run_blastp = subprocess.run(\
    # cmd,shell=True ,capture_output = True,text = True)
    # output Format 10 qseqid query (e.g. gene sequence id,\
    #  sseqid subject (e.g. reference genome) genome id)"""

    subject = SPECIES+subject  # ** later remove this
    cmd = ["blastp", "-query", "./query/query_"+str(parallel_num), "-subject", subject,
           "-matrix", blastp_matrix, "-outfmt", "10 qseqid sseqid score length"]

    run_blastp = subprocess.run(cmd, capture_output = True, text = True)
    # Only for version > 3.7
    if run_blastp.stderr:
        # "If Error rise log file will print and write into log file"
        print(run_blastp.stderr)
        with open(Log_file_name, 'a+') as err:
            err.write(str(run_blastp.stderr))
    ## !!_____if Error rise better exit loop here____
    return run_blastp.stdout


def same_species_forward_best_hit(blastp_score):
    """Search the forward best hit among the blastp scores of
    same species.Because there are duplicated genes in a same genome.
    blstp_score file is return value of run_blast() Function
    When the blastp score compare with blastp score of duplicate gene,
    if score and length are same,
    blasp score of duplicated gene is added to a second best score."""
    # print(blastp_score)
    bpscore_split_list = []
    temp_best_score = ['-1', '-1', '-1']
    second_temp_best_score = []
    best_score = []
    second_best_score = []
    bpscore_split = blastp_score.split("\n")
    # delete of [''] in the last index  ex) ['gi,gi,1,1','gi,gi,2,2','']
    del bpscore_split[-1]
    for _ in bpscore_split:
        bpscore_element = _.split(',')
        bpscore_split_list.append(bpscore_element)
    # ex) k is ['gi|15605613|ref|NP_212986.1|', 'gi|15605613|ref|NP_212986.1|', '3702', '699']
    for k in bpscore_split_list:
        # print(len(k)) ==4
        # example of k is
        # ["gi|15605616|ref|NP_212989.1","gi|15606180|ref|NP_213700.1","58","30"]
        #
        #time.sleep(5)
        if k[0] == k[1]:
            # same position gene will give high score,best_score
            best_score.append(k)
        elif k[0] != k[1]:
            if int(k[2]) > int(temp_best_score[2]):  # Compare score
                temp_best_score = k
            elif int(k[2]) == int(temp_best_score[2]):  # Is Equal
                if int(k[3]) > int(temp_best_score[3]):  # compare length
                    temp_best_score = k
                elif int(k[3]) == int(temp_best_score[3]):
                    second_temp_best_score.append(k)
    #   print ("############ temp best score ############", temp_best_score)
    second_best_score.append(temp_best_score)
    for j in second_temp_best_score:
        if j[2] == temp_best_score[2] and j[3] == temp_best_score[3]:
            second_best_score.append(j)
    for bscore in second_best_score:
        # m is replaced with bscore (best_score)
        # !! list index out of range ERROR
        # print("best_score",best_score)
        # print("bscore",bscore)
        # time.sleep(5)
        try:
            if (best_score[0][2] == bscore[2] and int(best_score[0][3]) <= int(bscore[3]))\
                or int(best_score[0][2]) < int(bscore[2]):
                # '104' < '23' is True because of string.
                # So the int function is used.
                best_score.append(bscore)
        except :
            pass
    #  print("best Score is", best_score)
    return best_score


def forward_best_hit(blastp_score):
    """ Search the forward best hit among the blastp scores of same species."""
    #  print("Rnunning GetForward BestHit")
    bpscore_split_list = []
    temp_best_score = ['-1', '-1', '-1']
    second_temp_best_score = []
    best_score = []
    bpscore_split = blastp_score.split("\n")
    # delete of ['']   ex) ['gi,gi,1,1','gi,gi,2,2','']
    del bpscore_split[-1]
    for _ in bpscore_split:
        bpscore_element = _.split(',')
        bpscore_split_list.append(bpscore_element)
    for k in bpscore_split_list:
        # ex) k is ['gi|15605613|ref|NP_212986.1|', 'gi|15605613|ref|NP_212986.1|', '3702', '699']
        # print ">>>>>>>>>>>>>>>forward_best_hit   k", k
        if int(k[2]) > int(temp_best_score[2]):  # Compare the Score of the Blast
            temp_best_score = k
        elif int(k[2]) == int(temp_best_score[2]):
            if int(k[3]) > int(temp_best_score[3]):  # Compare the Lenth of Sequence
                temp_best_score = k
            elif int(k[3]) == int(temp_best_score[3]):
                second_temp_best_score.append(k)
    best_score.append(temp_best_score)
    for j in second_temp_best_score:
        if (j[2] == temp_best_score[2]) and (j[3] == temp_best_score[3]):
            best_score.append(j)
    return best_score, bpscore_split_list


def division_parallel_query(query_v, query_division_value, cpu_count_, query_v_len):
    """query_v is list of queries, query_division_value is ,
        cpu_count_ is the Number of CPU User Entered and
        query_v_len is the length of the query"""
    parallel_query = []
    parallel_query_start = 0
    if query_v_len % cpu_count_ == 0:  # perfect division
        for i in range(cpu_count_):
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
        for i in range(cpu_count_):
            i += 1
            if parallel_query_start == 0:
                parallel_query.append(
                    query_v[int(parallel_query_start):i*int(query_division_value)])
                parallel_query_start += 1
            elif i < cpu_count_:
                parallel_query.append(query_v[int(
                    parallel_query_start)*int(query_division_value):i*int(query_division_value)])
                parallel_query_start += 1
            elif i == cpu_count_:
                parallel_query.append(
                    query_v[int(parallel_query_start)*int(query_division_value):])
                parallel_query_start += 1

    return parallel_query


def run_parallel_query(species_of_query, species_of_subject, query_v, parallel_num):
    """
    run_parallel_query(i ,k , query_v , cpu_count)
    i and k are user selected number in a list Format
    Run the following functions. write_query, run_blast,
    same_species_forward_best_hit, forward_best_hit
    Save the files which are oneway_ threshold_best_hit,
    second_oneway _threshold_best_hit,
    bpscore_split_list and raw_blastp_score(optional) by each species.
    parallel_num is the Number of CPU selected by the user if CPu 1 selected then 1
    # ! query_v = query_sequence(SELECTED_SPECIS[i])
    """
    # ! query_v ! = query_sequence(SELECTED_SPECIS[species_of_query]) why ??
    with alive_bar(len(query_v)) as bar_:
        # declare your set of items for loop, bar is blacklisted
        for j in query_v:
            # query_v is the list of the genes from genome
            write_query(j, parallel_num)
            # if 1 Added Here Add also to Run
        # This Function Only Write a file with j name and parallel_num i.e CPU Count
            blastp_score = run_blast(
                SELECTED_SPECIS[species_of_subject], parallel_num)
        # Return Byte File type, and this is only for one gene position query and Whole Genome
            if blastp_score != '':  # if blastp run successfully
                best_score, bpscore_split_list = forward_best_hit(
                    blastp_score)  # .decode() Convert in to str format
                # ex) AAE == AAE. It will save best_score without reversing run_blast.
                # The sample for blastp_score and bpscore_split_list is saved inside folder

                if species_of_query == species_of_subject:
                    same_species_fbs = same_species_forward_best_hit(
                        blastp_score)  # .decode()
                    for best_score_element in same_species_fbs:
                        # ex) [A1 of AAE, A1 of AAE, 30]
                        if best_score_element[0] == best_score_element[1]:
                            # if same species add write file in # !!!
                            with open(Score_file+SELECTED_SPECIS[species_of_query]\
                                + "_" +SELECTED_SPECIS[species_of_subject]\
                                + "_oneway_tbh_Score"\
                                + str(threshold_score), "a") as oneway_threshold_bh:
                                # best_hit is replaced with h
                                save_best_score = SELECTED_SPECIS[species_of_query]\
                                    + "_" + best_score_element[0]\
                                    .split(r"\s")[0]\
                                    + " "\
                                    + SELECTED_SPECIS[species_of_subject]\
                                    + "_"+best_score_element[1]\
                                    .split(r"\s")[0]+" " + best_score_element[2]+"\n"
                                # best_score_element[0].split("|")
                                # ==> ['gi','15642790', 'ref', 'NP_227831.1', '']
                                oneway_threshold_bh.write(
                                    save_best_score)
                        else:  # ex) [A1 of AAE, A2 of AAE, 30]
                            with open(Score_file+SELECTED_SPECIS[species_of_query]\
                                      +"_"\
                                      + SELECTED_SPECIS[species_of_subject]\
                                      + "_second_oneway_tbh_Score"\
                                      + str(threshold_score), "a") as second_oneway_threshold_bh:
                                      # best_hit is renamed with bh
                                second_save_best_score = SELECTED_SPECIS[species_of_query]\
                                    + "_" + best_score_element[0].split(r"\s")[0]\
                                    + " " + SELECTED_SPECIS[species_of_subject]\
                                    + "_" + best_score_element[1].split(r"\s")[0]\
                                    + " " + best_score_element[2]+"\n"
                                second_oneway_threshold_bh.write(
                                    second_save_best_score)
                else:
                    # If species_of_query not equal with species_of_subject, run reversing run_blast
                    for best_score_element in best_score:
                        if '-1' not in best_score_element:   # ! if not "-1" in best_score_element:
                            with open(Score_file+SELECTED_SPECIS[species_of_query]\
                                    + "_"\
                                    + SELECTED_SPECIS[species_of_subject]\
                                    + "_"+"best_score_S"\
                                    + str(threshold_score)\
                                    + "_"\
                                    + str(parallel_num), "a") as save_best_hit:
                                best_score_save = SELECTED_SPECIS[species_of_query]\
                                    + "_"\
                                    + best_score_element[0].split(r"\s")[0]\
                                    + " "+SELECTED_SPECIS[species_of_subject]\
                                    + "_"+best_score_element[1].split(
                                        r"\s")[0]\
                                    + " " + best_score_element[2]\
                                    + " " + best_score_element[3]+"\n"
                                save_best_hit.write(best_score_save)

                    for bpscore_split_list_element in bpscore_split_list:
                        with open(Score_file+SELECTED_SPECIS[species_of_query]+"_"
                                  + SELECTED_SPECIS[species_of_subject]
                                  + "_"+"bpscore_split_list_S"
                                  + str(threshold_score)
                                  + "_"+str(parallel_num), "a") as save_bpscore_split_list:

                            bpscore_split_list_save = SELECTED_SPECIS[species_of_query]\
                                    + "_"\
                                    + bpscore_split_list_element[0].split(\
                                        r"\s")[0]\
                                    + " "+SELECTED_SPECIS[species_of_subject]\
                                    + "_"+bpscore_split_list_element[1].split(r"\s")[0]\
                                    + " " + bpscore_split_list_element[2]\
                                    +" " + bpscore_split_list_element[3]+"\n"
                            save_bpscore_split_list.write(
                                bpscore_split_list_save)
            else:
                print("--------")
            bar_()  # Call after Consuming One Item
            if save_raw_blastp_score:
                with open(Score_file+SELECTED_SPECIS[species_of_query]
                          + "_"+SELECTED_SPECIS[species_of_subject]+"_S"
                          + str(threshold_score)+"_"+str(parallel_num), "a") as save_blastp:
                    # Only string is format is supported
                    save_blastp.write(blastp_score)  # .decode()
    # bar.finish() # progressing bar finish


def oneway_tbh(mode_):
    """ This Function accept the mode and Run program and return
     bbh_work_list mode 1 and 2 is supported
     python continue_loop """

    info = "Running the blastp & forward best hit searches "
    process_list = []
    bbh_work_list_ = []
    # !! need to varify global variable
    if "1" in mode_:
        # """We have 3 Mode 1 is for blastp,
        # Mode 2 is for BLASTP using precalcualted data and
        # Mode 3 is for clustering"""
        for i in USER_SELECTED_NUMBER:  # Select species to write query
            query_v = query_sequence(SELECTED_SPECIS[i])
            # query_v is a list Format  with a position gene id , Seq
            # User Selected  [1, 3, 5] is list of User input
            for k in USER_SELECTED_NUMBER:
                # ! changed by Krish
                # ! while can be used but create a infinity loop
                if k < i:
                    continue
                # nothing to do go to tp
                else:
                    print(info+" between %s genome and %s genome"
                          % (SELECTED_SPECIS[i], SELECTED_SPECIS[k]))
                    # this will Create a loop though 2 times
                    # length of Genome gene ID in Sequence File
                    query_v_len = len(query_v)
                    if cpu_count == 1:
                        # "No Parallel Computing while cpu count == 1"
                        blastp_time_start = time.time()
                        run_parallel_query(i, k, query_v, cpu_count)
                        # run_parallel_query doesnot return any value but it write
                        # i is first species k is second species
                        # number queryv is list file of i Position genome"
                        blastp_time_end = time.time()
                        print(info + "took %.2f minutes" % (
                            (blastp_time_end-blastp_time_start)/60))
                    else:
                        # If the number of query_v_len is less than cpu_count,
                        # Remark will select the number of query_v_len.
                        if query_v_len < cpu_count:
                            # the cpu_count will seprate list items from query_v file.
                            # this will not happened because length of query is always long"
                            blastp_time_start = time.time()
                            # 1 is query_division_value.
                            # Because query_v_len / query_v_len(=cpu_count) is 1.
                            parallel_query = division_parallel_query(
                                query_v, 1, query_v_len, query_v_len)
                            for m in range(query_v_len):
                                process = multiprocessing.Process(
                                    target=run_parallel_query, args=(i, k, parallel_query[m], m+1))
                                # args( i => species of query , k => species of subject,
                                # m+1 => cpu_count ex) 1, 2 ...)
                                process_list.append(process)
                                process.start()
                            for n in process_list:
                                n.join()
                            blastp_time_end = time.time()
                            print(info+"took %.2f minutes " % (
                                (blastp_time_end-blastp_time_start)/60))
                        else:
                            # "this wiil run if cpu_count is more than 1"
                            blastp_time_start = time.time()
                            query_division_value = query_v_len / cpu_count
                            parallel_query = division_parallel_query(
                                query_v, query_division_value, cpu_count, query_v_len)
                            for m in range(cpu_count):
                                process = multiprocessing.Process(
                                    target=run_parallel_query, args=(i, k, parallel_query[m], m+1))
                                # args(i=>species of query,k=>species of subject,
                                # m+1 => cpu_count ex) 1, 2 ...)
                                process_list.append(process)
                                process.start()
                            for n in process_list:
                                n.join()
                            blastp_time_end = time.time()
                            print(info+"took %.2f minutes " % (
                                (blastp_time_end-blastp_time_start)/60))
                if not i == k:
                    bbh_work_list_.append((i, k, query_v_len))
                    print("Backward")
    elif "2" in mode_:

        # If blast Files exist then this will pass the blast time other wise blast will run.
        # This will reduce the time for blast and increase the speed of the system.
        print("\nMode 2, Running from precalculated data sets which will save the Time")
        # print(len(precalculated_data_list),"oneway_tbh_ files are located")
        # # !!! Delete This row later
        for i in USER_SELECTED_NUMBER:  # Select species to write query
            query_v = query_sequence(SELECTED_SPECIS[i])
            for k in USER_SELECTED_NUMBER:  # Select of subject
                if BLASTP_DATA+SELECTED_SPECIS[i]+"_"+SELECTED_SPECIS[k]\
                    + "_oneway_tbh_Score"\
                        + str(threshold_score) in precalculated_data_list:
                    print("Skipped")
                    used_precalculated_data_list.append(
                        SELECTED_SPECIS[i]+"_"+SELECTED_SPECIS[k])
                    # used_precalcualted _data_list is created before calling this function
                    # if files exist then passed
                    # continue  # File will skkipped # !! Removed
                else:
                    if k < i:  # gene ====> query 1->1 1->2 1->3 2->2 2->3
                        continue
                    else:
                        # same as Mode_ 1, if Precalculated file doesnot exist
                        # print(info+ " between %s genome and %s genome" % (
                        #    SELECTED_SPECIS[i], SELECTED_SPECIS[k]))
                        query_v_len = len(query_v)
                        if cpu_count == 1:
                            blastp_time_start = time.time()
                            run_parallel_query(i, k, query_v, cpu_count)
                            blastp_time_end = time.time()
                            print(info + "took %.2f minutes" % (
                                (blastp_time_end-blastp_time_start)/60))
                        else:
                            # If the number of query_v_len is less than cpu_count,
                            # Remark will select the number of query_v_len.
                            if query_v_len < cpu_count:
                                blastp_time_start = time.time()
                                # 1 is query_division_value.
                                # Because query_v_len / query_v_len(=cpu_count) is 1.
                                parallel_query = division_parallel_query(
                                    query_v, 1, query_v_len, query_v_len)
                                for m in range(query_v_len):
                                    process = multiprocessing.Process(
                                        target=run_parallel_query,
                                        args=(i, k, parallel_query[m], m+1))
                                    # args( i => species of query , k => species of subject,
                                    # m+1 => cpu_count ex) 1, 2 ...)
                                    process_list.append(process)
                                    process.start()
                                for n in process_list:
                                    n.join()
                                blastp_time_end = time.time()
                                print(info+" took %.2f minutes" % (
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
                                    # args( i => species of query , k => species of subject,
                                    # m+1 => cpu_count ex) 1, 2 ...)
                                    process_list.append(process)
                                    process.start()
                                for n in process_list:
                                    n.join()
                                blastp_time_end = time.time()
                                print(info + " took %.2f minutes" % (
                                    (blastp_time_end-blastp_time_start)/60))
                        new_calculated_data_list.append(
                            SELECTED_SPECIS[i]+"_"+SELECTED_SPECIS[k])
                        if not i == k:
                            bbh_work_list_.append(
                                (i, k, query_v_len))
    return bbh_work_list_


def backward_best_hit(args):
    """ This will  calculate backward_best_hit with
    This function open 2 files from score_files folder
    Genome_genome_best_score_S_5_parallel_num and
    Genome_Genome_bpscore_split_list_S_5_parallel_num
    and convert last 2 digits to integer values.
    and again loopdone in
    forward_best_hit_score_list and
    bpscore_split_list

    """
    bck_info = "Running the backward_best_hit"
    # print(bck_info)
    species_of_query, species_of_subject, query_v_len = args
    start_time_bbh = time.time()
    forward_best_hit_score_list = []
    bpscore_split_list = []
    print(bck_info + " between %s genome %s genome" % (
        SELECTED_SPECIS[species_of_query], SELECTED_SPECIS[species_of_subject]))
    # If the number of query_v_len is less than cpu_count, the cpu_count is changed to query_v_len.
    if query_v_len < cpu_count:
        for parallel_num in range(query_v_len):
            parallel_num += 1
            with open(Score_file+SELECTED_SPECIS[species_of_query]+"_"
                      + SELECTED_SPECIS[species_of_subject]\
                      + "_"+"best_score_S"\
                      + str(threshold_score)\
                      +"_"+str(parallel_num), "r") as best_hit_score:
                    # AAE_gi|15605613|ref|NP_212986.1| ECO_gi|170082858|ref|YP_001732178.1| 3020 703
                    # AAE_gi|15605614|ref|NP_212987.1| ECO_gi|170083440|ref|YP_001732760.1| 1990 404
                    # sample of best_hit_score File
                for each_line in best_hit_score:
                    # the uppper line exmple is passed
                    split_each_line = each_line.split(" ")
                    # = ['AAE_gi|15605613|ref|NP_212986.1|', 'ECO_gi|170082858|ref|YP_001732178.1|', '3020', '703']

                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
                    forward_best_hit_score_list.append(split_each_line)

            # convert the blastp result score length and score to integer format
            with open(Score_file+SELECTED_SPECIS[species_of_query]+"_"
                      + SELECTED_SPECIS[species_of_subject] +
                      "_"+"bpscore_split_list_S"
                      + str(threshold_score)+"_"+str(parallel_num), 'r') as blastp_score:
                    # AAE_gi|15605613|ref|NP_212986.1| CAC_gi|15004754|ref|NP_149214.1| 80 55
                    # AAE_gi|15605613|ref|NP_212986.1| CAC_gi|15004707|ref|NP_149167.1| 69 14
                for each_line in blastp_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
                    bpscore_split_list.append(split_each_line)

    else:
        # if cpu length is smaller than queryV length
        for parallel_num in range(cpu_count):
            # cpu_count is the number of cpu used by the user
            # ! Error may raise if parallel Num is not in the file
            parallel_num += 1
            with open(Score_file+SELECTED_SPECIS[species_of_query]+"_"
                      + SELECTED_SPECIS[species_of_subject]+"_"+"best_score_S"
                      + str(threshold_score)+"_"+str(parallel_num), "r") as best_hit_score:
                for each_line in best_hit_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[3] = int(split_each_line[3])
                    forward_best_hit_score_list.append(split_each_line)
            with open(Score_file+SELECTED_SPECIS[species_of_query] + "_"
                      + SELECTED_SPECIS[species_of_subject]
                      + "_"+"bpscore_split_list_S"
                      + str(threshold_score)+"_"+str(parallel_num), 'r') as blastp_score:
                for each_line in blastp_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
                    bpscore_split_list.append(split_each_line)

    # "bar = Bar("Searching : "+SELECTED_SPECIS[species_of_query]
    # +"-"+SELECTED_SPECIS[species_of_subject], max = len(forward_best_hit_score_list)) ""

    for forward_best_hit_score_element in forward_best_hit_score_list:
        matching_list = []
        backward_best_score = ['-1', '-1', '-1']
        # length is 3 bcz the length is 3 in
        # example of bpscore_split_list
        # [AAE_gi|15605613|ref|NP_212986.1| ECO_gi|170082858|ref|YP_001732178.1| 3020 703,
        # AAE_gi|15605614|ref|NP_212987.1| ECO_gi|170083440|ref|YP_001732760.1| 1990 404 ]
        # botth bpscore_split_list and forward_best_hit_score_list shape are almost same
        # for element in bpscore_split_list:
        #    if element[1] == forward_best_hit_score_element[1]:
        #        matching_list.append(element)
        # ! The upper 3 for loop be converted list Comprehensive
        matching_list = [element for element in bpscore_split_list\
            if element[1]== forward_best_hit_score_element[1]]

        for element in matching_list:
            if int(element[2]) > int(backward_best_score[2]):
                backward_best_score = element
    # !! Future Purpose
    #        with open('./'+SELECTED_SPECIS[species_of_query]\
    # +"_"+SELECTED_SPECIS[species_of_subject]+'_subtraction
    # '+"_"+str(threshold_score), 'a') as subtraction :
    #            save_data = int(backward_best_score[2]) - int(forward_best_hit_score_element[2])
    #            subtraction.write(str(save_data)+"\n")

        if int(backward_best_score[2]) - int(forward_best_hit_score_element[2]) <= threshold_score:
            # the default threshold_score passed is 5
            with open(Score_file+SELECTED_SPECIS[species_of_query]\
                      + "_"+SELECTED_SPECIS[species_of_subject]\
                      + "_oneway_tbh_Score"+str(threshold_score), "a")\
                    as other_oneway_tbh:
                save_data = forward_best_hit_score_element[0]\
                    + " "+forward_best_hit_score_element[1]\
                    + " " + str(int(forward_best_hit_score_element[2]))+"\n"
                other_oneway_tbh.write(save_data)

    # bar.finish()
    finish_time_bbh = time.time()
    rbh_time = float((finish_time_bbh - start_time_bbh)/60)
    print(bck_info + " %s-%s took %.2f minutes" %
          (SELECTED_SPECIS[species_of_query],\
              SELECTED_SPECIS[species_of_subject], rbh_time))
    return rbh_time


def search_equal_bbh_data(target_a):
    """Search the equal backward best hit data. ex) AAE_AAE_backward_best_hit
 """
    # print("search_equal_bbh_data")
    # ! equal_BBH_data_dic is a blank dict defined in begining if mode = 3
    # target_a -->  [HIN_gi |68249031|ref|NP_149167.1,771]
    # [CAC_gi |15004707 |ref |NP_149167.1 , 771]
    # target_a is list_type
    # type of target_a is str
    # put_data = ["AAE_gi|15600000|ref|NP_213710.1",0]
    # target_a = AAE_gi|15600000|ref|NP_213710.1
    # erqual_BBH_data_dic =
    # {"AAE_gi|15600000|ref|NP_213710.1":["AAE_gi|15600000|ref|NP_213710.1",0]}
    # equal_BBH_data_dic is a big dictionary file

    try:
        put_data = equal_BBH_data_dic[target_a]
        # ! if target_a not found in equal_BBH_data will not create error
        # !!! returns value or default
        if put_data[1] != 0:
        # ! code changed if error raise check with old
            copy_put_data = copy.copy(put_data)
            copy_put_data.insert(0, target_a)
            results.put(copy_put_data)
            equal_BBH_data_dic[target_a][1] = 0
            # change the number value to 0
            for i in second_equal_BBH_data:
                # example of i
                # ["AAE_gi|15606581|ref|NP_213505.1","AAE_gi|15606877|ref|NP_214257.1",898]
                # ["AAE_gio|15606744|ref|NP_214124.1|","AAE_gi|15606128|ref|NP_213961.1",2170]
                if i[2] != 0 and (i[0] == target_a or i[1] == target_a):
                    # !! code changed if error raise check with old
                    copy_second_put_data = copy.copy(i)
                    # Don't put results as queue.
                    # Because the tasks will put copy_second_put_data to results as queue.
                    tasks.put(copy_second_put_data)
                    i[2] = 0
    except KeyError as err:
        print(err)
    # if target_a not found will return None value so need to remove


def search_unequal_bbh_data(target_b):
    """Search the unequal backward best hit data.
    ex) AAE_CAC_backward_best_hit
    No Return !! vlaues means Null return"""
    # print(unequal_BBH_data)
    for i in unequal_BBH_data:
        if i[2] != 0 and (target_b[0] == i[0] or target_b[0] == i[1]\
            or target_b[1] == i[0] or target_b[1] == i[1]):
            # sample of B  - ["AAE_gi|156|ref....","CAC_gi|1560..",85]
            copy_i = copy.copy(i)
            tasks.put(copy_i)
            unequal_BBH_data[unequal_BBH_data.index(i)][2] = 0


def matching_bbh(target):
    """ Match the backward best hit """
    # print("running matching_bbh")
    if target[2] == 0:
        return  # !!! Why Return Function

    else:
        copy_target = copy.copy(target)
        search_equal_bbh_data(copy_target[0])
        search_equal_bbh_data(copy_target[1])
        results.put(copy_target)
        unequal_BBH_data[unequal_BBH_data.index(target)][2] = 0

    for j in unequal_BBH_data:
        if j[2] == 0:
            pass  # ! Remove this later
        else:
            if copy_target[0] == j[0] or copy_target[0] == j[1] \
                or copy_target[1] == j[0] or copy_target[1] == j[1]:
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
    "This Function perform the matrix_clustering_ortholog"
    # matrix_clustering_ortholog(element_set, bar):
    # # this is Old Method old is generating is  removed
    # """ Generate the matrix of clustering ortholog. """
    # element_set is a list collection gene id and # ! score
    row_data = []
    col_data = []
    temp_results = queue.Queue()
    # Create a queue object with a given maximum size.
    # If maxsize is <=0, teh queue size is infinite
    for element in element_set:
        # example of element
        # ["AAE_gi|15605623|ref|NP_212996.1", "AAE_gi|15605623|ref|NP_212996.1",998]
        # ["AAE_gi|15004781|ref|NP_149241.1|","CAC_gi|15004781|ref|NP_149241.1|",9520]
        # if element[0] exist, returning the index in the row_data.
        if row_data.count(element[0]) > 0:
            # element[0] is data of row.
            # ['gi|15606057|ref|NP_213434.1|', 'gi|15606057|ref|NP_213434.1|', '3823\n']
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
    # np.zeros() can be used but as ndarray
    # will test later,matlib function return matrix,np.matrix(np.zeros((len,width)))

    while not temp_results.empty():
        get_temp_results = temp_results.get()
        row = get_temp_results[0]
        col = get_temp_results[1]
        score_matrix[row, col] = get_temp_results[2]
        score_matrix[col, row] = get_temp_results[2]

    # If the elements of row and col is less than 2, it is excluded.
    if len(row_data)*len(col_data) > 4:
        # The big size of matrix(bigger than 1000 X 1000)
        # will be computed by parallel_matrix_multiplication function.
        if len(row_data) > 1000 and cpu_count > 1:
            score_matrix = parallel_mcl(score_matrix)
        else:
            score_matrix = mcl(score_matrix)  # !!! Error
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
        for power in power_results:
            sum_matrix = power + sum_matrix

        # create a inflation_matrix(part 2)
        pool = multiprocessing.Pool(cpu_count)
        divide_results = pool.map(parallel_matrix_division,
                                  zip(power_results, repeat(sum_matrix)))
        pool.close()
        pool.join()

        # Make a Combined matrix for results of parallel_matrix_multiplication function.
        for i, _ in enumerate(divide_results):
            if i == 0:
                score_matrix = divide_results[i]
            else:
                score_matrix = np.concatenate(
                    (score_matrix, divide_results[i]), axis=0)
        # numpy.concatenate(a1,a2...)sequence of array_like,if axis not passed then

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
                  % ((mcl_time_finish - mcl_time_start)/60,\
                      count, score_matrix[0].size, score_matrix[0].size))
    return score_matrix


def mcl(score_matrix):
    """Input score_matrix is the user selected matrix llike BLOASUM45, BLOSUM62 or BLOSUM80"""
    #print("mcl funciton running")
    count = 0
    infinitesimal_value = 10**-10
    idempotent_matrix = numpy.matlib.ones((2, 2))
    # "idempotent_matrix = np.ones((2,2)) # i will later test with this"
    while idempotent_matrix.sum() > infinitesimal_value:  # > infinitesimal_value
        mcl_time_start = time.time()
        expansion_matrix = score_matrix ** 2
        # print("shape of input score_matrix", score_matrix.shape)
        # print(expansion_matrix)

        score_matrix = np.power(expansion_matrix, float(
            inflation_factor))  # !!! Eroor
        # print(score_matrix)
        score_matrix_sum = score_matrix.sum(axis=0)
        # print(score_matrix_sum,"line 983")
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
                  % ((mcl_time_finish - mcl_time_start)/60,\
                      count, score_matrix[0].size, score_matrix[0].size))
    return score_matrix


def clustering(row_data, col_data, score_matrix):
    """
    This Function will cluster the Matrix
    """
    global CLUSTER_COUNT, ortholog_count
    ortholog_temp_list = []
    for i, _ in enumerate(row_data):
        ortholog_list = []
        ortholog_sum = 0  # !! Not used
        ortholog = queue.Queue()  # It is Queue which is put the ortholog.
        # It is Queue which is put the ortholog having changed gene ID
        gene_id_queue = queue.Queue()
        for j, _ in enumerate(col_data):
            if score_matrix[i, j] >= 0.1:
                ortholog.put(col_data[j])
                ortholog_list.append(col_data[j])

        # If the ortholog queue  has the element of ortholog more than 2, it will be printed.
        if ortholog.qsize() >= 3:
            for element in ortholog_list:
                try:
                    ortholog_sum += ortholog_temp_list.index(element)+1

                except ValueError:
                    with open(CLUSTER_OUT+"_geneID_S"+str(threshold_score) +
                              "_"+str(inflation_factor), "a") as ortholog_list_save:
                        ortholog_print = "cluster "+str(CLUSTER_COUNT)+" :"
                        ortholog_list_save.write(ortholog_print)

                        while not ortholog.empty():
                            get_ortholog = ortholog.get()
                            ortholog_list_save.write("\t"+get_ortholog)

                            try:
                                # get_ortholog --> ECO_170082288,get_ortholog.split('_')
                                # --> ['ECO', '170082288']
                                get_ortholog_split = get_ortholog.split('_')
                                gene_id_queue.put(
                                    get_ortholog_split[0]+"_"+gene_id_dic[get_ortholog_split[1]])
                            except KeyError:
                                # If the gene_id_dic don't have get_ortholog,
                                # it will print the original ID(get_ortholog).
                                gene_id_queue.put(get_ortholog)
                        ortholog_list_save.write("\n")
                    with open(CLUSTER_OUT+"_KO_ID_S"\
                              + str(threshold_score)\
                              + "_"+str(inflation_factor), "a")\
                            as ortholog_list_geneID_save:
                        ortholog_list_geneid_print = "cluster " + \
                            str(CLUSTER_COUNT)+" :"
                        ortholog_list_geneID_save.write(
                            ortholog_list_geneid_print)
                        while not gene_id_queue.empty():
                            ortholog_list_geneID_save.write(
                                "\t"+gene_id_queue.get())
                            ortholog_count += 1
                        CLUSTER_COUNT += 1
                        ortholog_list_geneID_save.write("\n")
                    break
            ortholog_temp_list = operator.concat(
                ortholog_temp_list, ortholog_list)


def parallel_matrix_multiplication(data_):
    "Doing parallel_matrix_multiplication Using python Numpy"
    matrix_element, matrix = data_
    return matrix_element * matrix


def parallel_matrix_power(matrix_element):
    """This Function Compute Parallel_Matrix Power by using Numpy np.power() Function
    this will do element wise operation, power of each  element given by input
    """
    return np.power(matrix_element, float(inflation_factor))


def parallel_matrix_division(data_):
    """This Function Will Perform Matrix Division by Using Numpy library
    This Function is also do element wise operation on given matrix with sum data
    """
    matrix_element, sum_data = data_
    return np.divide(matrix_element, sum_data)


def read_species(pr=1):  # default 1 which will always shows the name of SPECIES
    """ If pr is 1, it will print "SPECIES_List"  Other wise only return the Value in dic Format """
    read_species_ = os.listdir(COMMAND_OPTIONS.SPECIES)
    selected_species_dic = {}
    backward_selected_species_ = {}
    for i, species in enumerate(sorted(read_species_), start=1):
        selected_species_dic[i] = species
        backward_selected_species_[species] = i
        if pr == 1:
            print(str(i)+".", species)
        number = i
    # number for total length of species
    return selected_species_dic, backward_selected_species_, number


def del_file(path, file):
    "This Function delete the file passed with path and files or folders\
        This will replaced with lambda Function "
    try:
        os.remove(path+file)
        print("File Successfully Removed")
        # "to remove all files inside the folder we need to use loop over the folder"
    except IOError as err:
        print("Check the File or Path", err)


def check_file(file):
    """When mode 3 is Selected.This function will run.
    file is user input value for clustering output.
    Check the file weather exists or not .
    if exist this will warn to use other name
    CLUSTER_OUT is file.CLUSTER_OUT(results)
    + _geneID_S+ str(threshold_score)
    +"_"+str(inflation_factor) or CLUSTER_OUT+"_KO_ID_S"+str(threshold_score)
    +"_"+str(inflation_factor))
    are files created by mcl algorithm so if these file exist system will exit
    """
    # Declare all variable as a globally Added
    # global CLUSTER_OUT, threshold_score, infinite_loop
    # list all related files in same directory local variable
    file_list = glob.glob(file + '*')
    if (CLUSTER_OUT+"_geneID_S"+str(threshold_score)
            + "_"+str(inflation_factor) or CLUSTER_OUT+"_KO_ID_S"
            + str(threshold_score)+"_"+str(inflation_factor)) in file_list:
        print("Please, set other name,%s is of output.line Number 966")
        sys.exit(2)
    else:
        print(file, "is ")


def read_equal_bbh(path):
    """Read_Equal BBH by user path Blast Best Hit
    the 'path = score' fiie location is ./score_file/ + species + species
    like ./score_file/A_C
    # !! Since the path i.e score file is always same so no neeed to pass path
    this function open genome_genome_oneway_tbh_score5 and
    genome_genome_second_ oneway_ threshold_best_hit_score5 and
    and convert these file to list inside list converting score value as int

    """
    with open(path+"_oneway_tbh_Score"
              + str(threshold_score), 'r') as equal_RBH:
        for _ in equal_RBH:
            split_data = _.split() # !!
            # Split Data sample
            # ['A_gi|15605613|ref|NP_212986.1|', 'C_gi|15004707|ref|NP_149167.1|', '51']
            split_data[2] = int(split_data[2])  # Converting to int format split_data[2] = 51
            equal_BBH_data.append(split_data)   # Equal_BBH is list created before running this
            equal_BBH_data_dic[split_data[0]] = split_data[1:]
    try:
        with open(path+"_second_oneway_tbh_Score"
                  + str(threshold_score), 'r') as second_equal_RBH:
            for _ in second_equal_RBH:
                split_data = _.split()
                # "['AAE_gi|15606128|ref|NP_213505.1|', 'AAE_gi|15606877|ref|NP_214257.1|', 898]
                split_data[2] = int(split_data[2])  # Change value to integer type
                # split_data[2] = 898 in the above case
                second_equal_BBH_data.append(split_data)
    except IOError as err:
        print(err,"Skipping")


def read_unequal_bbh(path):
    """ Read unequal BBH path passed by User.
    Since this Function will append value to unequal_BBH_data
    which created before running this function
    in case of read_unequal_bbh the only one file genome_genome_"""

    with open(path+"_oneway_tbh_Score"
              + str(threshold_score), 'r') as unequal_RBH:
        for j in unequal_RBH:
            split_data = j.split()
            split_data[2] = int(split_data[2])
            unequal_BBH_data.append(split_data)

def write_info():
    """  #
    with open("./function_samples/same_species_forward_best_hit(score_file).txt",'a+') as bp_score:
        bp_score.write(str(species_of_query)+str(species_of_subject))
        bp_score.write("\n\nTHis is Input file for same species_forward_best_hit()")
        bp_score.write("\n")
        bp_score.write(blastp_score)
        bp_score.write("*"*20+"\n")
        bp_score.write("This is same_species_ forward_best_score\n")
        bp_score.write("\n".join(str(x) for x in same_species_ forward_best_score))
        bp_score.write("\n"+"*"*20+"\n")
    """

print("="*82)
print("||\t****Default Variables(Values)****\t\t\t\t\t||")
print("||\tBlastp                 = %s\t\t\t\t\t\t||" % COMMAND_OPTIONS.BLASTP)
print("||\tBLASTP_DATA            = %s\t\t\t\t\t||" %
      COMMAND_OPTIONS.BLASTP_DATA)
print("||\tblstp_matrix           = %s\t\t\t\t\t||" %
      COMMAND_OPTIONS.blastp_matrix)
print("||\tcpu_count              = %s\t\t\t\t\t\t||" %
      COMMAND_OPTIONS.cpu_count)
print("||\tgenomes                = %s\t\t\t\t\t\t||" %
      COMMAND_OPTIONS.genomes)
print("||\tinfinite_loop          = %s\t\t\t\t\t\t||" %
      COMMAND_OPTIONS.infinite_loop)
print("||\tinflation_factor       = %s\t\t\t\t\t\t||" %
      COMMAND_OPTIONS.inflation_factor)
print("||\tMODE                   = %s\t\t\t\t\t\t||" % COMMAND_OPTIONS.MODE)
print("||\tCLUSTER_OUT            = %s\t\t\t\t\t||" %
      COMMAND_OPTIONS.CLUSTER_OUT)
print("||\tthreshold_score        = %s\t\t\t\t\t\t||" %
      COMMAND_OPTIONS.threshold_score)
print("||\tsave_raw_blastp_score  = %s\t\t\t\t\t\t||" %
      COMMAND_OPTIONS.save_raw_blastp_score)
print("||\tScore_file             = %s\t\t\t\t\t||" %
      COMMAND_OPTIONS.Score_file)
print("||\tSPECIES                = %s\t\t\t\t\t||" % COMMAND_OPTIONS.SPECIES)
print("||\tverbose                = %s\t\t\t\t\t\t||" %
      COMMAND_OPTIONS.verbose)
print("="*82)
print("\n")
print("Ortholog Detection Program Starts Now\n")

if not sys.argv[1:]:
    # If not Parameter Passed the Manual Process will Start,
    # in case of MODE the input parameter can be different
    # so we will not apply the below code function
    print("1. BLASTP. \n2. BLASTP using precalculated data. \n3. clustering.\n")
    MODE = input(">> Select a MODE or MODEs space seprater\n\
        (1 3 ***OR*** 2 3): ").split(" ")
    check_mode(MODE) # ! if 1 or 2 is not Passed the system will exist
    # ! Created by Krishna !!
    SELECTED_SPECIS, backward_selected_species, number_i = read_species(
        1)  # read_species(1) will return dictionary type
    SELECTED_NUMBER = input(
        ">> Select Genomes to detect Orthologs(e.g. 1 2 3 4 5 or 1-5) : ")
    # !! later we will use lambda function to conrol un usual number
    if SELECTED_NUMBER.find('-') > 0:
        # find() return the index position of first occurance
        SN = SELECTED_NUMBER.split("-")
        if int(SN[-1]) > number_i:
            # number_i is length of Genome file inside folder exit the process
            print("\nWrongInput\nInput must be less than", number_i)
            sys.exit(2)
        else:
            USER_SELECTED_NUMBER = range(int(SN[0]), int(SN[-1])+1)
            for j in USER_SELECTED_NUMBER:
                print(SELECTED_SPECIS[j], end=" ")  # loop in Dic
            print("Selected!!")

    else:
        USER_SELECTED_NUMBER = sorted(
            set([int(_) for _ in SELECTED_NUMBER.split()]))
        # Create a set (remove repeating)
        if int(USER_SELECTED_NUMBER[-1]) > number_i:
            print("\nWrongInput\nInput must be less than", number_i)
            sys.exit(2)
            # Greater than Genome list will system error
        else:
            for j in USER_SELECTED_NUMBER:
                print(SELECTED_SPECIS[j], end=" ")
            print("Selected!!")
    blastp_matrix = matrix_name()
    cpu_count = int(input("You can use %s Core.\nIf you input >= 2,\
        The Program will run a parallel computation for the blastp.\n"\
            % multiprocessing.cpu_count()\
            + "Enter the number of Core to use in this program (1 ~ %s): "\
                % multiprocessing.cpu_count()))

elif sys.argv[1:]:
    genomes = COMMAND_OPTIONS.genomes
    MODE = COMMAND_OPTIONS.MODE
    cpu_count = COMMAND_OPTIONS.cpu_count
    blastp_matrix = COMMAND_OPTIONS.blastp_matrix
    inflation_factor = COMMAND_OPTIONS.inflation_factor
    SELECTED_SPECIS, backward_selected_species, number_i = read_species()
    USER_SELECTED_NUMBER = [backward_selected_species[ele]
                            for ele in genomes]  # select by value of dictionary
    CLUSTER_OUT = COMMAND_OPTIONS.CLUSTER_OUT     # File Name to save
    check_mode(MODE) # ! If 1 or 2 not Passed system will exist

SPECIES = COMMAND_OPTIONS.SPECIES
BLASTP = COMMAND_OPTIONS.BLASTP
Score_file = COMMAND_OPTIONS.Score_file
BLASTP_DATA = COMMAND_OPTIONS.BLASTP_DATA
save_raw_blastp_score = COMMAND_OPTIONS.save_raw_blastp_score
threshold_score = COMMAND_OPTIONS.threshold_score
verbose = COMMAND_OPTIONS.verbose
infinite_loop = COMMAND_OPTIONS.infinite_loop
log_time = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")

# print(COMMAND_OPTIONS.SPECIES ,COMMAND_OPTIONS.Blastp ,COMMAND_OPTIONS.Score_file)

if "3" in MODE:
    # Mode 3 is for clustering the matrix data.
    # print("3 MODE clustering is selected")
    inflation_factor = input("Enter the inflation factor to cluster: ")
    CLUSTER_OUT = input("Set the name of clustering output Folder : ")
    """CLUSTER_OUT is user input value --> the name of clustering output <--results
    """
    check_file(CLUSTER_OUT)
    # if same file for user input species and cluester out and inflation factor is same
    # no need to run cluster again

# del_file(Score_file, "*")
# this Function will delete all files inside currently folder

if "3" in MODE:
    # for log file location
    Log_file_name = "./cluster_out/" + CLUSTER_OUT + "_S"\
        + str(threshold_score) + "_" \
        + str(inflation_factor) + log_time+"_log.txt"

elif "3" not in MODE:
    # if 3 is not passed the log file will store inside "log_files"
    Log_file_name = 'log_files/Log' + log_time+"_log.txt"

with open(Log_file_name, 'w') as log:
    log.write(str(datetime.datetime.now()))
    log.write("\nMODE :")
    for i in MODE:
        log.write(" " + i)
    log.write("\ngenomes : ")
    for i in USER_SELECTED_NUMBER:
        log.write(SELECTED_SPECIS[i] + " ")
    log.write("\ncpu_count : " + str(cpu_count))
    log.write("\nblastp matrix : " + blastp_matrix)
    if "3" in MODE:
        log.write("\ninflation_factor : " + str(inflation_factor))
        log.write("\nCluster out : " + CLUSTER_OUT)
    log.write("\nSPECIES : " + SPECIES)
    log.write("\nBLASTP : " + BLASTP)
    log.write("\nScore file : " + Score_file)
    log.write("\nBLASTP_DATA : " + BLASTP_DATA)
    log.write("\nsave rawblastp score : " + str(save_raw_blastp_score))
    log.write("\n")

start_time_OBH = time.time()
if "1" in MODE:    # 1 is for Blastp Run
    bbh_work_list = oneway_tbh(MODE)
    # backward _best_hit_work_list is changed to bbh_ work_list
    # print("Back Ward Best Hit result", bbh_work_list)
    if bbh_work_list !=[]:
        pool = multiprocessing.Pool(cpu_count)
        results = pool.map(backward_best_hit, bbh_work_list)
        # "function_name is backward_best_hit and the parameter is bbh_work_list"
        pool.close()
        pool.join()
    else:
        results = [0,0]

elif "2" in MODE:
    used_precalculated_data_list = []
    # the values are appended from oneway_tbh(MODE)
    new_calculated_data_list = []
    precalculated_data_list = glob.glob(
        BLASTP_DATA + "*oneway_tbh_Score" + str(threshold_score))
    # _S changed to Score file name and copy score file
    # from score_file to blastp_file folder manually to speed up the system
    # glob.glob function may not work properly in windows so  need to check the  result
    bbh_work_list = oneway_tbh(MODE)
    # This function also add elements to used_precalculated_data_list
    # If bbh_work_list is an empty list, pool instance can't finsh the work.
    if bbh_work_list != []:
        # If bbh_work_list is an empty list, pool instance can't finsh the work.
        # "The pool only run if backward best_hit_work_list has some value."
        pool = multiprocessing.Pool(cpu_count)
        # multiprocessing.pool() for Parallel
        results = pool.map(backward_best_hit, bbh_work_list)
        # backward_best_hit is the function and bbh_work_list is list parameter
        pool.close()
        pool.join()
    else:
        results = [0, 0]
# del_file("./", "query*")
# Delete all files start with query
# "Since the  query files will replaced by new files"

finish_time_OBH = time.time()
blastp_time_log = float(((finish_time_OBH - start_time_OBH)/60))
print("BLASTP searches + forward best Hit + backwardbest hit took %f minutes" %
      blastp_time_log)
# Del_file("./","query*") #delete all files of query
with open(Log_file_name, 'a') as log:
    log.write("backward_best_hit took " + str(max(results)) + "minutes\n")
    # results will not be empty because of pool not show any results then results = [0,0]
    log.write("BLASTP + Best_Hit + backward_best_hit searches took " +
              str(blastp_time_log) + " minutes\n")

if "3" in MODE:
    start_time_clustering = time.time()
    ##########################################################################################
    # generate matrix and calculate the matrix using mcl algorithm and cluster the ortholog."""
    """queue is a linear data structure that stores items in First In First Out(FIFO) manner.
    Create a queue object with a given maximum size.If maxsize is <=0, the queue size is infinite.
    the default
    queue.Queue(maxsize = 0)
    """
    print("\n>>>> Start mcl algorithm and clustering ortholog <<<<")
    # "the values of these blank list and dic is added inside the function which is very bad"
    equal_BBH_data = []  # store value as list inside list
    unequal_BBH_data = []  # store value as list inside list
    equal_BBH_data_dic = {}
    second_equal_BBH_data = []
    results = queue.Queue()
    tasks = queue.Queue()
    CLUSTER_COUNT = 1
    ortholog_count = 0
    gene_id_dic = {}

    with open("./db/myva=gb", "r") as id_read:
        # myvba=gb is a database
        for i in id_read:
            gene_name, gene_id = i.split()
            # remove "\n",and create a new dic
            gene_id_dic[gene_id.replace("\n", "")] = gene_name

    if "1" in MODE:
        # The file location is different in each MODE
        # in MODE 2 files are instide blastp_data but
        # in MODE 1 files are inside score_file

        for i in USER_SELECTED_NUMBER:
            for k in USER_SELECTED_NUMBER:
                if k < i:
                    # print(k,"is less than ", i)
                    pass
                elif i == k:
                    # If the both species is same, then read_equal_bbh() function will run
                    read_equal_bbh(
                        Score_file + SELECTED_SPECIS[i] + "_" + SELECTED_SPECIS[k])
                elif i != k:
                    # If both the species are different then read_unequal_bbh() function will run
                    read_unequal_bbh(
                        Score_file + SELECTED_SPECIS[i] + "_" + SELECTED_SPECIS[k])

    elif "2" in MODE:
        # The file location is different in each MODE
        # MODE 2 is for fast precalcualted data sets
        # and datas are inside the blastp_data"
        for used_data in used_precalculated_data_list:
            # used_precalculated_data_list is return from function one_way_threshold_best_hit(MODE)
            first, second = used_data.split("_")
            if first == second:
                # if Both genes are same
                read_equal_bbh(BLASTP_DATA + used_data)
            elif first != second:
                read_unequal_bbh(BLASTP_DATA + used_data)

        for new_data in new_calculated_data_list:
            # "new_calculated_data_list is lilst appended in
            # Function Oneway_threshold_best_hit(MODE =2) "
            first, second = new_data.split("_")
            if first == second:
                read_equal_bbh(Score_file + new_data)
            elif first != second:
                read_unequal_bbh(Score_file + new_data)

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
        # matrix_clustering_ortholog(data, bar)
    finish_time_clustering = time.time()
    MCL_TIME_LOG = float((finish_time_clustering - start_time_clustering)/60)
    REMARK_TIME_LOG = float((finish_time_clustering - start_time_OBH)/60)
    print("mcl algorithm and Ortholog clustering took %.2f minutes" % MCL_TIME_LOG)
    print("owPRemark program took %.2f minutes" % REMARK_TIME_LOG)

    if "3" in MODE:
        with open(Log_file_name, 'a') as log:
            log.write("Ortholog count : " + str(ortholog_count) + "," +
                      " Cluster count : " + str(CLUSTER_COUNT-1) + "\n")
            log.write("mcl algorithm and Ortholog clustering took " +
                      str(MCL_TIME_LOG) + " minutes\n")
            log.write("XXX program took " + str(REMARK_TIME_LOG) + "minutes\n")
print("All Code Success Fully executed")
# pylint: disable = import-self,invalid-name,unused-argument, too-many-lines,
