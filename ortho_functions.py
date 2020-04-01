#!/usr/bin/python3
import sys
import os  # For system operation
import glob
import subprocess
import time
#from progress.bar import Bar
from alive_progress import alive_bar
import multiprocessing
import queue   # For Multiprocessing
import numpy.matlib
import numpy as np     # For Mathmatical (Algebra) Operation
import copy
import operator
import pprint
from itertools import *
import argparse
import datetime
from sys import platform   # To verify the Operating Sys Win or Linux
from Bio import SeqIO  # For Bio Python Sequence Object


def Matrix_Name():
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


def Query_Sequence(genome):
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
            #print("Query_Sequence() Run Successfully")
            return gene_seq_list

    except IOError as err:
        print("IOError occurred in Query_Sequence function : " + str(err))


def Query_Seq(genome):
    """This Function Read Fasta (Genome) file and return as a list Format with gene Position and Sequence Developed by Krish"""
    try:
        return [str((seq_record.id+seq_record.seq)) for seq_record in SeqIO.parse(genome, "fasta")]
        # To understand this Function Bio Python library needs to be studied
    except IOError as err:
        print(str(err))


def Write_Query(query, parallel_num):
    "This Function Write Query with file Name query+ parallel_num in same directory and raise IO error if Error rises"
    # print("Write_Query(query,parallel_num)")
    try:
        with open("./query/query_"+str(parallel_num), "w") as write_query:
            write_query.write(query)
    except IOError as err:
        print("IOError occurred in Write_Query function : " + str(err))


def Run_Blast(subject, parallel_num):
    """By this Function it will create a Pipe line to run Blastp in Computer by input Parameter Subject is whole Genome
         and parallel_number is query file Created by early step"""

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


def Same_Species_Forward_Best_HIt(blastp_score):
    """Search the forward best hit among the blastp scores of same species.
    Because there are an duplicated genes in a same genome.
    blstp_score file is return value of Run_Blast() Function
    When the blastp score compare with blastp score of duplicate gene, if score and length are same, blasp score of duplicated gene is added to a second best score."""
    #print("Rnunning Get_Same_Sp0ecies_Forward_Best_Hit")
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
        if k[0] == k[1]:
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

    return best_score


def Forward_Best_HIt(blastp_score):
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
        # print ">>>>>>>>>>>>>>>Forward_Best_HIt   k", k
        if int(k[2]) > int(temp_best_score[2]):  # Compare the Score of the Blast
            temp_best_score = k
        elif int(k[2]) == int(temp_best_score[2]):
            if int(k[3]) > int(temp_best_score[3]):  # Compare the Lenth of Sequence
                temp_best_score = k
            elif int(k[3]) == int(temp_best_score[3]):
                second_temp_best_score.append(k)
 #    print "############ temp best score ############", temp_best_score
    best_score.append(temp_best_score)
    for j in second_temp_best_score:
        if j[2] == temp_best_score[2] and j[3] == temp_best_score[3]:
            best_score.append(j)
    return best_score, blastp_score_split_list


def Division_Parallel_Query(queryV, query_division_value, cpu_count, queryV_len):
    "queryV is list of queries, query_division_value is , cpu_count is the Number of CPU User Entered and queryV_len is the length of the query"
    parallel_query = []
    parallel_query_start = 0
    if queryV_len % cpu_count == 0:  # perfect division
        for i in range(cpu_count):
            i += 1
            if parallel_query_start == 0:
                parallel_query.append(
                    queryV[int(parallel_query_start):i*int(query_division_value)])
                parallel_query_start += 1
            else:
                parallel_query.append(queryV[int(
                    parallel_query_start)*int(query_division_value):i*int(query_division_value)])
                parallel_query_start += 1
    else:  # imperfect division
        for i in range(cpu_count):
            i += 1
            if parallel_query_start == 0:
                parallel_query.append(
                    queryV[int(parallel_query_start):i*int(query_division_value)])
                parallel_query_start += 1
            elif i < cpu_count:
                parallel_query.append(queryV[int(
                    parallel_query_start)*int(query_division_value):i*int(query_division_value)])
                parallel_query_start += 1
            elif i == cpu_count:
                parallel_query.append(
                    queryV[int(parallel_query_start)*int(query_division_value):])
                parallel_query_start += 1

    return parallel_query


def Run_Parallel_Query(species_of_query, species_of_subject, queryV, parallel_num):
    """ Run_Parallel_Query(i ,k , queryV , cpu_count) i and k are user selected number in a list Format
    Run the following functions. Write_Query, Run_Blast, Same_Species_Forward_Best_HIt, Forward_Best_HIt
    Save the files which are oneway_threshold_best_hit, second_oneway_threshold_best_hit, 
    blastp_score_split_list and raw_blastp_score (optional) by each species. 
    parallel_num is the Number of CPU selected by the user if CPu 1 selected then 1 """

    #print("Run_Parallel_Query Running")
    global selected_number, selected_species_dic

    # bar = Bar('Processing '+str(parallel_num), max = len(queryV)) #progressing bar setting , Creating a Object
    # bar is not Supported in Python 3    We use alive_bar instead
    with alive_bar(len(queryV)) as bar:  # declare your set of items for loop
        for j in queryV:
            i = 10
            # bar.next() #progressing bar not supported
            bar()  # Call after Consuming One Item
            Write_Query(j, parallel_num)  # if 1 Added Here Add also to Run
            Write_Query(j, i)
            i += 1  # !Only For checking Delete later
        # This Function Only Write a file with j name and parallel_num i.e CPU Count
            blastp_score = Run_Blast(
                selected_species_dic[species_of_subject], parallel_num)
        # Return Byte File type
            if blastp_score != '':  # Check whether blastp_score has the value
                best_score, blastp_score_split_list = Forward_Best_HIt(
                    blastp_score.decode())  # .decode() Convert in to str format
                # ex) AAE == AAE. It will save best_score without reversing Run_Blast.
                if species_of_query == species_of_subject:
                    same_species_forward_best_score = Same_Species_Forward_Best_HIt(
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
                else:  # If species_of_query not equal with species_of_subject, run reversing Run_Blast
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

            if save_raw_blastp_score:
                with open(Score_file+selected_species_dic[species_of_query]
                          + "_"+selected_species_dic[species_of_subject]+"_S"
                          + str(threshold_score)+"_"+str(parallel_num), "a") as save_blastp:
                    save_blastp.write(blastp_score)
    # bar.finish() # progressing bar finish
    return  # None


def Oneway_Threshold_Best_Hit(mode):
    """ This Function accept the mode and Run program and return backward_best_hit_work_list  """
    #print("OneWay_Threshold_Best_Hit running")
    global user_selected_number, cpu_count
    process_list = []
    backward_best_hit_work_list = []
    if "1" in mode:
        print("Blastp Mode is Selected")
        # We have 3 Mode 1 is for blastp
        # Mode 2 is for BLASTP using precalcualted data and Mode 3 is for Clustering
        for i in user_selected_number:  # Select species to write query
            queryV = Query_Sequence(selected_species_dic[i])
            # queryV is a list Format  with a position gene id , Seq
            # User Selected  [1, 3, 5] is list of User input
            for k in user_selected_number:
                if k < i:  # gene ====> query 1->1 1->2 1->3 2->2 2->3  Forward checking will skip the same or less
                    continue
                else:
                    print("Doing the blastp & forward best hit searches between %s genome and %s genome"
                          % (selected_species_dic[i], selected_species_dic[k]))
                    # this will Create a loop though 2 times

                    # length of Genome gene ID in Sequence File
                    queryV_len = len(queryV)
                    if cpu_count == 1:
                        print("cput_count =", cpu_count)
                        blastp_time_start = time.time()
                        Run_Parallel_Query(i, k, queryV, cpu_count)
                        blastp_time_end = time.time()
                        print("The blastp & forward best hit searches took %.2f minutes" % (
                            (blastp_time_end-blastp_time_start)/60))
                    else:
                        # If the number of queryV_len is less than cpu_count, Remark will select the number of queryV_len.
                        if queryV_len < cpu_count:
                            blastp_time_start = time.time()
                            # 1 is query_division_value. Because queryV_len / queryV_len(=cpu_count) is 1.
                            parallel_query = Division_Parallel_Query(
                                queryV, 1, queryV_len, queryV_len)
                            for m in range(queryV_len):
                                process = multiprocessing.Process(
                                    target=Run_Parallel_Query, args=(i, k, parallel_query[m], m+1))
                                # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                process_list.append(process)
                                process.start()
                            for n in process_list:
                                n.join()
                            blastp_time_end = time.time()
                            print("The blastp & forward best hit searches took %.2f minutes" % (
                                (blastp_time_end-blastp_time_start)/60))
                        else:
                            blastp_time_start = time.time()
                            query_division_value = queryV_len / cpu_count
                            parallel_query = Division_Parallel_Query(
                                queryV, query_division_value, cpu_count, queryV_len)
                            for m in range(cpu_count):
                                process = multiprocessing.Process(
                                    target=Run_Parallel_Query, args=(i, k, parallel_query[m], m+1))
                                # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                process_list.append(process)
                                process.start()
                            for n in process_list:
                                n.join()
                            blastp_time_end = time.time()
                            print("The blastp & forward best hit searches took %.2f minutes" % (
                                (blastp_time_end-blastp_time_start)/60))
                if not i == k:
                    backward_best_hit_work_list.append((i, k, queryV_len))
                    print("Backward")

    elif "2" in mode:
        for i in user_selected_number:  # Select species to write query
            queryV = Query_Sequence(selected_species_dic[i])
            for k in user_selected_number:  # Select of subject
                if Blastp_data+selected_species_dic[i]+"_"+selected_species_dic[k]\
                    + "_oneway_threshold_best_hit_Score"\
                        + str(threshold_score) in precalculated_data_list:
                    used_precalculated_data_list.append(
                        selected_species_dic[i]+"_"+selected_species_dic[k])
                    continue
                else:
                    if k < i:  # gene ====> query 1->1 1->2 1->3 2->2 2->3
                        continue

                    else:
                        print("Doing the blastp & forward best hit searches between %s genome and %s genome" % (
                            selected_species_dic[i], selected_species_dic[k]))

                        queryV_len = len(queryV)
                        if cpu_count == 1:
                            blastp_time_start = time.time()
                            Run_Parallel_Query(i, k, queryV, cpu_count)
                            blastp_time_end = time.time()
                            print("The blastp & forward best hit searches took %.2f minutes" % (
                                (blastp_time_end-blastp_time_start)/60))
                        else:
                            # If the number of queryV_len is less than cpu_count, Remark will select the number of queryV_len.
                            if queryV_len < cpu_count:
                                blastp_time_start = time.time()
                                # 1 is query_division_value. Because queryV_len / queryV_len(=cpu_count) is 1.
                                parallel_query = Division_Parallel_Query(
                                    queryV, 1, queryV_len, queryV_len)
                                for m in range(queryV_len):
                                    process = multiprocessing.Process(
                                        target=Run_Parallel_Query,
                                        args=(i, k, parallel_query[m], m+1))
                                    # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                    process_list.append(process)
                                    process.start()
                                for n in process_list:
                                    n.join()
                                blastp_time_end = time.time()
                                print("The blastp & forward best hit searches took %.2f minutes" % (
                                    (blastp_time_end-blastp_time_start)/60))
                            else:
                                blastp_time_start = time.time()
                                query_division_value = queryV_len / cpu_count
                                parallel_query = Division_Parallel_Query(
                                    queryV, query_division_value, cpu_count, queryV_len)
                                for m in range(cpu_count):
                                    process = multiprocessing.Process(
                                        target=Run_Parallel_Query,
                                        args=(i, k, parallel_query[m], m+1))
                                    # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                    process_list.append(process)
                                    process.start()
                                for n in process_list:
                                    n.join()
                                blastp_time_end = time.time()
                                print("The blastp & forward best hit searches took %.2f minutes" % (
                                    (blastp_time_end-blastp_time_start)/60))
                        new_calculated_data_list.append(
                            selected_species_dic[i]+"_"+selected_species_dic[k])
                        if not i == k:
                            backward_best_hit_work_list.append(
                                (i, k, queryV_len))
    return backward_best_hit_work_list


def Backward_Best_Hit(args):
    print("Running Backward_Best_Hit")
    species_of_query, species_of_subject, queryV_len = args
    start_time_BBH = time.time()
    forward_best_hit_score_list = []
    blastp_score_split_list = []
    print("Run the backward best hit between %s genome %s genome" % (
        selected_species_dic[species_of_query], selected_species_dic[species_of_subject]))
    # If the number of queryV_len is less than cpu_count, the cpu_count is changed to queryV_len.
    if queryV_len < cpu_count:
        for parallel_num in range(queryV_len):
            parallel_num += 1
            with open(Score_file+selected_species_dic[species_of_query]+"_"
                      + selected_species_dic[species_of_subject] + "_"+"best_score_S"
                      + str(threshold_score)+"_"+str(parallel_num), "r") as best_hit_score:
                for each_line in best_hit_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
        #            print split_each_line
                    forward_best_hit_score_list.append(split_each_line)
            with open(Score_file+selected_species_dic[species_of_query]+"_"
                      + selected_species_dic[species_of_subject] +
                      "_"+"blastp_score_split_list_S"
                      + str(threshold_score)+"_"+str(parallel_num), 'r') as blastp_score:
                for each_line in blastp_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
        #            print split_each_line
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
        #            print split_each_line
                    forward_best_hit_score_list.append(split_each_line)
            with open(Score_file+selected_species_dic[species_of_query] + "_"
                      + selected_species_dic[species_of_subject]
                      + "_"+"blastp_score_split_list_S"
                      + str(threshold_score)+"_"+str(parallel_num), 'r') as blastp_score:
                for each_line in blastp_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
        #            print split_each_line
                    blastp_score_split_list.append(split_each_line)

    #bar = Bar("Searching : "+selected_species_dic[species_of_query]+"-"+selected_species_dic[species_of_subject], max = len(forward_best_hit_score_list))

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
    finish_time_BBH = time.time()
    RBH_time = float((finish_time_BBH - start_time_BBH)/60)
    print("BackwardBestHit of %s-%s took %.2f minutes" %
          (selected_species_dic[species_of_query], selected_species_dic[species_of_subject], RBH_time))
    return RBH_time


def Search_Equal_BBH_Data(target_A):
    """Search the equal backward best hit data. ex) AAE_AAE_backward_best_hit """
    # print("Search_Equal_BBH_Data")
    put_data = equal_BBH_data_dic[target_A]
    if put_data[1] == 0:
        pass
    else:
        copy_put_data = copy.copy(put_data)
        copy_put_data.insert(0, target_A)
        results.put(copy_put_data)
        equal_BBH_data_dic[target_A][1] = 0
    #        print "---put zero in Search_Equal_BBH_Data----- ", equal_BBH_data_dic[target_A]
        for i in second_equal_BBH_data:
            if i[2] == 0:
                pass
            else:
                if i[0] == target_A or i[1] == target_A:
                    copy_second_put_data = copy.copy(i)
                    # Don't put results as queue. Because the tasks will put copy_second_put_data to results as queue.
                    tasks.put(copy_second_put_data)
                    i[2] = 0
    #                    print "---put zero in second_equal_BBH_data---", i
    return


def Search_Unequal_BBH_Data(target_B):
    """Search the unequal backward best hit data. ex) AAE_CAC_backward_best_hit"""
    print("Search_Unequal_BBH_Data")
    for i in unequal_BBH_data:
        if i[2] == 0:
            pass
        else:
            if target_B[0] == i[0] or target_B[0] == i[1] or target_B[1] == i[0] or target_B[1] == i[1]:
                copy_i = copy.copy(i)
                tasks.put(copy_i)
                unequal_BBH_data[unequal_BBH_data.index(i)][2] = 0
    return # None


def Matching_BBH(target):

    """ Match the backward best hit """
    #print("running Matching_BBH")
    if target[2] == 0:
        return

    else:
        copy_target = copy.copy(target)
        Search_Equal_BBH_Data(copy_target[0])
        Search_Equal_BBH_Data(copy_target[1])
        results.put(copy_target)
        unequal_BBH_data[unequal_BBH_data.index(target)][2] = 0

    for j in unequal_BBH_data:
        if j[2] == 0:
            pass

        else:
            if copy_target[0] == j[0] or copy_target[0] == j[1] or copy_target[1] == j[0] or copy_target[1] == j[1]:
                copy_j = copy.copy(j)
                # print "targ_get , j = %s %s" % (copy_target, j)
                unequal_BBH_data[unequal_BBH_data.index(j)][2] = 0
                tasks.put(copy_j)

    while not tasks.empty():
        get_task = tasks.get()
        Search_Equal_BBH_Data(get_task[0])
        Search_Equal_BBH_Data(get_task[1])
        results.put(get_task)
        Search_Unequal_BBH_Data(get_task)

    return #None 


def Generating_Matrix_Clustering_Ortholog(element_set):
    # Generating_Matrix_Clustering_Ortholog(element_set, bar): # this is Old Method
    """ Generate the matrix of clustering ortholog. """
    print("Generating_Matrix_Clustering_Ortholog")
    row_data = []
    col_data = []
    temp_results = queue.Queue()
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
        # The big size of matrix(bigger than 1000 X 1000) will be computed by Parallel_Matrix_Multiplication function.
        if len(row_data) > 1000 and cpu_count > 1:
            score_matrix = Parallel_MCL(score_matrix)
        else:
            score_matrix = MCL(score_matrix)
        Clustering(row_data, col_data, score_matrix)


def Parallel_MCL(score_matrix):
    # print("Parallel_MCL")
    count = 0
    infinitesimal_value = 10**-10
    idempotent_matrix = numpy.matlib.ones((2, 2))

    while idempotent_matrix.sum() > infinitesimal_value:  # > infinitesimal_value
        MCL_time_start = time.time()
        pool = multiprocessing.Pool(cpu_count)  # create a expansion_matrix
        multiplication_results = pool.map(Parallel_Matrix_Multiplication,
                                          zip(score_matrix, repeat(score_matrix)))
        pool.close()
        pool.join()

        # create a inflation_matrix(part 1)
        pool = multiprocessing.Pool(cpu_count)
        power_results = pool.map(
            Parallel_Matrix_Power, multiplication_results)
        pool.close()
        pool.join()

        sum_matrix = 0
        for i in power_results:
            sum_matrix = i + sum_matrix

        # create a inflation_matrix(part 2)
        pool = multiprocessing.Pool(cpu_count)
        divide_results = pool.map(Parallel_Matrix_Division,
                                  zip(power_results, repeat(sum_matrix)))
        pool.close()
        pool.join()

        # Make a Combined matrix for results of Parallel_Matrix_Multiplication function.
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
        if count > infinite_loop:  # It will prevent the infinite loop of MCL algorithm.
            break
        MCL_time_finish = time.time()
        if verbose:
            print(" MCL time : %f, count : %d, matrix size : %d * %d"
                  % ((MCL_time_finish - MCL_time_start)/60, count, score_matrix[0].size, score_matrix[0].size))
    return score_matrix


def MCL(score_matrix):
    count = 0
    infinitesimal_value = 10**-10
    idempotent_matrix = numpy.matlib.ones((2, 2))
    #idempotent_matrix = np.ones((2,2))
    while idempotent_matrix.sum() > infinitesimal_value:  # > infinitesimal_value
        MCL_time_start = time.time()
        expansion_matrix = score_matrix ** 2
        score_matrix = np.power(expansion_matrix, inflation_factor)
        score_matrix_sum = score_matrix.sum(axis=0)
        # create a inflation_matrix
        score_matrix = np.divide(score_matrix, score_matrix_sum)
        # identify whether inflation_matrix is idempotent matrix or not.
        idempotent_matrix = abs(score_matrix - expansion_matrix)
        count += 1
        if count > infinite_loop:  # It will prevent the infinite loop of MCL algorithm.
            break
        MCL_time_finish = time.time()
        if verbose:
            print(" MCL time : %f, count : %d, matrix size : %d * %d"
                  % ((MCL_time_finish - MCL_time_start)/60, count, score_matrix[0].size, score_matrix[0].size))
    return score_matrix


def Clustering(row_data, col_data, score_matrix):
    global cluster_count
    global ortholog_count
    ortholog_temp_list = []
    for i in range(len(row_data)):
        ortholog_list = []
        ortholog_sum = 0
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
                        ortholog_list_geneID_print = "cluster " + \
                            str(cluster_count)+" :"
                        ortholog_list_geneID_save.write(
                            ortholog_list_geneID_print)
                        while not gene_id_queue.empty():
                            ortholog_list_geneID_save.write(
                                "\t"+gene_id_queue.get())
                            ortholog_count += 1
                        cluster_count += 1
                        ortholog_list_geneID_save.write("\n")
                    break
            ortholog_temp_list = operator.concat(
                ortholog_temp_list, ortholog_list)


def Parallel_Matrix_Multiplication(data):
    "Doing Parallel_Matrix_Multiplication Using Numpy"
    matrix_element, matrix = data
    result = matrix_element * matrix
    return result


def Parallel_Matrix_Power(matrix_element):
    "This Function Compute Parallel_Matrix Power by using Numpy np.power() Function"
    power_matrix_element = np.power(matrix_element, inflation_factor)
    return power_matrix_element


def Parallel_Matrix_Division(data):
    "This Function Will Perform Matrix Division by Using Numpy library"
    matrix_element, sum_data = data
    return np.divide(matrix_element, sum_data)


def Read_Species(pr=1):  # default 1 which will always shows the name of Species
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


def Del_File(path, file):
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
    except:
        print("Check the File or Path")


def Check_File(file):
    "Check the file weather exist or not  This will run after Mode 3 is selected"
    # Declare all variable as a globally Added
    global Cluster_out, threshold_score, infinite_loop
    file_list = glob.glob(file+'*') #list all related files
    #file_list = os.listdir(file)
    if (Cluster_out+"_geneID_S"+str(threshold_score)
        + "_"+str(inflation_factor) or Cluster_out+"_KO_ID_S"
            + str(threshold_score)+"_"+str(inflation_factor)) in file_list:
        print("Please, set other name of output.")
        sys.exit(2)


def Read_Equal_BBH(path):
    "Read_Equal BBH by user path Blast Best Hit"
    global threshold_score

    with open(path+"_oneway_threshold_best_hit_Score"
              + str(threshold_score), 'r') as equal_RBH:
        for j in equal_RBH:
            split_data = j.split()
            split_data[2] = int(split_data[2])
            equal_BBH_data.append(split_data)
            equal_BBH_data_dic[split_data[0]] = split_data[1:]
    try:
        with open(path+"_second_oneway_threshold_best_hit_Score"
                  + str(threshold_score), 'r') as second_equal_RBH:
            for j in second_equal_RBH:
                split_data = j.split()
                split_data[2] = int(split_data[2])
                second_equal_BBH_data.append(split_data)
    except:
        pass


def Read_Unequal_BBH(path):
    "Read unequal BBH path passed by User"
    with open(path+"_oneway_threshold_best_hit_Score"
              + str(threshold_score), 'r') as unequal_RBH:
        for j in unequal_RBH:
            split_data = j.split()
            split_data[2] = int(split_data[2])
            unequal_BBH_data.append(split_data)



if __name__ == "__main__":
    Matrix_Name()                                   #1
    Query_Sequence()                                  #2
    Write_Query()                                        #3
    Run_Blast()                                          #4
    Same_Species_Forward_Best_Hit()                 #5
    Forward_BestHit()                                 #6
    Division_Parallel_Query()                             #7
    Parallel_Query()                                  #8
    Oneway_Threshold_Best_Hit()                         #10
    Backward_Best_Hit()                                 #11
    Search_Equal_BBH_Data()                             #12
    Search_Unequal_BBH_Data()                           #13
    Matching_BBH()                                      #14
    Generating_Matrix_Clustering_Ortholog()             #15
    Parallel_MCL()                                      #16
    MCL()                                               #17
    Clustering()                                        #18
    Parallel_Matrix_Multiplication()        #19
    Parallel_Matrix_Division()                #20
    Read_Species(1)                                 #21
    Del_File()                                    #22
    Check_File()                                        #23
    Read_Equal_BBH()                                    #2
    Read_Unequal_BBH()                                  #26
