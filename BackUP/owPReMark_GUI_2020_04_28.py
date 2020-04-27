
import random
import webbrowser # To browse help Menu
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from time import strftime
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

selected_species = ""

# !! All Functions for this Program
# ! Function for Browsing Folder
def browse_folder(arg):
    """Allow user to select a directory and store it in global variable variable_name --folder_path."""
    fold_name = filedialog.askdirectory()
    arg.set(fold_name+"/")
    return arg

# ! Functrion for browsing File
def browse_file(arg):
    "Allow user to select a file path and store in a global variable "
    file_name = filedialog.askopenfilename() #(filetyopes = [("Fastaa file),".faa"])
    arg.set(file_name.split("/")[-1]) # ! For windows

# ! Function for Help Menu
def openURL(x):
    """Open the help menu of Git Hub."""
    new = 2
    webbrowser.open(x,new = new)

# ! Function to display the information about Company & Developer
def info():
    """Show the information about company & developer."""
    tk.messagebox.showinfo("info","Theragen Genome Care, South Korea,\n\
                           concept developed by : Kim\n\
                        Developer: Adhikari Krishna,\n\
                        Supporting bak, okg")

# ! Function to exit Program with warning
def quit():
    """Exit the program."""
    msg = tk.messagebox.askquestion("Exit Application ?","Are you sure to exit the application ?",icon ="warning")
    if msg == "yes":
        window.destroy()
    else:
        tk.messagebox.showinfo("Returning to original stage","Values for variable didnot changed")

# ! Function to Reset all the Variables
def reset():
    """Reset the defaults parameters."""
    global CLUSTER_COUNT, ortholog_count
    msg = tk.messagebox.askquestion("Restore defaults ?",\
        "Are you sure to set default values for all variables ?",icon ="warning")
    if msg =="yes":
        blastp_data. set("./blastp_folder/")
        score_file.set("./score_file/")
        threshold_score_.set(5)
        inflation_factor_.set(1.5)
        log_file_.set("log")
        cluster_out.set("results")
        cpu_count_.set(1)
        var_in.set("Click Check Variable")
        species_folder.set("")
        blastp_folder_.set("")
        subject_file_.set(" ")
        query_file_.set(" ")
        result_display.set(">>>>>>")
        CLUSTER_COUNT = 1
        ortholog_count = 0
    else:
        tk.messagebox.showinfo("Return","you will now return to the application screen")

def display():
    """Display the information about variable."""
    selected = "Threshold Score is -->"+str(threshold_score_.get())\
        +"\nInflation Factor is -->"+str(inflation_factor_.get())\
        +"\nScoreFile is -->"+str(score_file.get())\
        +"\nLog File is -->"+log_file_.get()\
        +"\nCPU selected -->"+str(cpu_count_.get())\
        +"\nMetrix selection-->"+ str(matrix_selection)\
        +"\nSpecies Folder-->"+species_folder.get()\
        +"\nFlastp_files_folder-->"+blastp_folder_.get()\
        +"\nSubject_Genome-->"+subject_file_.get()\
        +"\nQuery_genome-->"+query_file_.get()
    var_in.set(selected)

def show_choice():
    global matrix_selection
    matrix_selection = bp_matrix.get()


# ! Functions for  Blast running
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
    global SPECIES
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
    return run_blastp.stdout


def same_species_forward_best_hit(blastp_score):
    """Search the forward best hit among the blastp scores of
    same species.Because there are duplicated genes in a same genome.
    blstp_score file is return value of run_blast() Function
    When the blastp score compare with blastp score of duplicate gene,
    if score and length are same,
    blasp score of duplicated gene is added to a second best score."""
    # print(blastp_score)
    blastp_score_split_list = []
    temp_best_score = ['-1', '-1', '-1']
    second_temp_best_score = []
    best_score = []
    second_best_score = []
    blastp_score_split = blastp_score.split("\n")
    # delete of [''] in the last index  ex) ['gi,gi,1,1','gi,gi,2,2','']
    del blastp_score_split[-1]
    for _ in blastp_score_split:
        blastp_score_element = _.split(',')
        blastp_score_split_list.append(blastp_score_element)
    # ex) k is ['gi|15605613|ref|NP_212986.1|', 'gi|15605613|ref|NP_212986.1|', '3702', '699']
    for k in blastp_score_split_list:
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
    """Search the forward best hit among the blastp scores of same species."""
    #  print("Rnunning GetForward BestHit")
    blastp_score_split_list = []
    temp_best_score = ['-1', '-1', '-1']
    second_temp_best_score = []
    best_score = []
    blastp_score_split = blastp_score.split("\n")
    # delete of ['']   ex) ['gi,gi,1,1','gi,gi,2,2','']
    del blastp_score_split[-1]
    for _ in blastp_score_split:
        blastp_score_element = _.split(',')
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
    best_score.append(temp_best_score)
    for j in second_temp_best_score:
        if (j[2] == temp_best_score[2]) and (j[3] == temp_best_score[3]):
            best_score.append(j)
    return best_score, blastp_score_split_list


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
    Save the files which are oneway_threshold_best_hit,
    second_oneway_threshold_best_hit,
    blastp_score_split_list and raw_blastp_score(optional) by each species.
    parallel_num is the Number of CPU selected by the user if CPu 1 selected then 1

    """

    with alive_bar(len(query_v)) as bar_:
        # declare your set of items for loop, bar is blacklisted
        for j in query_v:
            # query_v is the list of the genes from genome
            write_query(j, parallel_num)
            # if 1 Added Here Add also to Run
        # This Function Only Write a file with j name and parallel_num i.e CPU Count
            blastp_score = run_blast(
                selected_species[species_of_subject], parallel_num)
        # Return Byte File type, and this is only for one gene position query and Whole Genome
            if blastp_score != '':  # if blastp run successfully
                best_score, blastp_score_split_list = forward_best_hit(
                    blastp_score)  # .decode() Convert in to str format
                # ex) AAE == AAE. It will save best_score without reversing run_blast.
                # The sample for blastp_score and blastp_score_split_list is saved inside folder

                if species_of_query == species_of_subject:
                    same_species_forward_best_score = same_species_forward_best_hit(
                        blastp_score)  # .decode()
                    for best_score_element in same_species_forward_best_score:
                        # ex) [A1 of AAE, A1 of AAE, 30]
                        if best_score_element[0] == best_score_element[1]:
                            # if same species add write file in # !!!
                            with open(Score_file + selected_species[species_of_query]\
                                + "_" + selected_species[species_of_subject]\
                                + "_oneway_threshold_best_hit_Score"\
                                + str(threshold_score), "a") as oneway_threshold_bh:
                                # best_hit is replaced with h
                                save_best_score = selected_species[species_of_query]\
                                    + "_" + best_score_element[0]\
                                    .split(r"\s")[0]\
                                    + " "\
                                    + selected_species[species_of_subject]\
                                    + "_"+best_score_element[1]\
                                    .split(r"\s")[0]+" " + best_score_element[2]+"\n"
                                # best_score_element[0].split("|")
                                # ==> ['gi','15642790', 'ref', 'NP_227831.1', '']
                                oneway_threshold_bh.write(
                                    save_best_score)
                        else:  # ex) [A1 of AAE, A2 of AAE, 30]
                            with open(Score_file + selected_species[species_of_query]\
                                      +"_"\
                                      + selected_species[species_of_subject]\
                                      + "_second_oneway_threshold_best_hit_Score"\
                                      + str(threshold_score), "a") as second_oneway_threshold_bh:
                                      # best_hit is renamed with bh
                                second_save_best_score = selected_species[species_of_query]\
                                    + "_"+best_score_element[0].split(r"\s")[0]\
                                    + " " + selected_species[species_of_subject]\
                                    + "_" + best_score_element[1].split(r"\s")[0]\
                                    + " "+best_score_element[2]+"\n"
                                second_oneway_threshold_bh.write(
                                    second_save_best_score)
                else:
                    # If species_of_query not equal with species_of_subject, run reversing run_blast
                    for best_score_element in best_score:
                        if '-1' not in best_score_element:   # ! if not "-1" in best_score_element:
                            with open(Score_file + selected_species[species_of_query]\
                                    + "_"\
                                    + selected_species[species_of_subject]\
                                    + "_"+"best_score_S"\
                                    + str(threshold_score)\
                                    +"_"\
                                    +str(parallel_num), "a") as save_best_hit:
                                best_score_save = selected_species[species_of_query]\
                                    + "_"\
                                    + best_score_element[0].split(r"\s")[0]\
                                    + " "+selected_species[species_of_subject]\
                                    + "_"+best_score_element[1].split(
                                        r"\s")[0]\
                                    + " " + best_score_element[2]\
                                    + " " + best_score_element[3]+"\n"
                                save_best_hit.write(best_score_save)

                    for blastp_score_split_list_element in blastp_score_split_list:
                        with open(Score_file+selected_species[species_of_query]+"_"
                                  + selected_species[species_of_subject]
                                  + "_"+"blastp_score_split_list_S"
                                  + str(threshold_score)
                                  + "_"+str(parallel_num), "a") as save_blastp_score_split_list:

                            blastp_score_split_list_save = selected_species[species_of_query]\
                                    + "_"\
                                    + blastp_score_split_list_element[0].split(\
                                        r"\s")[0]\
                                    + " "+ selected_species[species_of_subject]\
                                    + "_"+blastp_score_split_list_element[1].split(r"\s")[0]\
                                    + " " + blastp_score_split_list_element[2]\
                                    +" " + blastp_score_split_list_element[3]+"\n"
                            save_blastp_score_split_list.write(
                                blastp_score_split_list_save)
            else:
                print("--------")
            bar_()  # Call after Consuming One Item
            if save_raw_blastp_score:
                with open(Score_file+selected_species[species_of_query]
                          + "_"+selected_species[species_of_subject]+"_S"
                          + str(threshold_score)+"_"+str(parallel_num), "a") as save_blastp:
                    # Only string is format is supported
                    save_blastp.write(blastp_score)  # .decode()
    # bar.finish() # progressing bar finish


def oneway_threshold_best_hit():
    """ This Function accept the mode and Run program and return
     bbh_work_list mode 1 and 2 is supported
     python continue_loop """

    info = "Running the blastp & forward best hit searches "
    process_list = []
    bbh_work_list_ = []
    print("Running oneway threshold best_hit")
    for i,genome_i in enumerate(selected_species):  # Select species to write query
        query_v = query_sequence(genome_i) # name changed as species

        for k,_ in enumerate(selected_species):
            # ! changed by Krish
            # ! while can be used but create a infinity loop
            if k < i:
                continue
            else:
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

    return bbh_work_list_


def oneway_tbh():
    global precalculated_data_list , result_display
    print("Running oneway_tbh")
    process_list = []
    bbh_work_list_ = []
    info = "Running the blastp & forward best hit searches "
    print("\nMode 2, Running from precalculated data sets which will save the Time")
    # print(len(precalculated_data_list),"oneway_threshold_best_hit_ files are located")
    # # !!! Delete This row later
    for i, genome_i in enumerate(selected_species):  # Select species to write query
        query_v = query_sequence(genome_i)
        for k,genome_k in enumerate(selected_species):  # Select of subject
            if BLASTP_DATA + genome_i +"_"+ genome_k\
                + "_oneway_threshold_best_hit_Score"\
                    + str(threshold_score) in precalculated_data_list:
                print("Skipped")
                used_precalculated_data_list.append(
                    genome_i + "_"+ genome_k)
            else:
                result_display.set(result_display.get()+"Precalcualted Blast File Not Found\nRunning Blast again\Which is time Consuming")
                if k < i:  # gene ====> query 1->1 1->2 1->3 2->2 2->3
                    continue
                else:
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

                                process_list.append(process)
                                process.start()
                            for n in process_list:
                                n.join()
                            blastp_time_end = time.time()
                            print("**"*20)
                            print(type(info) , type(blastp_time_end) , type(blastp_time_end))
                            print("**"*20)
                            print(info + " took %.2f minutes" % ((blastp_time_end - blastp_time_start)/60))
                    new_calculated_data_list.append(
                        genome_i+"_"+ genome_k)
                    if not i == k:
                        bbh_work_list_.append(
                            (i, k, query_v_len))
    return bbh_work_list_


def backward_best_hit(args):
    """ This will  calculate backward_best_hit with
    This function open 2 files from score_files folder
    Genome_genome_best_score_S_5_parallel_num and
    Genome_Genome_blastp_score_split_list_S_5_parallel_num
    and convert last 2 digits to integer values.
    and again loopdone in
    forward_best_hit_score_list and
    blastp_score_split_list

    """
    bck_info = "Running the backward_best_hit"
    # print(bck_info)
    species_of_query, species_of_subject, query_v_len = args
    start_time_bbh = time.time()
    forward_best_hit_score_list = []
    blastp_score_split_list = []
    print(bck_info + " between %s genome %s genome" % (
        selected_species[species_of_query], selected_species[species_of_subject]))
    # If the number of query_v_len is less than cpu_count, the cpu_count is changed to query_v_len.
    if query_v_len < cpu_count:
        for parallel_num in range(query_v_len):
            parallel_num += 1
            with open(Score_file + selected_species[species_of_query]+"_"
                      + selected_species[species_of_subject]\
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
            with open(Score_file+selected_species[species_of_query]+"_"
                      + selected_species[species_of_subject] +
                      "_"+"blastp_score_split_list_S"
                      + str(threshold_score)+"_"+str(parallel_num), 'r') as blastp_score:
                    # AAE_gi|15605613|ref|NP_212986.1| CAC_gi|15004754|ref|NP_149214.1| 80 55
                    # AAE_gi|15605613|ref|NP_212986.1| CAC_gi|15004707|ref|NP_149167.1| 69 14
                for each_line in blastp_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
                    blastp_score_split_list.append(split_each_line)

    else:
        # if cpu length is smaller than queryV length
        for parallel_num in range(cpu_count):
            # cpu_count is the number of cpu used by the user
            # ! Error may raise if parallel Num is not in the file
            parallel_num += 1
            with open(Score_file+selected_species[species_of_query]+"_"
                      + selected_species[species_of_subject]+"_"+"best_score_S"
                      + str(threshold_score)+"_"+str(parallel_num), "r") as best_hit_score:
                for each_line in best_hit_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[3] = int(split_each_line[3])
                    forward_best_hit_score_list.append(split_each_line)
            with open(Score_file+selected_species[species_of_query] + "_"
                      + selected_species[species_of_subject]
                      + "_"+"blastp_score_split_list_S"
                      + str(threshold_score)+"_"+str(parallel_num), 'r') as blastp_score:
                for each_line in blastp_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
                    blastp_score_split_list.append(split_each_line)

    # "bar = Bar("Searching : "+SELECTED_SPECIS_DIC[species_of_query]
    # +"-"+SELECTED_SPECIS_DIC[species_of_subject], max = len(forward_best_hit_score_list)) ""

    for forward_best_hit_score_element in forward_best_hit_score_list:
        matching_list = []
        backward_best_score = ['-1', '-1', '-1']
        # length is 3 bcz the length is 3 in
        # example of blastp_score_split_list
        # [AAE_gi|15605613|ref|NP_212986.1| ECO_gi|170082858|ref|YP_001732178.1| 3020 703,
        # AAE_gi|15605614|ref|NP_212987.1| ECO_gi|170083440|ref|YP_001732760.1| 1990 404 ]
        # botth blastp_score_split_list and forward_best_hit_score_list shape are almost same
        # for element in blastp_score_split_list:
        #    if element[1] == forward_best_hit_score_element[1]:
        #        matching_list.append(element)
        # ! The upper 3 for loop be converted list Comprehensive
        matching_list = [element for element in blastp_score_split_list\
            if element[1]== forward_best_hit_score_element[1]]

        for element in matching_list:
            if int(element[2]) > int(backward_best_score[2]):
                backward_best_score = element


        if int(backward_best_score[2]) - int(forward_best_hit_score_element[2]) <= threshold_score:
            # the default threshold_score passed is 5
            with open(Score_file+selected_species[species_of_query]\
                      + "_"+selected_species[species_of_subject]\
                      + "_oneway_threshold_best_hit_Score"+str(threshold_score), "a")\
                    as other_oneway_threshold_best_hit:
                save_data = forward_best_hit_score_element[0]\
                    + " "+forward_best_hit_score_element[1]\
                    + " " + str(int(forward_best_hit_score_element[2]))+"\n"
                other_oneway_threshold_best_hit.write(save_data)

    # bar.finish()
    finish_time_bbh = time.time()
    rbh_time = float((finish_time_bbh - start_time_bbh)/60)
    print(bck_info + " %s-%s took %.2f minutes" %
          (selected_species[species_of_query],\
              selected_species[species_of_subject], rbh_time))
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

    global results , equal_BBH_data_dic, second_equal_BBH_data

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
    global unequal_BBH_data
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
    global unequal_BBH_data , results
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
    global results
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
    global CLUSTER_COUNT, ortholog_count, gene_id_dic
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
                                    get_ortholog_split[0]+"_" + gene_id_dic[get_ortholog_split[1]])
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
    this function open genome_genome_oneway_threshold_best_hit_score5 and
    genome_genome_second_oneway_threshold_best_hit_score5 and
    and convert these file to list inside list converting score value as int

    """
    with open(path+"_oneway_threshold_best_hit_Score"
              + str(threshold_score), 'r') as equal_RBH:
        for _ in equal_RBH:
            split_data = _.split() # !!
            # Split Data sample
            # ['A_gi|15605613|ref|NP_212986.1|', 'C_gi|15004707|ref|NP_149167.1|', '51']
            split_data[2] = int(split_data[2])  # Converting to int format split_data[2] = 51
            equal_BBH_data.append(split_data)   # Equal_BBH is list created before running this
            equal_BBH_data_dic[split_data[0]] = split_data[1:]
    try:
        with open(path+"_second_oneway_threshold_best_hit_Score"
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

    with open(path+"_oneway_threshold_best_hit_Score"
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
        bp_score.write("This is same_species_forward_best_score\n")
        bp_score.write("\n".join(str(x) for x in same_species_forward_best_score))
        bp_score.write("\n"+"*"*20+"\n")
    """


def main(mode_):
    global selected_species, blastp_matrix, cpu_count, blastp_folder, SPECIES, BLASTP, Score_file, BLASTP_DATA, threshold_score,\
        inflation_factor, CLUSTER_OUT, Log_file_name
    global equal_BBH_data , unequal_BBH_data, equal_BBH_data_dic, second_equal_BBH_data, matched_BBH_data,\
        matched_BBH_element_data_set, tasks, results, gene_id_dic, CLUSTER_COUNT, ortholog_count, used_precalculated_data_list,\
            precalculated_data_list, new_calculated_data_list

    #print( query_file_.get(), selected_species_, blastp_matrix, SPECIES, Score_file, BLASTP_DATA, threshold_score,\
    #    blastp_matrix, blastp_folder, inflation_factor)
    selected_species = [subject_file_.get(), query_file_.get()]
    blastp_matrix = bp_matrix.get()
    cpu_count = cpu_count_.get()
    blastp_folder = blastp_folder_.get() # Path for Precalculated Blastp Folder
    SPECIES = species_folder.get()
    BLASTP = "blastp"
    Score_file = "./score_file/"
    BLASTP_DATA = blastp_folder_.get()
    threshold_score = threshold_score_.get()
    log_time = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
    inflation_factor = inflation_factor_.get()
    CLUSTER_OUT = cluster_out.get()
    Log_file_name = "./log_files/"+CLUSTER_OUT + "_S"+ str(threshold_score) + "_"\
        + str(inflation_factor) + "_" + log_time[6:] + log_file_.get() + ".txt"

    # log file name is so long so reducing the length
    with open(Log_file_name, 'w+') as log:
        log.write(str(datetime.datetime.now()))
        log.write(query_file_.get())
        log.write("\ngenomes : ")
        log.write("\n"+str(selected_species))
        log.write("\ncpu_count : " + str(cpu_count))
        log.write("\nblastp matrix : " + blastp_matrix)
        log.write("\ninflation_factor : " + str(inflation_factor))
        log.write("\nCluster out : " + CLUSTER_OUT)
        log.write("\nSPECIES : " + SPECIES)
        log.write("\nBLASTP : " + BLASTP)
        log.write("\nScore file : " + Score_file)
        log.write("\nBLASTP_DATA : " + BLASTP_DATA)
        log.write("\nsave rawblastp score : " + blastp_folder)
        log.write("\n")
    # ! XXXXXXX
    #input_ = input("NUMber to exit")
    # if input_ == "1":
    #    sys.exit()
    start_time_OBH = time.time()
    if mode_ in  ["blastp_only","blastp+mcl"]:
        # blastp is common for all mode
        bbh_work_list = oneway_threshold_best_hit() # blast running
        if bbh_work_list !=[]:
            pool = multiprocessing.Pool(cpu_count)
            results = pool.map(backward_best_hit, bbh_work_list)
        # "function_name is backward_best_hit and the parameter is bbh_work_list"
            pool.close()
            pool.join()
        else:
            results = [0,0] # not passing blank value so
    elif mode_ == "mcl_blastp_data":
        # print("you clicked mcl with blast data")
        used_precalculated_data_list = [] # ! global variable
        new_calculated_data_list = []
        precalculated_data_list = glob.glob(BLASTP_DATA + "*_oneway_threshold_best_hit_Score"+str(threshold_score))
        if precalculated_data_list:
            print("data exists")
        # !!! XXXX
        print(len(precalculated_data_list))
        print("*"*20)
        # ! XXXX
        bbh_work_list = oneway_tbh()  # oneway_threshold_best_hit(MODE)
        if bbh_work_list != []:
            pool = multiprocessing.Pool(cpu_count)
            results = pool.map(backward_best_hit, bbh_work_list)
            pool.close()
            pool.join()
        else:
            results = [0,0]
        # ! del_file("./query/", "query*")
    # if mode_ == "blastp_only":
    #    print("Blastp run success")
    #    sys.exit(2)

    finish_time_OBH = time.time()
    blastp_time_log = float(((finish_time_OBH - start_time_OBH)/60))
    print("BLASTP searches + forward best Hit + backwardbest hit took %f minutes" %
        blastp_time_log)

    #! Del_file("./","query*") #delete all files of query
    with open(Log_file_name, 'a') as log:
        log.write("backward_best_hit took " + str(max(results)) + "minutes\n")
        # results will not be empty because of pool not show any results then results = [0,0]
        log.write("BLASTP + Best_Hit + backward_best_hit searches took " +
                str(blastp_time_log) + " minutes\n")

    start_time_clustering = time.time()
    results = queue.Queue()
    tasks = queue.Queue()

    if mode_ in ["mcl_blastp_data","blastp+mcl"]:
        with open("./db/myva=gb", "r") as id_read:
            for i in id_read:
                gene_name, gene_id = i.split()
                # remove "\n",and create a new dic
                gene_id_dic[gene_id.replace("\n", "")] = gene_name

    if mode_ == "blastp+mcl":
            for index_i,i in enumerate(selected_species):
                for index_k,k in enumerate(selected_species):
                    if index_k < index_i:
                        pass
                    elif index_i == index_k:
                        read_equal_bbh(
                                Score_file + selected_species[index_i] + "_" + selected_species[index_k])
                    elif i != k:
                        read_unequal_bbh(
                            Score_file + selected_species[index_i] + "_" + selected_species[index_k])

    if mode_ == "mcl_blastp_data":
        for used_data in used_precalculated_data_list:
            first, second = used_data.split("_")
            if first == second :
                read_equal_bbh(BLASTP_DATA + used_data)
            elif first != second:
                read_unequal_bbh(BLASTP_DATA + used_data)
        for new_data in new_calculated_data_list:
            first, second = new_data.split("_")
            if first == second :
                read_equal_bbh(Score_file + new_data)
            elif first != second:
                read_unequal_bbh(Score_file + new_data)

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
    print("MCL algorithm and Ortholog clustering took %.2f minutes" % MCL_TIME_LOG)
    print("owPRemark program took %.2f minutes" % REMARK_TIME_LOG)

    with open(Log_file_name, 'a') as log:
        temp_ = "Ortholog count : " + str(ortholog_count) + "," +\
            " Cluster count : " + str(CLUSTER_COUNT-1) + "\n"
        result_display.set(temp_)
        temp_ = "MCL algorithm and Ortholog clustering took " +\
                    str(MCL_TIME_LOG)[:10] + " minutes\n"
        log.write(temp_)
        result_display.set(result_display.get()+temp_)
        temp_ = "owPReMark program took " + str(REMARK_TIME_LOG)[:10] + "minutes\n"
        log.write(temp_)
        result_display.set(result_display.get()+temp_)

    result_display.set(result_display.get() + "\n\n\nAll Code Success Fully Run\
        \nHave a great day")

## !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
window = tk.Tk()
window.title("owPReMark_Orthologs_Detection_Software")
window.geometry("1000x800+100+20")
window.resizable(0,0)
# photo = tk.PhotoImage(file = "./icons/python.gif")
# window.tk.call("wm","iconphoto",window._w,photo)
window.title("owPReMark, Software to detect Orthologes among Genomes")

# !! ______________Creating Different Frames_better_view_______________
tk.Frame(window, width=1000,height =10,bg="green").pack(side="top")
# top and buttom Green line
tk.Frame(window,width=1000,height =10,bg="green").pack(side="bottom")
f1 = tk.Frame(window,width = 800,height = 100,bg = "yellow")
f1.pack(side = "bottom") # This hold Click Buttons

# !!!  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Frame 2 for extra Features
tk.Frame(window,width = 500,height = 200,bg= "#FEF8DD",bd =15,).place(x=490,y=510)
tk.Frame(window,width = 500,height = 350,bg= "#8DE4FD").place(x=490,y=180)

# !!!!  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
tk.Label(window,font = ("arial",40,"bold"),text = "owPReMark V1.0",\
    fg = "white",bg="#093145").pack()

# !! Current Time
local_time=time.asctime(time.localtime(time.time()))
tk.Label(window,font = ("arial",18,"bold"),\
                      text ="Started-->  "+local_time,fg = "#093145", bd =10).pack()
local_time=time.asctime(time.localtime(time.time()))
tk.Label(window,font = ("arial",18,"bold"),text ="Total --> "+local_time,fg = "#093145",bd =10).pack() # time_label
# !time_label will display Time Consumed after Function clicked

# !_______________________help & Credit Button__________________
tk.Button(window,text="Help....", width=10,padx=5,pady=8,bd=8,font=('arial',16,"bold"),bg="powder blue",\
    command= lambda : openURL("https://www.google.com/")).place(x =800,y = 20)

tk.Button(window,text = "Credit", width=10,padx=1,pady=8,bd=8,font=('arial',16,"bold"),bg="powder blue",command = info).place(x =20,y = 20)

## !!XXXXXXXXXXX__Global_Variables__XXXXXXXXXXXXXXXXXXXXXXX
# Global variables
cpu_count_ = tk.IntVar() ; cpu_count_.set(1)
threshold_score_ =  tk.DoubleVar() ; threshold_score_.set(5)
inflation_factor_ = tk.DoubleVar() ; inflation_factor_.set(1.5)
blastp_data = tk.StringVar() ; blastp_data. set("./blastp_data/")
score_file  = tk.StringVar() ; score_file.set("./score_file/")
subject_file_ = tk.StringVar()
log_file_ = tk.StringVar() ; log_file_.set("log")
cluster_out = tk.StringVar() ; cluster_out.set("results")
query_file_ = tk.StringVar()
matrix_selection = tk.StringVar()
species_folder = tk.StringVar()
blastp_folder_ = tk.StringVar(); blastp_folder_.set("./blastp_data/")
var_in = tk.StringVar() ; var_in.set("...........")
bp_matrix = tk.StringVar(); bp_matrix.set("BLOSUM45")
result_display = tk.StringVar(); result_display.set(">>>>")

# !! Main Variables XXXXXX
selected_species_ = [subject_file_.get(), query_file_.get()] 
blastp_matrix = "blosum45"
cpu_count = cpu_count_.get()
blastp_folder = blastp_folder_.get() # Path for Precalculated Blastp Folder
SPECIES = species_folder.get()
BLASTP = "blastp"
Score_file = "./score_file/"
BLASTP_DATA = blastp_folder_.get()
threshold_score = threshold_score_.get()
save_raw_blastp_score = ".score_file/"
log_time = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
inflation_factor = inflation_factor_.get()
CLUSTER_OUT = "result"
infinite_loop = 60
verbose = False
Log_file_name = ""
equal_BBH_data = []
unequal_BBH_data = []  # store value as list inside list
equal_BBH_data_dic = {}
second_equal_BBH_data = []
matched_BBH_data = []
matched_BBH_element_data_set = []
CLUSTER_COUNT = 1
ortholog_count = 0
gene_id_dic = {}
tasks = ""
results = ""
used_precalculated_data_list = []
precalculated_data_list = []
new_calculated_data_list = []

# !_____________Create 3 button inside Fram1_

tk.Button(f1,padx=16,pady=15,bd=10,fg="green",\
    font=('arial',12,"bold"),width=10,text="Blastp_Only",\
    activebackground = "#63ace5",bg="powder blue", command = lambda : main("blastp_only")).pack(side = "left")

tk.Button(f1,padx=16,pady=15,bd=10,bg = "#0057e7",fg="white",\
    font=('arial',12,"bold"),width=12,text="Blastp+MCL",activebackground = "#63ace5",\
    highlightcolor="blue",command = lambda : main("blastp+mcl")).pack(side = "left")


tk.Button(f1,padx=16,pady=15,bd=10,fg="white",\
    font=('arial',12,"bold"),width=20,text="MCL_with_blastp_data",\
        bg="#0057e7",activebackground = "#63ace5", command = lambda : main("mcl_blastp_data")).pack(side = "left")


tk.Button(f1,padx=16,pady=15,bd=10,fg="red2",font=('arial',12,"bold"),\
    width=8,text="Reset **",bg="powder blue",command=reset,activebackground = "#63ace5").pack(side = "left")


tk.Button(f1,padx=16,pady=15,bd=10,fg="red2",font=('arial',12,"bold"),\
    width=8,text="Exit !!",bg="powder blue",command=quit,activebackground = "#fe8a71").pack(side = "left")

## !!____________________________________labels & Entry_______________________________________________

# ! For Threshold_score
tk.Label(window,text = "Threshold Score",font=('arial',16,"bold")).place(x=20,y=200)
tk.Entry(window,bd =5,bg = "#CAF1DE", textvariable = threshold_score_,\
                    font=('arial',12,"bold"), width=16).place(x=300,y=200,height = 35)

# ! For Inflation Factor
tk.Label(window,text = "Inflation factor value",font=('arial',16,"bold")).place(x=20,y=250)
tk.Entry(window,bd =5,bg = "#CAF1DE",textvariable = inflation_factor_,\
                    font=('arial',12,"bold"), width=16).place(x=300,y=250,height = 35)

# ! Clustering File
tk.Label(window,text = "Clustering File Name",font=('arial',16,"bold")).place(x=20,y=300)
tk.Entry(window,bd =5,bg = "#CAF1DE",textvariable = cluster_out,font=('arial',12,"bold"),\
width=16).place(x=300,y=300,height = 35)

# ! For Log Files
tk.Label(window,text = "Log File Name",font=('arial',16,"bold")).place(x=20,y=350)
tk.Entry(window,bd =5,bg = "#CAF1DE", textvariable = log_file_,font=('arial',12,"bold"),\
            width=16).place(x=300,y=350,height = 35)

# ! for CPU _information
tk.Label(window, text = "Cpu to use",font=('arial',16,"bold")).place(x=20,y=400)
tk.Label(window,text = "   <="+str(os.cpu_count()),font=('arial',12,"bold"),fg="red").place(x=200,y=400)
tk.Entry(window,bd =5, textvariable = cpu_count_,bg = "#CAF1DE", font=('arial',12,"bold"),\
            width=16).place(x=300,y=400,height = 35)

# !!__________________________Browse_Button__________________________________________

tk.Label(window,textvariable = species_folder).place(x=100,y=455)
tk.Button(window,text="Species_folder",font = ("arial",12,"bold"),bd =2,
          bg ="#8DE4FD",command =lambda: browse_folder(species_folder)).place(x=20,y=455)


tk.Label(window,textvariable = blastp_folder_).place(x=100,y=490)
tk.Button(window,text="blastp_folder",font = ("arial",12,"bold"),bd =2,
          bg ="#8DE4FD",command =lambda: browse_folder(blastp_folder_)).place(x=20,y=490)

# !  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
tk.Canvas(window,bg = "#C85A17",height=3,width=480).place(x=0,y = 445)
tk.Canvas(window,bg = "#C85A17",height=3,width=480).place(x=0,y = 530)
tk.Canvas(window,bg = "#C85A17",height=3,width=480).place(x=0,y = 612)
tk.Canvas(window,bg = "#C85A17",height=3,width=480).place(x=0,y = 690)

# !lebel and button for subject species
tk.Label(window,textvariable = subject_file_).place(x=100,y=540)
tk.Button(window,text="subject",font=("arial",12,"bold"),bd =2,
          bg ="#8DE4FD",command =lambda: browse_file(subject_file_)).place(x=20,y=540)

# ! lebel and button for query File
tk.Label(window,textvariable = query_file_).place(x=100,y=590)
tk.Button(window,text="query",font = ("arial",12,"bold"),bd =2,
          bg ="#8DE4FD",command =lambda: browse_file(query_file_)).place(x=20,y=575)

## !! _____________________Display input_parameter inside frame 3 with button

# ! Display Result in the top right side
tk.Label(window,text = "All Final Results will publish here ", bg = "gold").place(x=630,y=180)
tk.Label(window,textvariable = result_display,font=('arial',10,"bold"),fg="white" , bg = "#8A0651").place(x=490,y=200)


tk.Label(window,textvariable = var_in,bg = "#98DBC6").place(x=490,y=530) # info_lbl
tk.Button(window,text = "Check Variable", command = display).place(x=490,y=500,bordermode = "outside")
#show_variable

# !!___________________________________Radio_Button_for Matrix_selection_________________________________

tk.Label(window,text = "Choose a Matrix for blastp",font=('arial',14,"bold")).place(x=50,y=630)

tk.Radiobutton(window,text = "BLOSUM45",variable=bp_matrix,value = "BLOSUM45", \
               command = show_choice,font=('arial',10,"bold")).place(x=10,y=660)
tk.Radiobutton(window,text = "BLOSUM62",variable = bp_matrix,value="BLOSUM62",\
               command = show_choice,font=('arial',10,"bold")).place(x=120,y=660)
tk.Radiobutton(window,text = "BLOSUM80",variable = bp_matrix,value="BLOSUM80",command = show_choice,font=('arial',10,"bold")).place(x=240,y=660)

window.mainloop()
