#!/usr/bin/python
import glob
import sys
import subprocess
import pdb
import time
from progress.bar import Bar
import multiprocessing
import Queue
import numpy.matlib
import numpy
import copy
import operator
import pprint
from itertools import *
import argparse
import datetime

print "="*77
print "\t"*2+" "*4+"owPReMark for Linux version 1.0 (C) 2013"
print "\t"*6+" "*4+"JungWon Park, Sunshin Kim"
print "\t"*6+" "*2+"GPLv3 Licence (see LICENCE)"
print "="*77


parser = argparse.ArgumentParser(description = 'The owPReMark is a program to detect ortholog between the protein sequences from different genomes and to cluster orthologs to ortholog groups. ', 
                                 add_help=True, prefix_chars='-+')
essential_group = parser.add_argument_group('essential arguments')
blastp_group = parser.add_argument_group('blastp optional arguments')
mcl_group = parser.add_argument_group('clustering optional arguments')
essential_group.add_argument('-M','--mode', 
                             action='store', 
                             nargs='+',
                             dest='mode',
                             choices=['1','2','3'],
                             help="1 : BLASTP. 2 : BLASTP using precalculated data. 3 : Clustering. ex) 1 or 2 or 1 3 or 2 3")
essential_group.add_argument('-g','--genomes', 
                             nargs='+',
                             dest='genomes', 
                             help="select an genomes to detect orthologs.")
parser.add_argument('-c','--cpu', 
                             action='store', 
                             default='1', 
                             dest='cpu_count', 
                             type=int, 
                             help='set the number of CPU to use in program. The default is "1". This system has '+str(multiprocessing.cpu_count())+" CPUs")
blastp_group.add_argument('+r', 
                             action='store_true', 
                             default=False,
                             dest='save_raw_blastp_score', 
                             help='save the raw score of blastp to the score_file dirctory.')
blastp_group.add_argument('-t', '--threshold', 
                             action='store', 
                             default=0, 
                             type=int,
                             dest='threshold_score', 
                             help='set the threshold score. The threshold score is an allowable range to be an ortholog.'+
                             'If the score difference between backward best hit score and forward best hit score is less than threshold score, it is considered the forward best hit pair to be the ortholog.'+
                             'That is  called the "one-way threshold best hit" by us. The default value is "0". "backward best hit score - forward best hit score <= threshold score"')
blastp_group.add_argument('-m','--matrix', 
                             action='store', 
                             default='BLOSUM62', 
                             dest='blastp_matrix', 
                             choices=['BLOSUM45', 'BLOSUM62', 'BLOSUM82'],
                             help='select the matrix to use in blastp. The default value is "BLOSUM62".')
blastp_group.add_argument('-s','--species', 
                             action='store', 
                             default = './species/', 
                             dest='Species', 
                             help='set the path of species directory. The default is  "./species/".')
blastp_group.add_argument('-b','--blastp', 
                             action='store', 
                             default ='/usr/bin/blastp', 
                             dest='Blastp', 
                             help='set the path of blastp file to run the blastp program. The default is "blastp".')
blastp_group.add_argument('-F','--scorefile', 
                             action='store', 
                             default='./score_file/',
                             dest='Score_file', 
                             help='The default path is "./score_file/".')
blastp_group.add_argument('-B','--blastp_data', 
                             action='store', 
                             default='./blastp_data/', 
                             dest='Blastp_data', 
                             help='set the path of precalculated blastp data. The default path is "./blastp_data/".')
mcl_group.add_argument('-i' ,'--IF', 
                             action='store', 
                             default=1.4, 
                             dest='inflation_factor', 
                             type=float, 
                             help ='The default value is "1.4".')
mcl_group.add_argument('-l', '--loop', 
                             action='store', 
                             default=100, 
                             type=int,
                             dest='infinite_loop', 
                             help='prevent the infinite loop of MCL algorithm. The default value is "60".')
mcl_group.add_argument('-o','--out', 
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

class OneWayBlastp:
    def GetMatrixNumber(self):
        print "\n1. BLOSUM45\n2. BLOSUM62\n3. BLOSUM82\n4. Quit"
        get_matrix_number = raw_input("\nEnter a matrix number: ")
        if get_matrix_number == None:
            matrix_number = "BLOSUM62"
        elif get_matrix_number == "1":
            matrix_number = "BLOSUM45"
        elif get_matrix_number == "2":
            matrix_number = "BLOSUM62"
        elif get_matrix_number == "3":
            matrix_number = "BLOSUM82"
        elif get_matrix_number == "4":
            sys.exit(0)
        else :
            print "Wrong typing! Try again!"
            sys.exit(2)
        print "\nMatrix number is : " + matrix_number + "\n"
        return matrix_number 

    def GetQuerySequence(self, genome):
        gene_sequence=""
        gene_sequence_list= []
        try:
            with open(Species+genome) as gene:
                for each_line in gene:
                    if ">" in each_line:
                        if gene_sequence != "":                        
                            gene_sequence_list.append(gene_sequence)
                            gene_sequence = ""
                    gene_sequence = gene_sequence+each_line
                gene_sequence_list.append(gene_sequence)                   
            return gene_sequence_list                        

        except IOError as err:
            print "IOError occurred in GetQuerySequence function : " + str(err)    

    def WriteQuery(self, query, parallel_num):
        try:
            with open("./query"+"_"+str(parallel_num), "w") as write_query:
                write_query.write(query)
        except IOError as err:
            print "IOError occurred in WriteQuery function : " + str(err)

    def RunBlast(self, subject, parallel_num):
        subject = Species+subject       
        run_blastp =subprocess.Popen([Blastp, "-query", "./query"+"_"+str(parallel_num), "-subject", subject,
                                    "-matrix", blastp_matrix, "-outfmt", "10 qseqid sseqid score length"],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        run_blastp_stream = run_blastp.communicate()
        run_blastp_output_stream = run_blastp_stream[0]
        run_blastp_error_stream = run_blastp_stream[1]          
        return run_blastp_output_stream

    def Get_Same_Species_Forward_Best_Hit(self, blastp_score): 
        """Search the forward best hit among the blastp scores of same species. Because there are an duplicated genes in a same genome.
        When the blastp score compare with blastp score of duplicate gene, if score and length are same, blasp score of duplicated gene is added to a second best score.""" 
        blastp_score_split_list = []                   
        temp_best_score =['-1','-1','-1']
        second_temp_best_score = []
        best_score = [] 
        second_best_score = []       
        blastp_score_split = blastp_score.split("\n") 
        del blastp_score_split[-1] # delete of ['']   ex) ['gi,gi,1,1','gi,gi,2,2','']
        for i in blastp_score_split:       
            blastp_score_element = i.split(',') 
            blastp_score_split_list.append(blastp_score_element)
        for k in blastp_score_split_list: # ex) k is ['gi|15605613|ref|NP_212986.1|', 'gi|15605613|ref|NP_212986.1|', '3702', '699']
            if k[0] == k[1]:
                best_score.append(k)        
            elif k[0] != k[1]:            
                if int(k[2]) > int(temp_best_score[2]) : #Compare score
                    temp_best_score = k
                elif int(k[2]) == int(temp_best_score[2]):
                    if int(k[3]) > int(temp_best_score[3]): #compare length
                        temp_best_score = k
                    elif int(k[3]) == int(temp_best_score[3]):
                        second_temp_best_score.append(k) 
    #    print "############ temp best score ############", temp_best_score
        second_best_score.append(temp_best_score)
        for j in second_temp_best_score:
            if j[2] == temp_best_score[2] and j[3] == temp_best_score[3]:
                second_best_score.append(j)
        for m in second_best_score:
            if (best_score[0][2] == m[2] and int(best_score[0][3]) <= int(m[3])) or int(best_score[0][2]) < int(m[2]): # '104' < '23' is True because of string. So the int function is used.
                best_score.append(m)

        return best_score

    def GetForwardBestHit(self, blastp_score):
        """Search the forward best hit among the blastp scores of same species."""
        blastp_score_split_list = []
        temp_best_score =['-1','-1','-1']
        second_temp_best_score = []
        best_score = []
        blastp_score_split = blastp_score.split("\n") 
        del blastp_score_split[-1] # delete of ['']   ex) ['gi,gi,1,1','gi,gi,2,2','']
        for i in blastp_score_split:       
            blastp_score_element = i.split(',') 
            blastp_score_split_list.append(blastp_score_element)

        for k in blastp_score_split_list: # ex) k is ['gi|15605613|ref|NP_212986.1|', 'gi|15605613|ref|NP_212986.1|', '3702', '699']
    #        print ">>>>>>>>>>>>>>>GetForwardBestHit   k", k
            if int(k[2]) > int(temp_best_score[2]) : #Compare score
                temp_best_score = k
            elif int(k[2]) == int(temp_best_score[2]):
                if int(k[3]) > int(temp_best_score[3]): #compare length
                    temp_best_score = k
                elif int(k[3]) == int(temp_best_score[3]):
                    second_temp_best_score.append(k) 
    #    print "############ temp best score ############", temp_best_score
        best_score.append(temp_best_score)
        for j in second_temp_best_score:
            if j[2] == temp_best_score[2] and j[3] == temp_best_score[3]:
                best_score.append(j)

        return best_score, blastp_score_split_list

    def DivisionParallelQuery(self, queryV, query_division_value, cpu_count, queryV_len):
        parallel_query = []
        parallel_query_start = 0

        if queryV_len % cpu_count == 0 : # perfect division 
            for i in range(cpu_count):
                i += 1
                if parallel_query_start == 0 :
                    parallel_query.append(queryV[parallel_query_start:i*query_division_value])
                    parallel_query_start += 1
                else :
                    parallel_query.append(queryV[parallel_query_start*query_division_value:i*query_division_value])
                    parallel_query_start += 1
        else : #imperfect division
            for i in range(cpu_count):
                i += 1
                if parallel_query_start == 0 :
                    parallel_query.append(queryV[parallel_query_start:i*query_division_value])
                    parallel_query_start += 1
                elif i < cpu_count :
                    parallel_query.append(queryV[parallel_query_start*query_division_value:i*query_division_value])
                    parallel_query_start += 1
                elif i == cpu_count :
                    parallel_query.append(queryV[parallel_query_start*query_division_value:])
                    parallel_query_start += 1       

        return parallel_query

    def RunParallelQuery(self, species_of_query, species_of_subject,queryV, parallel_num):    
        """ Run the following functions. WriteQuery, RunBlast, Get_Same_Species_Forward_Best_Hit, GetForwardBestHit
        Save the files which are oneway_threshold_best_hit, second_oneway_threshold_best_hit, blastp_score_split_list and raw_blastp_score (optional) by each species. """
        bar = Bar('Processing '+str(parallel_num), max = len(queryV)) #progressing bar setting    
        OWBP = OneWayBlastp()       
        for j in queryV:
            bar.next() #progressing bar print
            OWBP.WriteQuery(j,parallel_num)        
            blastp_score  = OWBP.RunBlast(selected_species_dic[species_of_subject], parallel_num)        
            if blastp_score != '': # Check whether blastp_score has the value

                best_score, blastp_score_split_list = OWBP.GetForwardBestHit(blastp_score)

                if species_of_query == species_of_subject: # ex) AAE == AAE. It will save best_score without reversing RunBlast. 
                    same_species_forward_best_score = OWBP.Get_Same_Species_Forward_Best_Hit(blastp_score)
                    for best_score_element in same_species_forward_best_score:
                        if best_score_element[0] == best_score_element[1]:  #  ex) [A1 of AAE, A1 of AAE, 30]                                           

                            with open(Score_file+selected_species_dic[species_of_query]+"_"+selected_species_dic[species_of_subject]+"_oneway_threshold_best_hit_S"+str(threshold_score), "a") as oneway_threshold_best_hit:

                                save_best_score = selected_species_dic[species_of_query]+"_"+best_score_element[0].split("|")[1]+" "+selected_species_dic[species_of_subject]+"_"+best_score_element[1].split("|")[1]+" "+best_score_element[2]+"\n" # best_score_element[0].split("|") ==> ['gi', '15642790', 'ref', 'NP_227831.1', '']
                                oneway_threshold_best_hit.write(save_best_score)

                        else: # ex) [A1 of AAE, A2 of AAE, 30]
                            with open(Score_file+selected_species_dic[species_of_query]+"_"+selected_species_dic[species_of_subject]+"_second_oneway_threshold_best_hit_S"+str(threshold_score), "a") as second_oneway_threshold_best_hit:
                                second_save_best_score = selected_species_dic[species_of_query]+"_"+best_score_element[0].split("|")[1]+" "+selected_species_dic[species_of_subject]+"_"+best_score_element[1].split("|")[1]+" "+best_score_element[2]+"\n"
                                second_oneway_threshold_best_hit.write(second_save_best_score)

                else: # If species_of_query not equal with species_of_subject, run reversing RunBlast
                    for best_score_element in best_score:
                        with open(Score_file+selected_species_dic[species_of_query]+"_"+selected_species_dic[species_of_subject]+"_"+"best_score_S"+str(threshold_score)+"_"+str(parallel_num), "a") as save_best_hit:
                            best_score_save = selected_species_dic[species_of_query]+"_"+best_score_element[0].split("|")[1]+" "+selected_species_dic[species_of_subject]+"_"+best_score_element[1].split("|")[1]+" "+best_score_element[2]+" "+best_score_element[3]+"\n"
                            save_best_hit.write(best_score_save)

                    for blastp_score_split_list_element in blastp_score_split_list:
                        with open(Score_file+selected_species_dic[species_of_query]+"_"+selected_species_dic[species_of_subject]+"_"+"blastp_score_split_list_S"+str(threshold_score)+"_"+str(parallel_num), "a") as save_blastp_score_split_list:                            
                            blastp_score_split_list_save = selected_species_dic[species_of_query]+"_"+blastp_score_split_list_element[0].split("|")[1]+" "+selected_species_dic[species_of_subject]+"_"+blastp_score_split_list_element[1].split("|")[1]+" "+blastp_score_split_list_element[2]+" "+blastp_score_split_list_element[3]+"\n"
                            save_blastp_score_split_list.write(blastp_score_split_list_save)

            if save_raw_blastp_score :
                 with open(Score_file+selected_species_dic[species_of_query]+"_"+selected_species_dic[species_of_subject]+"_S"+str(threshold_score)+"_"+str(parallel_num), "a") as save_blastp:
                     save_blastp.write(blastp_score)
        bar.finish() # progressing bar finish 
        return 



    def Search_Equal_BBH_Data(self, target_A):
        """Search the equal backward best hit data. ex) AAE_AAE_backward_best_hit """        
        put_data = equal_BBH_data_dic[target_A]
        if put_data[1] == 0:
            pass
        else :
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
                        tasks.put(copy_second_put_data)  # Don't put results as queue. Because the tasks will put copy_second_put_data to results as queue.
                        i[2] = 0
    #                    print "---put zero in second_equal_BBH_data---", i
        return

    def Search_Unequal_BBH_Data(self, target_B):
        """Search the unequal backward best hit data. ex) AAE_CAC_backward_best_hit"""
        for i in unequal_BBH_data:
            if i[2] == 0 :
                pass
            else:
                if target_B[0] == i[0] or target_B[0] == i[1] or target_B[1] == i[0] or target_B[1] == i[1]:
                    copy_i = copy.copy(i)
                    tasks.put(copy_i)
                    unequal_BBH_data[unequal_BBH_data.index(i)][2] = 0
        return

    def Matching_BBH(self, target):
        """ Match the backward best hit """ 
        OWBP = OneWayBlastp()
        if target[2] == 0 :
            return 

        else:
            copy_target = copy.copy(target)
            OWBP.Search_Equal_BBH_Data(copy_target[0])
            OWBP.Search_Equal_BBH_Data(copy_target[1])
            results.put(copy_target)
            unequal_BBH_data[unequal_BBH_data.index(target)][2] = 0

        for j in unequal_BBH_data:
            if j[2] == 0 :
                pass

            else:
                if copy_target[0]==j[0] or copy_target[0]==j[1] or copy_target[1]==j[0] or copy_target[1]==j[1]:
                    copy_j = copy.copy(j)
    #                print "targ_get , j = %s %s" % (copy_target, j)
                    unequal_BBH_data[unequal_BBH_data.index(j)][2] = 0
                    tasks.put(copy_j)

        while not tasks.empty():              
            get_task = tasks.get()
            OWBP.Search_Equal_BBH_Data(get_task[0])
            OWBP.Search_Equal_BBH_Data(get_task[1])
            results.put(get_task)
            OWBP.Search_Unequal_BBH_Data(get_task)

            
    def Read_Equal_BBH(self, path):            
        with open(path+"_oneway_threshold_best_hit_S"+str(threshold_score), 'r') as equal_RBH:                          
            for j in equal_RBH:                    
                split_data = j.split()
                split_data[2] = int(split_data[2])
                equal_BBH_data.append(split_data)                   
                equal_BBH_data_dic[split_data[0]] = split_data[1:]
        try :
            with open(path+"_second_oneway_threshold_best_hit_S"+str(threshold_score), 'r') as second_equal_RBH:                
                for j in second_equal_RBH:
                    split_data = j.split()
                    split_data[2] = int(split_data[2])                        
                    second_equal_BBH_data.append(split_data)       
        except : 
            pass

    def Read_Unequal_BBH(self, path):
        with open(path+"_oneway_threshold_best_hit_S"+str(threshold_score), 'r') as unequal_RBH:
            for j in unequal_RBH:
                split_data = j.split()
                split_data[2] = int(split_data[2])        
                unequal_BBH_data.append(split_data)
        

def Oneway_Threshold_Best_Hit(mode):
    OWBP = OneWayBlastp()    
    process_list = []
    backward_best_hit_work_list = []
    if "1" in mode:
        for i in user_selected_number: #Select species to write query    
            queryV = OWBP.GetQuerySequence(selected_species_dic[i])
            for k in user_selected_number: #Select of subject
                if k < i: # gene ====> query 1->1 1->2 1->3 2->2 2->3  
                    continue

                else :
                    print "Doing the blastp & forward best hit searches between %s genome and %s genome" % (selected_species_dic[i], selected_species_dic[k])

                    queryV_len = len(queryV)
                    if cpu_count == 1:
                        blastp_time_start = time.time()
                        OWBP.RunParallelQuery(i, k, queryV, cpu_count)
                        blastp_time_end = time.time()
                        print "The blastp & forward best hit searches took %.2f minutes" % ((blastp_time_end-blastp_time_start)/60)
                    else :
                        if queryV_len < cpu_count: #If the number of queryV_len is less than cpu_count, Remark will select the number of queryV_len.
                            blastp_time_start = time.time()
                            parallel_query = OWBP.DivisionParallelQuery(queryV, 1, queryV_len, queryV_len) # 1 is query_division_value. Because queryV_len / queryV_len(=cpu_count) is 1.
                            for m in range(queryV_len):
                                process = multiprocessing.Process(target=OWBP.RunParallelQuery, args=(i, k, parallel_query[m], m+1))
                                # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                process_list.append(process)
                                process.start()              
                            for n in process_list:
                                n.join()
                            blastp_time_end = time.time()
                            print "The blastp & forward best hit searches took %.2f minutes" % ((blastp_time_end-blastp_time_start)/60)
                        else :
                            blastp_time_start = time.time()
                            query_division_value = queryV_len / cpu_count
                            parallel_query = OWBP.DivisionParallelQuery(queryV, query_division_value, cpu_count, queryV_len)
                            for m in range(cpu_count):
                                process = multiprocessing.Process(target=OWBP.RunParallelQuery, args=(i, k, parallel_query[m], m+1))
                                # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                process_list.append(process)
                                process.start()              
                            for n in process_list:
                                n.join()
                            blastp_time_end = time.time()
                            print "The blastp & forward best hit searches took %.2f minutes" % ((blastp_time_end-blastp_time_start)/60)
                    if not i == k :                         
                        backward_best_hit_work_list.append((i,k, queryV_len))
    elif "2" in mode:
        for i in user_selected_number: #Select species to write query    
            queryV = OWBP.GetQuerySequence(selected_species_dic[i])
            for k in user_selected_number: #Select of subject                               
                if Blastp_data+selected_species_dic[i]+"_"+selected_species_dic[k]+"_oneway_threshold_best_hit_S"+str(threshold_score) in precalculated_data_list :
                    used_precalculated_data_list.append(selected_species_dic[i]+"_"+selected_species_dic[k])
                    continue
                else :                              
                    if k < i: # gene ====> query 1->1 1->2 1->3 2->2 2->3  
                        continue

                    else :
                        print "Doing the blastp & forward best hit searches between %s genome and %s genome" % (selected_species_dic[i], selected_species_dic[k])

                        queryV_len = len(queryV)
                        if cpu_count == 1 :
                            blastp_time_start = time.time()
                            OWBP.RunParallelQuery(i, k, queryV, cpu_count)
                            blastp_time_end = time.time()
                            print "The blastp & forward best hit searches took %.2f minutes" % ((blastp_time_end-blastp_time_start)/60)
                        else :
                            if queryV_len < cpu_count: #If the number of queryV_len is less than cpu_count, Remark will select the number of queryV_len.
                                blastp_time_start = time.time()
                                parallel_query = OWBP.DivisionParallelQuery(queryV, 1, queryV_len, queryV_len) # 1 is query_division_value. Because queryV_len / queryV_len(=cpu_count) is 1.
                                for m in range(queryV_len):
                                    process = multiprocessing.Process(target=OWBP.RunParallelQuery, args=(i, k, parallel_query[m], m+1))
                                    # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                    process_list.append(process)
                                    process.start()              
                                for n in process_list:
                                    n.join()
                                blastp_time_end = time.time()
                                print "The blastp & forward best hit searches took %.2f minutes" % ((blastp_time_end-blastp_time_start)/60)
                            else :
                                blastp_time_start = time.time()
                                query_division_value = queryV_len / cpu_count
                                parallel_query = OWBP.DivisionParallelQuery(queryV, query_division_value, cpu_count, queryV_len)
                                for m in range(cpu_count):
                                    process = multiprocessing.Process(target=OWBP.RunParallelQuery, args=(i, k, parallel_query[m], m+1))
                                    # args( i => species of query , k => species of subject, m+1 => cpu_count ex) 1, 2 ...)
                                    process_list.append(process)
                                    process.start()              
                                for n in process_list:
                                    n.join()
                                blastp_time_end = time.time()
                                print "The blastp & forward best hit searches took %.2f minutes" % ((blastp_time_end-blastp_time_start)/60)
                        new_calculated_data_list.append(selected_species_dic[i]+"_"+selected_species_dic[k])
                        if not i == k :                               
                            backward_best_hit_work_list.append((i,k, queryV_len))                                                          
    return backward_best_hit_work_list

def Backward_Best_Hit(args):
    species_of_query, species_of_subject, queryV_len = args    
    start_time_BBH = time.time()
    forward_best_hit_score_list = []
    blastp_score_split_list = []    
    print "Start to run the backward best hit between %s genome %s genome" % (selected_species_dic[species_of_query], selected_species_dic[species_of_subject])
    if queryV_len < cpu_count : #If the number of queryV_len is less than cpu_count, the cpu_count is changed to queryV_len.
        for parallel_num in range(queryV_len):
            parallel_num += 1                                  
            with open(Score_file+selected_species_dic[species_of_query]+"_"+selected_species_dic[species_of_subject]+"_"+"best_score_S"+str(threshold_score)+"_"+str(parallel_num), "r") as best_hit_score:
                for each_line in best_hit_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
        #            print split_each_line
                    forward_best_hit_score_list.append(split_each_line)
            with open(Score_file+selected_species_dic[species_of_query]+"_"+selected_species_dic[species_of_subject]+"_"+"blastp_score_split_list_S"+str(threshold_score)+"_"+str(parallel_num), 'r') as blastp_score :
                for each_line in blastp_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
        #            print split_each_line
                    blastp_score_split_list.append(split_each_line)   
    else :            
        for parallel_num in range(cpu_count):
            parallel_num += 1
            with open(Score_file+selected_species_dic[species_of_query]+"_"+selected_species_dic[species_of_subject]+"_"+"best_score_S"+str(threshold_score)+"_"+str(parallel_num), "r") as best_hit_score:
                for each_line in best_hit_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
        #            print split_each_line
                    forward_best_hit_score_list.append(split_each_line)
            with open(Score_file+selected_species_dic[species_of_query]+"_"+selected_species_dic[species_of_subject]+"_"+"blastp_score_split_list_S"+str(threshold_score)+"_"+str(parallel_num), 'r') as blastp_score :
                for each_line in blastp_score:
                    split_each_line = each_line.split(" ")
                    split_each_line[2] = int(split_each_line[2])
                    split_each_line[3] = int(split_each_line[3])
        #            print split_each_line
                    blastp_score_split_list.append(split_each_line)    

    bar = Bar("Searching : "+selected_species_dic[species_of_query]+"-"+selected_species_dic[species_of_subject], max = len(forward_best_hit_score_list))

    for forward_best_hit_score_element in forward_best_hit_score_list:
        matching_list = []        
        backward_best_score = ['-1','-1','-1']
        bar.next()
        for element in blastp_score_split_list:
            if element[1] == forward_best_hit_score_element[1]:                
                matching_list.append(element) 

        for element in matching_list:
            if int(element[2]) > int(backward_best_score[2]):
                backward_best_score = element
#        with open('./'+selected_species_dic[species_of_query]+"_"+selected_species_dic[species_of_subject]+'_subtraction'+"_"+str(threshold_score), 'a') as subtraction :
#            save_data = int(backward_best_score[2]) - int(forward_best_hit_score_element[2])
#            subtraction.write(str(save_data)+"\n")

        if int(backward_best_score[2]) - int(forward_best_hit_score_element[2]) <= threshold_score :           
            with open(Score_file+selected_species_dic[species_of_query]+"_"+selected_species_dic[species_of_subject]+"_oneway_threshold_best_hit_S"+str(threshold_score), "a") as other_oneway_threshold_best_hit:
                save_data = forward_best_hit_score_element[0]+" "+forward_best_hit_score_element[1]+" " +str(int(forward_best_hit_score_element[2]))+"\n"#           
                other_oneway_threshold_best_hit.write(save_data)

    bar.finish()
    finish_time_BBH = time.time()
    RBH_time = float((finish_time_BBH - start_time_BBH)/60)
    print  "BackwardBestHit of %s-%s took %.2f minutes" % (selected_species_dic[species_of_query], selected_species_dic[species_of_subject], RBH_time )
    return RBH_time





class MclDetect:
    def Generating_Matrix_Clustering_Ortholog(self, element_set, bar):  
        """ Generate the matrix of clustering ortholog. """

        MD = MclDetect()

        row_data = []
        col_data = []
        temp_results = Queue.Queue()  
        bar.next()
        for element in element_set:    
            if row_data.count(element[0]) > 0 : # if element[0] exist, returning the index in the row_data.
                row = row_data.index(element[0]) # element[0] is data of row. ['gi|15606057|ref|NP_213434.1|', 'gi|15606057|ref|NP_213434.1|', '3823\n']

            else :
                row = len(row_data) 
                row_data.append(element[0])
                if col_data.count(element[0]) < 1: # if element[0] doesn't exist, appending the element[0] to the col_data.
                    col_data.append(element[0])                    

            if col_data.count(element[1]) > 0:
                col = col_data.index(element[1]) #element[1] is data of col.

            else:
                col = len(col_data) 
                col_data.append(element[1])    
                if row_data.count(element[1]) < 1:
                    col_data.append(element[1])

            temp_results.put([row,col,element[2]])            

        score_matrix = numpy.matlib.zeros((len(row_data),len(col_data)), dtype=numpy.float) # create a new matrix of given shape(the size_resuls) and type, filled with zeros.              
        while not temp_results.empty():
            get_temp_results = temp_results.get()
            row = get_temp_results[0]
            col = get_temp_results[1]
            score_matrix[row,col] = get_temp_results[2]
            score_matrix[col,row] = get_temp_results[2]

        if len(row_data)*len(col_data) > 4 :  # If the elements of row and col is less than 2, it is excluded.
            if len(row_data) > 1000 and cpu_count > 1 : #The big size of matrix(bigger than 1000 X 1000) will be computed by Parallel_Matrix_Multiplication_Using_Numpy function.            
                score_matrix = MD.Parallel_MCL(score_matrix)
            else :            
                score_matrix = MD.MCL(score_matrix)
            MD.Clustering(row_data, col_data, score_matrix)    

    def Parallel_MCL(self, score_matrix):    
        count = 0
        infinitesimal_value = 10**-10
        idempotent_value = 1   

        while idempotent_value > infinitesimal_value: # > infinitesimal_value 
            MCL_time_start = time.time()        
            pool = multiprocessing.Pool(cpu_count) # create a expansion_matrix  
            multiplication_results = pool.map(Parallel_Matrix_Multiplication_Using_Numpy, izip(score_matrix,repeat(score_matrix)))
            pool.close()
            pool.join()

            pool = multiprocessing.Pool(cpu_count) # create a inflation_matrix(part 1)   
            power_results = pool.map(Parallel_Matrix_Power_Using_Numpy, multiplication_results)
            pool.close()
            pool.join()

            sum_matrix = 0
            for i in power_results: 
                sum_matrix = i + sum_matrix

            pool = multiprocessing.Pool(cpu_count) # create a inflation_matrix(part 2)   
            divide_results = pool.map(Parallel_Matrix_Divide_Using_Numpy, izip(power_results, repeat(sum_matrix)))
            pool.close()
            pool.join()

            score_matrix = numpy.concatenate(divide_results)                                

            idempotent_value = abs(numpy.sum(score_matrix) - numpy.sum(multiplication_results)) # identify whether inflation_matrix is idempotent matrix or not.

            count += 1
            if count > infinite_loop : # It will prevent the infinite loop of MCL algorithm.       
                break     
            MCL_time_finish = time.time()
            if verbose :
                print " MCL time : %f, count : %d, matrix size : %d * %d" % ((MCL_time_finish - MCL_time_start)/60, count, score_matrix[0].size, score_matrix[0].size)                     
        return score_matrix    

    def MCL(self, score_matrix):    
        count = 0
        infinitesimal_value = 10**-10
        idempotent_matrix = numpy.matlib.ones((2,2))    
        while idempotent_matrix.sum() > infinitesimal_value: # > infinitesimal_value 
            MCL_time_start = time.time()        
            expansion_matrix = score_matrix ** 2      
            score_matrix = numpy.power(expansion_matrix, inflation_factor)
            score_matrix_sum =  score_matrix.sum(axis = 0)
            score_matrix = numpy.divide(score_matrix, score_matrix_sum) # create a inflation_matrix        
            idempotent_matrix =abs(score_matrix - expansion_matrix) # identify whether inflation_matrix is idempotent matrix or not.        
            count += 1
            if count > infinite_loop : # It will prevent the infinite loop of MCL algorithm.       
                break     
            MCL_time_finish = time.time()
            if verbose :
                print " MCL time : %f, count : %d, matrix size : %d * %d" % ((MCL_time_finish - MCL_time_start)/60, count, score_matrix[0].size, score_matrix[0].size)            
        return score_matrix

    def Clustering(self, row_data, col_data, score_matrix):
        global cluster_count
        global ortholog_count
        ortholog_temp_list = []            
        for i in range(len(row_data)):
            ortholog_list = []
            ortholog_sum = 0
            ortholog = Queue.Queue() # It is Queue which is put the ortholog.
            gene_id_queue = Queue.Queue() # It is Queue which is put the ortholog having changed gene ID 
            for j in range(len(col_data)):
                if 0.1 <= score_matrix[i,j]:            
                    ortholog.put(col_data[j])
                    ortholog_list.append(col_data[j])

            if ortholog.qsize() >= 3: # If the ortholog queue  has the element of ortholog more than 2, it will be printed.
                for element in ortholog_list:
                    try :
                        ortholog_sum += ortholog_temp_list.index(element)+1

                    except ValueError:                                                                     
                        with open(Cluster_out+"_geneID_S"+str(threshold_score)+"_"+str(inflation_factor), "a") as ortholog_list_save:
                            ortholog_print = "cluster "+str(cluster_count)+" :" 
                            ortholog_list_save.write(ortholog_print)                                                           

                            while not ortholog.empty():
                                get_ortholog = ortholog.get()
                                ortholog_list_save.write("\t"+get_ortholog)

                                try:
                                    get_ortholog_split = get_ortholog.split('_') # get_ortholog --> ECO_170082288,  get_ortholog.split('_') --> ['ECO', '170082288']
                                    gene_id_queue.put(get_ortholog_split[0]+"_"+gene_id_dic[get_ortholog_split[1]])                                
                                except KeyError:
                                    gene_id_queue.put(get_ortholog)  # If the gene_id_dic don't have get_ortholog, it will print the original ID(get_ortholog).
                            ortholog_list_save.write("\n")
                        with open(Cluster_out+"_KO_ID_S"+str(threshold_score)+"_"+str(inflation_factor), "a") as ortholog_list_geneID_save:                                 
                            ortholog_list_geneID_print = "cluster "+str(cluster_count)+" :" 
                            ortholog_list_geneID_save.write(ortholog_list_geneID_print)
                            while not gene_id_queue.empty():
                                ortholog_list_geneID_save.write("\t"+gene_id_queue.get())
                                ortholog_count += 1                                                   
                            cluster_count += 1
                            ortholog_list_geneID_save.write("\n")
                        break
                ortholog_temp_list = operator.concat(ortholog_temp_list, ortholog_list) 



    def Parallel_Matrix_Multiplication_Using_Numpy(self, data):
        matrix_element , matrix = data     
        result = matrix_element * matrix   
        return result

    def Parallel_Matrix_Power_Using_Numpy(self, matrix_element):
        power_matrix_element = numpy.power(matrix_element, inflation_factor)
        return power_matrix_element

    def Parallel_Matrix_Divide_Using_Numpy(self, data):
        matrix_element, sum_data = data
        return numpy.divide(matrix_element, sum_data)

def Read_Species_List(pr=0):
    """ If pr is 1, it will print "Species_List". """ 
    read_species =glob.glob(command_options.Species+"*")
    global selected_species_dic
    global backward_selected_species_dic   
    selected_species_dic = {}
    backward_selected_species_dic = {}
    for i, species in enumerate(sorted(read_species), start=1):
        selected_species_dic[i] = species.split('/')[-1]
        backward_selected_species_dic[species.split('/')[-1]] = i 
        if pr == 1 :
            print str(i)+".", species.split('/')[-1]
        number = i
    return selected_species_dic, backward_selected_species_dic, number

def Del_File(path, file):
    del_file =subprocess.Popen(["rm "+path+file], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    del_file_stream = del_file.communicate()   
    if not del_file_stream[1]:
        print "Done to del "+path+file
    elif del_file_stream[1]:
        print del_file_stream[1]

def Check_File(file):
    file_list = glob.glob(file+'*')    
    if (Cluster_out+"_geneID_S"+str(threshold_score)+"_"+str(inflation_factor) or Cluster_out+"_KO_ID_S"+str(threshold_score)+"_"+str(inflation_factor)) in file_list:
        print "Please, set other name of output."
        sys.exit(2)
        

print "Blastp = ", command_options.Blastp
print "Blastp_data = ", command_options.Blastp_data
print "blstp_matrix = ", command_options.blastp_matrix
print "cpu_count = ", command_options.cpu_count
print "genomes", command_options.genomes
print "infinite_loop = ", command_options.infinite_loop
print "inflation_factor = ", command_options.inflation_factor
print "mode : ", command_options.mode
print "Cluster_out : ", command_options.Cluster_out
print "threshold_score : ", command_options.threshold_score
print "save_raw_blastp_score : ", command_options.save_raw_blastp_score
print "Score_file : ", command_options.Score_file
print "Species : ", command_options.Species
print "verbose : ", command_options.verbose
print 

    
def main(argv):

    global Cluster_out, threshold_score, inflation_factor, user_selected_number, Species, Blastp_data
    global blastp_matrix, Blastp, precalculated_data_list, used_precalculated_data_list, equal_BBH_data
    global equal_BBH_data_dic, unequal_BBH_data, results, second_equal_BBH_data, tasks, infinite_loop
    global verbose, cluster_count, gene_id_dic, ortholog_count, cpu_count, Score_file, save_raw_blastp_score
    global threshold_score, mode

    OWBP = OneWayBlastp()

    if not sys.argv[1:]:    
        print "1. BLASTP. \n2. BLASTP using precalculated data. \n3. Clustering.\n"
        mode = raw_input(">> Select a mode or modes (1 or 2 or 1 3 or 2 3): ")    
        selected_species_dic, backward_selected_species_dic, number_i = Read_Species_List(pr=1)

        selected_number= raw_input(">> Select an genomes to detect orthologs(e.g. 1 2 3 4 5 or 1-5) : ")

        if selected_number.find('-') > 0:
            SN=selected_number.split("-")
            if int(SN[-1]) > number_i:
                print "\nWrong typing! Try again!\n"
                sys.exit(2)
            else :
                user_selected_number=range(int(SN[0]),int(SN[-1])+1)
                for j in user_selected_number:
                    print selected_species_dic[j],
                print "are selected!!"

        else :        
            user_selected_number = sorted(set([int(read_species) for read_species in selected_number.split()]))                      
            if int(user_selected_number[-1]) > number_i:
                print "\nWrong typing! Try again!\n"
                sys.exit(2)
            else :
                for j in user_selected_number:
                    print selected_species_dic[j],
                print "are selected!!"       

        blastp_matrix = OWBP.GetMatrixNumber()

        cpu_count = input("The processors of CPU are detected. You can use %s processors.\nIf you input >= 2, The Remark will run a parallel computation for the blastp.\n" % multiprocessing.cpu_count()
                              + "Enter the number of process to use in this program (1 ~ %s): " % multiprocessing.cpu_count())
        if "3" in mode :
            inflation_factor = input("Enter the inflation factor to cluster: ")
            Cluster_out = raw_input("Set the name of clustering output : ")



    elif sys.argv[1:] :
        genomes = command_options.genomes        
        mode = command_options.mode
        cpu_count = command_options.cpu_count
        blastp_matrix = command_options.blastp_matrix   
        inflation_factor = command_options.inflation_factor
        selected_species_dic, backward_selected_species_dic, number_i = Read_Species_List()
        user_selected_number = [backward_selected_species_dic[ele] for ele in genomes]
        Cluster_out = command_options.Cluster_out   


    Species = command_options.Species
    Blastp = command_options.Blastp
    Score_file = command_options.Score_file
    Blastp_data = command_options.Blastp_data
    save_raw_blastp_score = command_options.save_raw_blastp_score
    threshold_score = command_options.threshold_score
    verbose = command_options.verbose
    infinite_loop = command_options.infinite_loop

    if "3" in mode :
         Check_File(Cluster_out)

    Del_File(Score_file, "*")

    if "3" in mode :
        Log_file_name = Cluster_out+"_S"+str(threshold_score)+"_"+str(inflation_factor)+".log"
    elif not "3" in mode :
        Log_file_name = 'Log'   

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
   

    if "1" in mode :            
        backward_best_hit_work_list = Oneway_Threshold_Best_Hit(mode)    
        pool = multiprocessing.Pool(cpu_count)    
        results = pool.map(Backward_Best_Hit,backward_best_hit_work_list)
        pool.close()
        pool.join()   

    elif "2" in mode:
        pdb.set_trace()
        used_precalculated_data_list=[]
        new_calculated_data_list=[]
        precalculated_data_list = glob.glob(Blastp_data+"*oneway_threshold_best_hit_S"+str(threshold_score))
        print precalculated_data_list
        backward_best_hit_work_list = Oneway_Threshold_Best_Hit(mode)
        if not backward_best_hit_work_list == []: # If backward_best_hit_work_list is an empty list, pool instance can't finsh the work.                
            pool = multiprocessing.Pool(cpu_count)
            results = pool.map(Backward_Best_Hit,backward_best_hit_work_list)
            pool.close()
            pool.join()
        else : results = [0,0]

    Del_File("./", "query*")
    finish_time_OBH = time.time()
    blastp_time_log = float(((finish_time_OBH - start_time_OBH)/60))
    print  "BLASTP searches + forward best Hit + backwardbest hit took %f minutes" % blastp_time_log
    with open(Log_file_name, 'a') as log:
        log.write("Backward_Best_Hit took "+str(max(results))+" minutes\n")
        log.write("BLASTP + Best_Hit + backward_best_hit searches took "+str(blastp_time_log)+" minutes\n")

    if "3" in mode :
        start_time_clustering = time.time()
        ##########################################################################################
        #generate matrix and calculate the matrix using MCL algorithm and cluster the ortholog."""    
        print "\n>>>> Start MCL algorithm and Clustering ortholog <<<<"    
        equal_BBH_data = []
        unequal_BBH_data = []
        equal_BBH_data_dic ={}
        second_equal_BBH_data = []
        results = Queue.Queue()
        tasks = Queue.Queue()
        cluster_count = 1
        ortholog_count = 0
        gene_id_dic = {}    

        with open("myva=gb", "r") as id_read:
            for i in id_read:
                gene_name, gene_id=i.split()
                gene_id_dic[gene_id.replace("\n","")] = gene_name # remove "\n"

        if "1" in mode :
            for i in user_selected_number:
                for k in user_selected_number:       
                    if k < i :
                        pass
                    elif i == k :
                        OWBP.Read_Equal_BBH(Score_file+selected_species_dic[i]+"_"+selected_species_dic[k])
                    elif i != k :                              
                        OWBP.Read_Unequal_BBH(Score_file+selected_species_dic[i]+"_"+selected_species_dic[k])

        elif "2" in mode :   
            for used_data in used_precalculated_data_list :
                first, second = used_data.split("_")
                if first == second :
                    OWBP.Read_Equal_BBH(Blastp_data+used_data)                           
                elif first != second :
                    OWBP.Read_Unequal_BBH(Blastp_data+used_data)   

            for new_data in new_calculated_data_list :
                first, second = new_data.split("_")
                if first == second :
                    OWBP.Read_Equal_BBH(Score_file+new_data)
                elif first != second :
                    OWBP.Read_Unequal_BBH(Score_file+new_data)          

        matched_BBH_data = []
        matched_BBH_element_data_set = []

        for unequal_RBH_element in unequal_BBH_data:   
            OWBP.Matching_BBH(unequal_RBH_element)  
            temp_results_list = []
            if results._qsize()  != 0: # return the number of results as Queue.      
                while not results.empty():
                    get_results = results.get()
                    temp_results_list.append(get_results)
                matched_BBH_data.append(temp_results_list)

        MD = MclDetect()

        bar = Bar("processing ", max = len(matched_BBH_data))
        for data in matched_BBH_data :
            MD.Generating_Matrix_Clustering_Ortholog(data, bar) 
        bar.finish()
        finish_time_clustering = time.time()
        mcl_time_log = float((finish_time_clustering - start_time_clustering)/60)
        remark_time_log = float((finish_time_clustering - start_time_OBH)/60)
        print "MCL algorithm and Ortholog Clustering took %.2f minutes" % mcl_time_log
        print "owPRemark program took %.2f minutes" % remark_time_log

        if "3" in mode :
            with open(Log_file_name, 'a') as log:
                log.write("Ortholog count : "+str(ortholog_count)+","+" Cluster count : "+str(cluster_count-1)+"\n")    
                log.write("MCL algorithm and Ortholog Clustering took "+str(mcl_time_log)+" minutes\n")
                log.write("Remark program took "+str(remark_time_log)+ "minutes\n")



if __name__ == "__main__":
   main(sys.argv[1:])
