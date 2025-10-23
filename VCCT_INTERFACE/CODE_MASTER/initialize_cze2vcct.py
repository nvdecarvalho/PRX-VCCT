# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 13:58:24 2022

@author: ndecarva
"""
import numpy as np
import re

def read_elements(filename):
    #reads the real elements from _VCCT part
   
    all_els_dict = {}
    read_elements = False
    part_VCCT = False
    with open(filename) as f:
        for line in f:
            line = line.strip() #added to remove \r\n
            if line.find('_VCCT') != -1:
                part_VCCT = True       
            if part_VCCT == True:                                                 
                if line[0:8].lower() == '*element' and line.lower().find('u308') != -1:              
                    read_elements = True
                else:
                    if read_elements == True:          
                        if (line[0] != '*'):
                                s = np.fromstring(line, dtype=int, count=-1, sep=',')
                                all_els_dict[s[0]] = {'nodes':s[1:5]}
                        else: 
                            break
                            
    return all_els_dict

def read_nodes(filename,nodes_vcct):
    #reads the real nodes from _VCCT part
    all_nodes_dict = {}
    read_nodes = False
    part_VCCT = False
    with open(filename) as f:
        for line in f:
            line = line.strip() #added to remove \r\n
            if line.find('_VCCT') != -1:
                part_VCCT = True       
            if part_VCCT == True:                                                 
                if line[0:5].lower() == '*node':
                    read_nodes = True
                else: 
                    if read_nodes == True:
                        if (line[0] != '*'):
                            s = np.fromstring(line, dtype=float, count=-1, sep=',')
                            node = int(s[0])
                            if node in nodes_vcct:
                                all_nodes_dict[int(s[0])] = {'coords':s[1:],'d':0}
                        else:
                            break
               
    return all_nodes_dict


def read_mat_list_elset(inpfile):

    inputfile= open(inpfile,'r')
    params=[[]]
    uel_mat_prop=[[]]
    uel_el_set_name=[[]]
    total_ele_set_name=[[]]
    total_ele_set=[[]]
    uel_fatig_prop_name=[[]]
    uel_fatig_prop_list=[[]]
    uel_fatig_rcurve_name=[[]]
    uel_fatig_rcurve_list=[[]]
#
    templist=[]
    node3=[]
    node31=[]
    node32=[]
    node33=[]
    node34=[]
    node4=[]
    node41=[]
    node42=[]
    node43=[]
    node44=[]
    write=10
    counter=99
    bcounter=0
    lines = [x.strip() for x in inputfile.readlines() if x.strip()]
    for i in range(0,len(lines)):
        line=lines[i]
        if i < len(lines)-2:
            linenext=lines[i+1]
        else:
            linenext='EOF'
        sentence1p=re.split('=|\,|\s',line)
        sentence2p=[x for x in sentence1p if x]
        sentence=sentence2p
        sentencenext=linenext.split()
        nextline0=sentencenext[0][0]
        if  len(sentencenext[0])<=1:
            sentencenext[0]=sentencenext[0]+'*'
            nextline1=sentencenext[0][1]
        else:
            nextline1=sentencenext[0][1]
        if len(sentence)!=0:
            if sentence[0].lower()=='*parameter':
                write=0
            if sentence[0].upper()=='*UEL':
                write=0
            if sentence[0].lower()=='*elset':
                write=0
            if sentence[0].lower()=='**vcct_dadn':
                write=0
            if sentence[0].lower()=='**vcct_r_curve':
                write=0
            if write==0 and len(sentence)!=0:   
                if sentence[0].upper()=='*UEL':
                    write=1
                if sentence[0].lower()=='*parameter':
                    write=1
                if sentence[0].lower()=='*elset':
                    if sentence[2][0:4].lower() == 'vcct':
                        write=1
                if sentence[0].lower()=='**vcct_dadn':
                    write=1
                if sentence[0].lower()=='**vcct_r_curve':
                    write=1
                if sentence[0].lower()=='*parameter' and bcounter==0:
                    counter=0
                    bcounter=bcounter+1
                if sentence[0].upper()=='*UEL':
                    counter=1
                    word2s=sentence[3]
                    word3s=word2s.split(",")
                    uel_el_set_name[counter-1].append(word3s)
                if sentence[0].lower()=='*elset' :
                    if sentence[2][0:4].lower()=='vcct':
                        counter=2
                        word4s=sentence[2]
                        word5s=word4s.split(",")
                        total_ele_set_name[counter-2].append(word5s)
                if sentence[0].lower()=='**vcct_dadn':
                    counter=3
                    word6s=sentence[2]
                    word7s=word6s.split(",")
                    uel_fatig_prop_name[counter-3].append(word7s)
                if sentence[0].lower()=='**vcct_r_curve':
                    counter=4
                    word8s=sentence[2]
                    word9s=word8s.split(",")
                    uel_fatig_rcurve_name[counter-4].append(word8s)
#   
            elif write==1 and len(sentence)!=0: 
                if sentence[0].lower()!='*Step':
                    if counter==0 :
                        if sentence[0][0]!='*':
                            for word in sentence:
                                word1=word.split(",")
                                node2 = [e for e in word1 if e]
                                node32=node32+node2
                                node42=node32
                                params[-1].append(node42)
                                node32=[]
                                node42=[]
                        if nextline0=='*' and nextline1.isalpha==True:
                            write=2
                            node32=[]
                            node42=[]
                            counter=99
                    if counter==1:
                        if sentence[0][0]!='*':
                            for word in sentence:
                                word1=word.split(",")
                                node2 = [e for e in word1 if e]
                                node3=node3+node2
                                node4=node3
                        if sentencenext[0][0]=='*':
                            write=2
                            uel_mat_prop[-1].append(node4)
                            node3=[]
                            node4=[]
                            counter=99
                    if counter==2 :
                        if sentence[0][0]!='*':
                            for word in sentence:
                                word1=word.split(",")
                                node2 = [e for e in word1 if e]
                                node31=node31+node2
                                node41=node31
                        if sentencenext[0][0]=='*':
                            write=2
                            total_ele_set[-1].append(node41)
                            node31=[]
                            node41=[]
                            counter=99
                    if counter==3 :
                        if sentence[1][0]!='*':
                            for word in sentence:
                                word1=word.split("**")
                                word1=filter(None, word1)
                                node2 = [e for e in word1 if e != '**']
                                node33=node33+node2
                                node43=node33
                        if len(sentencenext[0])==2:
                            write=2
                            uel_fatig_prop_list[-1].append(node43)
                            node33=[]
                            node43=[]
                            counter=99
                    if counter==4 :
                        if sentence[1][0]!='*':
                            for word in sentence:
                                word1=word.split("**")
                                word1=filter(None, word1)
                                node2 = [e for e in word1 if e != '**']
                                node34=node34+node2
                                node44=node34
                        if len(sentencenext[0])==2:
                            write=2
                            uel_fatig_rcurve_list[-1].append(node44)
                            node34=[]
                            node44=[]
                            counter=99

#           if counter==99:
            if counter==5:
                break
#
    inputfile.close()
    params1=[[]]
    uel_mat_prop1=[[]]
    total_ele_set1=[[]]
    modf_ele_set=[[]]
    param2_arr=[[]]
    param2_arrn=[[]]
#    
    for index, value in enumerate(params):
        params1[index]= [x for x in value if len(x) > 0 ]   
    for index, value in enumerate(uel_mat_prop):
        uel_mat_prop1[index]= [x for x in value if len(x) > 0 ]   
    for index, value in enumerate(total_ele_set):
        total_ele_set1[index]= [x for x in value if len(x) > 0 ]   
#
    ltp=len(params1[0])
    param2_arr = params1[0][1:ltp:2]
    param2_arrn = params1[0][0:ltp:2]
#
#
    for index, value in enumerate(param2_arrn):
        for index1, value1 in enumerate(uel_mat_prop1):
            for index2, value2 in enumerate(value1):
                for index3, value3 in enumerate(value2):
                    if value3[0] == '<'  and value3[-1] == '>':
                        value31=value3[1:-1]
                    if value3[0] != '<'  and value3[-1] != '>':
                        value31=value3
                    if value[0].lower() == value31.lower():
                        temp=param2_arr[index]
                        uel_mat_prop1[index1][index2][index3]=temp[0]
#
    for indexn, valuen in enumerate(uel_el_set_name[0]):
        for indexn1, valuen1 in enumerate(total_ele_set_name[0]):
            if valuen[0].lower() == valuen1[0].lower():
                modf_ele_set[-1].append(total_ele_set1[0][indexn1])
#                modf_ele_set[-1].append(total_ele_set1[0][indexn1])
#
    uel_mat_list = uel_mat_prop1[0]
    uel_ele_set_list = modf_ele_set[0]
#
    mat_list_arr = []
    set_list_arr = []
    fatig_prop_list_arr = []
    fatig_prop_rcurve_list_arr = []

    for i in range(len(uel_ele_set_list)):
        set_list_arr.append(np.asarray(uel_ele_set_list[i],dtype=int))
    for i in range(len(uel_mat_list)):
        mat_list_arr.append(np.asarray(uel_mat_list[i],dtype=float))
    for i in range(len(uel_fatig_prop_list[0])):
        fatig_prop_list_arr.append(np.asarray(uel_fatig_prop_list[0][i],dtype=float))
    for i in range(len(uel_fatig_rcurve_list[0])):
        fatig_prop_rcurve_list_arr.append(np.asarray(uel_fatig_rcurve_list[0][i],dtype=float))
#
    fatig_prop_list_arr = reshape_table(fatig_prop_list_arr,6)
    fatig_prop_rcurve_list_arr = reshape_table(fatig_prop_rcurve_list_arr,3)
        
    return mat_list_arr,set_list_arr,fatig_prop_list_arr,fatig_prop_rcurve_list_arr

def reshape_table(table_in,no_columns):
    
    no_tables = len(table_in)
    table_out = []
    for i in range(no_tables):
        table = table_in[i]
        no_rows = int(len(table)/no_columns)
        table_out.append(table.reshape((no_rows,no_columns)))

    return table_out



ind_ini_dstat = 0

#indexes for orientation corresponding to theta i in the material list (*uel list)
ind_ini_theta_T = 1
ind_ini_theta_B = 2

#indexes to fracture toughness values
ind_ini_GIC_B   = 10
ind_ini_GIIC_B  = 11
ind_ini_BK_B   = 12

ind_ini_GIC_T    = 13
ind_ini_GIIC_T  = 14
ind_ini_BK_T   = 15

#indexes to strength values (static/miners)
ind_ini_Z_T = 6
ind_ini_Z_S = 8

#indexes to r-curves in the material list (*uel): currently only enabled for the interface
ind_ini_r_curve_B = 17
ind_ini_r_curve_T = 19
#indexes to dadn laws in the material list (*uel)
ind_ini_dadn_B = 16
ind_ini_dadn_T = 18

