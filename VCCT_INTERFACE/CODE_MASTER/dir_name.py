# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 09:47:02 2016

@author: ndecarva
"""
import os
import numpy as np
import sys
import platform

def read_info_from_python(cdirectory,input_file,variable):

    python_file = input_file + '.py'
    filename_in = cdirectory + '/'+python_file
    with open(filename_in,'U') as f_in:
        character_pos = 0
        for line in f_in:
            if line.find(str(variable)) != -1:
                if line.find('True') != -1:
                    variable = True
                elif line.find('False') != -1:
                    variable = False
                else:
                    for ch in line:
                        if ch == '=':
                            variable = line[character_pos+1:]
                        else:
                            character_pos += 1
                            break        
    return variable
                     


def read_file_name_from_input_stat():
    """obtains input file_name from input_stat.sh"""
    filename_in = cdirectory+'/input_stat.sh'
    with open(filename_in,'U') as f_in:
        character_pos = 0
        for line in f_in:
            if line.find('j=') != -1:
                for j in range(len(line)):
                    if line[j] == 'j': #begining of j=aksnksnd
                        for k in range(j,len(line)):
                            if np.all(line[k:k+4] != 'user'): #begining of user=...
                                character_pos +=1
                            else:
                                break
                        word_list = line.split()
                        input_file_n = word_list[2]
    input_file_n = input_file_n[6:] #removes input=
    
    return input_file_n

def writes_file_name_into_python(cdirectory,python_file,python3_bool):
    
    filename_in =  cdirectory + '/'+python_file
    filename_out = cdirectory + '/'+python_file+'New'
    f_out = open(filename_out,'wb') #writing in binry mode to avoid spurious ^m characters depending on the system
    stat = 0
    with open(filename_in,'U') as f_in: #opened in universal mode (to avoid ^M character)           
        for line in f_in:
            if line.find('import') != -1 and stat == 0:
                s_dummy = 'import '+input_file_n+' as inp_vars\n'
                f_out.write(s_dummy.encode('UTF-8'))
                stat = 1
            f_out.write(line.encode('UTF-8'))
        stat = 0
        if len(sys.argv)>2:
            if sys.argv[2] == 'abaverify':
                if line.find('###End NVCCT') != -1 and stat == 0:  
                    f_out.write('###Recall verts_path variable, calculate effective crack length and write to file "eff_da_length_vs_N.txt"\n'.encode('utf-8'))
                    f_out.write('    verts_path = inp_vars.verts_path\n'.encode('utf-8'))
                    f_out.write('    eff_da_length=eff_crack_length(all_nodes,all_elements,failed_els_I,verts_path)\n'.encode('utf-8'))
                    f_out.write('    no_chars = len(current_dir)\n'.encode('utf-8'))
                    if os.path.isdir(os.getcwd() +'/testOutput'):
                        #write this if 'testOutput' folder exists
                        f_out.write('    jobName = current_dir + "testOutput/" + input_file[no_chars:] + ".txt"\n'.encode('utf-8'))
                    else:
                        #write this if there is no 'testOutput' folder
                        f_out.write('    jobName = current_dir + input_file[no_chars:] + ".txt"\n'.encode('utf-8'))
                    f_out.write('    effdatofile(jobName,ntotal,eff_da_length)\n'.encode('utf-8'))
                    f_out.write('\n'.encode('utf-8'))
                    stat = 1 
    f_out.close()            
    os.remove(filename_in)
    os.rename(filename_out,filename_in)
    
    return
    
def writes_into_fortran(cdirectory,fortran_file,no_characters,nlgeom,contact):

    if nlgeom :
        nlgeom_s = '.True.'
    else:
        nlgeom_s = '.False.'

    if contact:
        contact = '.True.'
    else:
        contact = '.False.'

    
    filename_in =  cdirectory + '/'+fortran_file 
    filename_out = cdirectory + '/'+fortran_file + 'New'
    f_out = open(filename_out,'wb') #writing in binry mode to avoid spurious ^m characters depending on the system
    with open(filename_in,'U') as f_in: #opened in universal mode (to avoid ^M character)
            for line in f_in: 
                if line.find('nlgeom=') != -1:
                    for j in range(len(line)):
                        if line[j] == '=':
                            s_dummy = line[:j+1]+nlgeom_s+'\n'
                            f_out.write(s_dummy.encode('UTF-8'))   
                elif line.find('contact=') != -1:
                    for j in range(len(line)):
                        if line[j] == '=':
                            s_dummy = line[:j+1]+contact+'\n'                              
                            f_out.write(s_dummy.encode('UTF-8'))                                              
                        
                                                
                else:
                    f_out.write(line.encode())

                        
    f_out.close() 
    os.remove(filename_in)
    os.rename(filename_out,filename_in)
    return


def updates_k_c_user_in_fortran(cdirectory,fortran_file,k_c_user):

    filename_in =  cdirectory + '/'+fortran_file 
    filename_out = cdirectory + '/'+fortran_file + 'New'
    f_out = open(filename_out,'wb') #writing in binry mode to avoid spurious ^m characters depending on the system
    with open(filename_in,'U') as f_in: #opened in universal mode (to avoid ^M character)
            for line in f_in: 
                if line.find('k_c_user =') != -1:
                    for j in range(len(line)):
                        if line[j] == '=':
                            character_pos = j+1
                    f_out.write(line[:character_pos]+str(k_c_user)+'\n')
                else:
                    f_out.write(line)
    f_out.close() 
    os.remove(filename_in)
    os.rename(filename_out,filename_in)
    return
    

def updates_vis_output_freq_in_fortran(cdirectory,fortran_file,output_freq_v):
    
    filename_in =  cdirectory + '/'+fortran_file 
    filename_out = cdirectory + '/'+fortran_file + 'New'
    f_out = open(filename_out,'wb') #writing in binry mode to avoid spurious ^m characters depending on the system
    with open(filename_in,'U') as f_in: #opened in universal mode (to avoid ^M character)
            for line in f_in:
                if line.lower().replace(' ','').find('output_freq_vis=') != -1:
                    for j in range(len(line)):
                        if line[j] == '=':
                            character_pos = j+1
                    f_out.write(line[:character_pos]+str(output_freq_v)+'\n')
                else:
                    f_out.write(line)
    f_out.close() 
    os.remove(filename_in)
    os.rename(filename_out,filename_in)
    return


def updates_python_input_file(cdirectory,input_file_n):
    
        #updates pyhton inpt file "input_file_n.py"              
    filename_in =  cdirectory + "/"+input_file_n+".py"
    filename_out = cdirectory + "/"+input_file_n+"new.py" 
    f_out = open(filename_out,'wb') #writing in binry mode to avoid spurious ^m characters depending on the system
    with open(filename_in,'U') as f_in: #opened in universal mode (to avoid ^M character)            
        for line in f_in:
            character_pos = 0
            if line.find('current_dir =') != -1:
                for j in range(len(line)):
                    if line[j] == '=':
                        character_pos = j+1
                        s_dummy = line[:character_pos]+"'"+cdirectory[:]+"/'\n"
                        f_out.write(s_dummy.encode('UTF-8'))  
            elif line.find('input_file =') != -1:
                for j in range(len(line)):
                    if line[j] == '=':
                       character_pos = j+1
                       s_dummy = line[:character_pos]+"'"+input_file_n+"'\n"
                       f_out.write(s_dummy.encode('UTF-8'))  
            else:
                f_out.write(line.encode())
                #f_out.write(line)
    f_out.close()              
    os.remove(filename_in)
    os.rename(filename_out,filename_in)   
    
    return


def updates_shell_input_file(cdirectory,shell_update,python_file):
    
        #updates shell file: update_elements_cz.sh
    filename_in = cdirectory +"/"+shell_update
    filename_out = cdirectory +"/"+shell_update+'New'
    f_out = open(filename_out,'wb')
    with open(filename_in,'U') as f_in: #opened in universal mode (to avoid ^M character)            
        for line in f_in:
            if line[:6].find('python') != -1:
                s_dummy = line[:7]+cdirectory[:]+"/"+python_file+"\n"
                f_out.write(s_dummy.encode('UTF-8'))  
            else:
                f_out.write(line.encode())
    f_out.close()
    
    os.remove(filename_in)
    os.rename(filename_out,filename_in)
    
    return
    
cdirectory=os.getcwd()
    
if __name__ == '__main__':

    cdirectory=os.getcwd()
    
    try:
        input_file_n = sys.argv[1]
    except:
        input_file_n = read_file_name_from_input_stat()

    python3_bool = False
    if int(platform.python_version()[0]) >= 3:
        python3_bool = True
        
    writes_file_name_into_python(cdirectory,'NVCCT_I.py',python3_bool)

    no_characters=len(cdirectory)
                   
    updates_python_input_file(cdirectory,input_file_n)     
    
    inp_vars = __import__(input_file_n)

    
    

    #reads nlgeom info from python, nlgeom:
    nlgeom = False    
    if hasattr(inp_vars,'nlgeom'):
        nlgeom = inp_vars.nlgeom
    print('nlgeom',nlgeom)

    #reads contact info from python, interface:
    contact = False    
    if hasattr(inp_vars,'contact'):
        contact = inp_vars.contact
    print('contact',contact)

    writes_into_fortran(cdirectory,'UEL3D_INTERFACE_NVCCT.f',no_characters,nlgeom,contact)
        
    if hasattr(inp_vars, 'output_vis_freq'):
        updates_vis_output_freq_in_fortran(cdirectory,'UEL3D_INTERFACE_NVCCT.f',inp_vars.output_vis_freq)
    if hasattr(inp_vars, 'k_c_user'):
        updates_k_c_user_in_fortran(cdirectory,'UEL3D_INTERFACE_NVCCT.f',inp_vars.k_c_user)
        
    updates_shell_input_file(cdirectory,'update_elements_nvcct.sh','NVCCT_I.py')
                        
                       
