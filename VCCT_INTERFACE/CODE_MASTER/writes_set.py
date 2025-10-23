# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 10:19:20 2022

@author: ndecarva
"""
import numpy as np
import os

def read_elements_cze(filename):
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
                                all_els_dict[s[0]] = {'nodes':s[1:]}
                        else: 
                            break
                            
    return all_els_dict



def writeset(filename,settype,setname,members,instance_name):
    with open(filename,'w') as f:
        if len(instance_name) == 0:
            f.write('*'+settype+','+ settype+'='+setname+'\n') 
        else:
            f.write('*'+settype+','+ settype+'='+setname+',instance='+instance_name+'\n') 
        count = 0
        for i in members:
            count = count + 1
            if count == 1:
                f.write(str(i))
            else:
                f.write(','+str(i))
            if count == 16:
                f.write('\n')
                count = 0
        if count != 16:
            f.write('\n') 

def splitline(i,maxnopos,maxnochar):    
    leni_vect = len(i) #lenght of the vector
    if leni_vect > maxnopos:
        iout = i[0:maxnopos]
        irest = i[maxnopos:leni_vect]
    else:
        iout = i
        irest = []
         
    leniout_str = len(','.join(['%d' % num for num in iout])+',\n')
    leniout_vect = len(iout)
    if leniout_str > maxnochar:
        #estimate an approximate ratio for decreasing the size of the vector 
        ratio = leniout_str/float(maxnochar) 
        leniout_vect = leniout_vect*ratio
        #integer 
        leniout_vect = int(leniout_vect)
        iout = iout[0:leniout_vect] #updates the vector size
        leniout_str = len(','.join(['%d' % num for num in iout])+',\n') #forms a new string
        
        while leniout_str > maxnochar: #determines whether the size of the new string is less than maxnochr, there could be a difference since the ration was determined based on the size of the string but applied to the vector components
            iout = iout[0:leniout_vect-1]  #decreases the size of the vector by onw
            irest = i[leniout_vect-1:leni_vect] #updates the "rest" vector
            leniout_str = len(','.join(['%d' % num for num in iout])+',\n') #computes the new string size from the new vector
            leniout_vect = len(iout) #computes the new vector lenggh
            
    return iout, irest


            
def writeels (filename,elements,eltype,elset = ''):
    
    maxnopos = 16 #maximum number of positions (integers - nodes and element number) that can be written per line    
    maxnochar = 80 #maximum number of characters perline

    with open(filename,'w') as f:
        f.write('*Element,'+eltype+','+elset+'\n')
        for i in elements:
            iout = i
            while (len(iout) > maxnopos) or (len(','.join(['%d' % num for num in iout])) > maxnochar):
                iout,irest =splitline(iout,maxnopos,maxnochar)
                f.write(','.join(['%d' % num for num in iout])+',\n')
                if len(irest) != 0:
                    iout = irest
            f.write(','.join(['%d' % num for num in iout])+'\n')
    
    
def write_viz_file(current_dir,real_el_no,ghost_el_no):
    filename = current_dir + '\\viz.uel'
    with open(filename,'w') as f:
        f.write(str(int(real_el_no))+'\n')
        f.write(str(int(ghost_el_no))+'\n')
    
    return 

if __name__ == "__main__":
    
    
    x = input("To generate elset insert (1) to generate ghost elements insert (2): \n") 
    if x == '1':
        start_el_no = int(input('insert start el no: \n'))
        end_el_no  = int(input('insert end el no: \n'))
        print('set written to file elset.txt\n')
        
        elset = np.arange(start_el_no, end_el_no)
        filename = 'elset.txt'
        writeset(filename,'*ELSET','','',elset,[])
        
    elif x == '2':
   
        f_name = input('insert filename including .inp: \n')
        max_part_el_no = int(input('insert maximum element number in part: \n'))
        
        current_dir  = os.getcwd()
        filename = current_dir + '\\' + f_name
        
        
        all_els_dict = read_elements_cze(filename)
        
        
        all_els_np_arr = np.array([np.concatenate((np.array([el_no]),nodes['nodes'])) for el_no, nodes in all_els_dict.items()])
        
        all_els_np_arr = all_els_np_arr[np.argsort(all_els_np_arr[:,0]),:]
        
        min_vcct_el_no =  all_els_np_arr[0,0]
        max_vcct_el_no = all_els_np_arr[-1,0]
        all_els_np_arr[:,0] += (max_part_el_no -min_vcct_el_no)+1
        
        writeels(filename[:-4] + 'ghost.elements',all_els_np_arr,'type=','elset=ghost')
        
        print('ghost elements created and written to file: '+filename[:-4] + 'ghost.elements\n')
        
        write_viz_file(current_dir,min_vcct_el_no,max_part_el_no+1)
            
        
        
        
        