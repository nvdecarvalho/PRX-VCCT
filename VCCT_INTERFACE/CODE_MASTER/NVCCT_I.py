# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 14:39:06 2022

@author: ndecarva
"""
import time
import warnings
import sys
import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection 
import initialize_cze2vcct as ini
#import test_NVCCT_I as tn

class line2D:
    
    def __init__(self):
        self.n = np.array([0.0,0.0])
        self.point = np.array([0.0,0.0])
    
    def line_two_points(self,x0,x1):   
        self.n = x1-x0
        self.point = x0
        
    def slope(self):
        #calculates the slope of the line "m" in y = mx+b
        if self.n[0] != 0:
            m = self.n[1]/self.n[0]
        else:
            m = np.tan(np.pi/2.0)         
        return m
    
    def y_axis_int(self):
        #calculates the intersection with the y axis
        p_x = self.point[0]
        p_y = self.point[1]
        n_x = self.n[0]
        n_y = self.n[1]
        y_int = np.array([])
        if n_x != 0:
            t = -p_x/n_x
            y_int = p_y + n_y*t
        
        return y_int
    
    def intersect_w_line(self,line1):
        co_linear  = False
        x0  = self.point[0]
        y0  = self.point[1]
        n0x = self.n[0]
        n0y = self.n[1]
        x1  = line1.point[0]
        y1  = line1.point[1]
        n1x = line1.n[0]
        n1y = line1.n[1]
        t0 = None    
        cross_p = -n0x*n1y+n1x*n0y
        if cross_p == 0:
            #lines are parallel
            int_l = np.array([])
             
            v0x = x0-x1
            v0y = y0-y1
            
            cross_pv = -v0x*n1y+n1x*v0y
            if cross_pv == 0:
                #lines are coincident (all points are concident)
                co_linear = True
                
        else:
            t0 = (-(x1-x0)*n1y+n1x*(y1-y0))/cross_p
            #t1 = (-(x1-x0)*n0y+n0x*(y1-y0))/cross_p
            
            x_int0 = x0 + n0x*t0
            y_int0 = y0 + n0y*t0
            
#            x_int1 = x1 + n1x*t1
#            y_int1 = y1 + n1y*t1
            #print (x_int0,x_int1,y_int0,y_int1)
 
            int_l = np.array([x_int0,y_int0])
            
        return int_l,t0,co_linear


def point_in_element_asc(el_no,point,all_els_dict,all_nodes_dict):
    """determines whether a given point (2d) is inside an element,
    if the point is located over a corner/edge is considered to be inside
    test case: point_in_element"""
    
    tol = 1.0e-4 #defined by abaqus' node coordinate precision
    tol2 = 1e-16
    tol_arr = np.array([tol,tol])

    nodes_el = all_els_dict[el_no]['nodes']
    a_arr = np.array([all_nodes_dict[node]['coords'][:2] for node in nodes_el])

    #checks corners:
    for i,corner in enumerate(a_arr):
        if np.all(abs(point - corner) < tol_arr):
            return True
    
    #checks inside
    for i,corner in enumerate(a_arr):
        edge = a_arr[index_roll(i+1,4),:] - corner
        edge_n = edge/np.linalg.norm(edge)
        v_trial = corner-point
        v_trial_n = v_trial/np.linalg.norm(v_trial)    
        if np.cross(v_trial_n,edge_n) < -tol:
            return False
        if abs(np.dot(edge_n,-v_trial_n) - 1.0) <tol2:
            if np.linalg.norm(v_trial)/np.linalg.norm(edge) <= 1.0-tol:
                return True
    return True
    
def point_in_element(el_no,point,all_els_dict,all_nodes_dict):
    """determines whether a given point (2d) is inside an element,
    if the point is located over a corner/edge is considered to be inside
    test case: point_in_element"""
    
    tol = 1.0e-4 #defined by abaqus' node coordinate precision
    tol2 = 1e-16
    tol_arr = np.array([tol,tol,tol])

    nodes_el = all_els_dict[el_no]['nodes']
    a_arr = np.array([all_nodes_dict[node]['coords'] for node in nodes_el])

    #checks corners:
    for i,corner in enumerate(a_arr):
        if np.all(abs(point - corner) < tol_arr):
            return True
    
    #checks inside
    for i,corner in enumerate(a_arr):
        edge = a_arr[index_roll(i+1,4),:] - corner
        edge_n = edge/np.linalg.norm(edge)
        v_trial = corner-point
        v_trial_n = v_trial/np.linalg.norm(v_trial)    
        if np.cross(v_trial_n,edge_n) < -tol:
            return False
        if abs(np.dot(edge_n,-v_trial_n) - 1.0) <tol2:
            if np.linalg.norm(v_trial)/np.linalg.norm(edge) <= 1.0-tol:
                return True
    return True


def project_point_on_plane(p_normal,p_ref_point,point):
    
    v = point - p_ref_point
    v_n = np.dot(v,p_normal)
    proj_point = point-v_n*p_normal
    
    return proj_point

def calc_normal2el(el,all_els_dict,all_nodes_dict):
    
    nodes = all_els_dict[el]['nodes']
    nodes_coords = np.array([all_nodes_dict[n]['coords'] for n in nodes]) 
    v1 = nodes_coords[1,:] -  nodes_coords[0,:]
    v2 = nodes_coords[2,:] - nodes_coords[1,:]
    
    normal2el = np.cross(v1,v2)
    normal2el = normal2el/np.linalg.norm(normal2el)
    
    return normal2el


def inverse_iso_p_asc(el,xp,all_els_dict,all_nodes_dict):
    """determines the natural coordinates given the real coordinates 
      via newton raphson   
    """
    nodes = all_els_dict[el]['nodes']
    xy_nodes = np.array([all_nodes_dict[n]['coords'][:2] for n in nodes])        
      
    #f(nat) = xp - N(nat)x_nodes
    #find f(nat) = 0
    #implements newton-raphson
    #nat = nat_prev - inv(df(nat_prev)/dnat)*f(nat_prev)
    nat_trial = np.array([0,0])
    tol = 1.0E-4
    # if el == 1573:
    #     print('here')
            
    for i in range(100):
        n_shape = shape_f(nat_trial,'quad_linear')
        xp_nat_trial = np.array([np.dot(n_shape,xy_nodes[:,0]),np.dot(n_shape,xy_nodes[:,1])])
        #evaluates f(nat)
        f_nat = xp - xp_nat_trial
        #if larger than tolerance for zero:
        if np.linalg.norm(f_nat) > tol:
            
            #calculates df(nat_prev)/dnat = -dx_dn (via shape functions)
            n_dshape = dshape_f(nat_trial,'quad_linear')
            dx_dn = -np.array([[np.dot(n_dshape[:,0],xy_nodes[:,0]),np.dot(n_dshape[:,1],xy_nodes[:,0])],
                              [np.dot(n_dshape[:,0],xy_nodes[:,1]),np.dot(n_dshape[:,1],xy_nodes[:,1])]])           
            dn_dx = np.linalg.inv(dx_dn)
            
            #calculates next trial: nat = nat_prev - inv(df(nat_prev)/dnat)*f(nat_prev) 
            next_nat_trial = nat_trial - np.dot(dn_dx,f_nat)
            nat_trial = next_nat_trial
            #nat_trial = np.clip(nat_trial,-1.0,1.0)
        else:
            return nat_trial
        
    print("problem sampling")
    return None


def inverse_iso_p(el,xp,all_els_dict,all_nodes_dict):
    """determines the natural coordinates given the real coordinates 
      via newton raphson   
    """
    nodes = all_els_dict[el]['nodes']
    xy_nodes = np.array([all_nodes_dict[n]['coords'] for n in nodes])        
      
    #f(nat) = xp - N(nat)x_nodes
    
    #find f(nat) = 0
    
    #implements newton-raphson
    #nat = nat_prev - inv(df(nat_prev)/dnat)*f(nat_prev)
    nat_trial = np.array([0,0])
    tol = 1.0E-4
    # if el == 1573:
    #     print('here')
            
    for i in range(100):
        n_shape = shape_f(nat_trial,'quad_linear')
        xp_nat_trial = np.array([np.dot(n_shape,xy_nodes[:,0]),np.dot(n_shape,xy_nodes[:,1]),np.dot(n_shape,xy_nodes[:,2])])
        #evaluates f(nat)
        f_nat = xp - xp_nat_trial
        #if larger than tolerance for zero:
        if np.linalg.norm(f_nat) > tol:
            
            #calculates df(nat_prev)/dnat = -dx_dn (via shape functions)
            n_dshape = dshape_f(nat_trial,'quad_linear')
            dx_dn = -np.array([[np.dot(n_dshape[:,0],xy_nodes[:,0]),np.dot(n_dshape[:,1],xy_nodes[:,0])],
                              [np.dot(n_dshape[:,0],xy_nodes[:,1]),np.dot(n_dshape[:,1],xy_nodes[:,1])],
                              [np.dot(n_dshape[:,0],xy_nodes[:,2]),np.dot(n_dshape[:,1],xy_nodes[:,2])]])     
            dx_dn_t = np.transpose(dx_dn)
            dn_dx = np.matmul(np.linalg.inv(np.matmul(dx_dn_t,dx_dn)),dx_dn_t)
            
            #dn_dx = np.linalg.inv(dx_dn)
            
            #calculates next trial: nat = nat_prev - inv(df(nat_prev)/dnat)*f(nat_prev) 
            next_nat_trial = nat_trial - np.dot(dn_dx,f_nat)
            nat_trial = next_nat_trial
            #nat_trial = np.clip(nat_trial,-1.0,1.0)
        else:
            return nat_trial
        
    print("problem sampling")
    return None



def shape_f(xnat,s_type):
    shape_v = None
    if s_type == 'quad_linear':
        shape_v = np.array([0.25*(1-xnat[0])*(1-xnat[1]),
                 0.25*(1+xnat[0])*(1-xnat[1]),
                 0.25*(1+xnat[0])*(1+xnat[1]),
                 0.25*(1-xnat[0])*(1+xnat[1])])
    
    return shape_v

def dshape_f(xnat,s_type):
    
    dshape_m = None
    if s_type == 'quad_linear':
        dshape_m = np.array([[-0.25*(1-xnat[1]),-0.25*(1-xnat[0])],
                          [0.25*(1-xnat[1]),-0.25*(1+xnat[0])],
                          [0.25*(1+xnat[1]),0.25*(1+xnat[0])],
                          [-0.25*(1+xnat[1]),0.25*(1-xnat[0])]])

    return dshape_m

def plot_elements_dict(ax,all_els_dict,all_nodes_dict,coords_indexes,color):
    """plots a mesh"""
    vertices = []
    for el,nodes in all_els_dict.items():
        v_el = np.array([all_nodes_dict[node]['coords'][coords_indexes]  for node in nodes['nodes']]) 
        vertices += [v_el]
    coll = PolyCollection(vertices,facecolors = color,edgecolors = 'black')
    ax.add_collection(coll)
    return coll

#



def bisect_vectors_asc(vect1,vect2,ref_vect):
    """obtains the vector bisecting two vectors vect1,vect2 with direction along
    ref vect"""
    
    norm_vect1 = np.linalg.norm(vect1)
    norm_vect2 = np.linalg.norm(vect2)
    
    b_vect = norm_vect2*vect1 + norm_vect1*vect2
    b_vect_norm = np.linalg.norm(b_vect)
    
    if b_vect_norm == 0:
        if norm_vect1 != 0:
            b_vect = np.array([-vect1[1],vect1[0]])/norm_vect1           
        elif norm_vect2 != 0:
            b_vect = np.array([-vect2[1],vect2[0]])/norm_vect2
        else:
            b_vect = vect1
    else:
        b_vect = b_vect/b_vect_norm
        if np.dot(b_vect,ref_vect) < 0:
            b_vect = -b_vect

    return b_vect


def bisect_vectors(vect1,vect2,ref_vect):
    """obtains the vector bisecting two vectors vect1,vect2 with direction along
    ref vect"""
    
    norm_vect1 = np.linalg.norm(vect1)
    norm_vect2 = np.linalg.norm(vect2)
    
    b_vect = norm_vect2*vect1 + norm_vect1*vect2
    b_vect_norm = np.linalg.norm(b_vect)
    if b_vect_norm == 0:
        if norm_vect1 != 0:
            c_vect = np.cross(ref_vect,vect1)
            b_vect = np.cross(vect1,c_vect)
            b_vect_norm = np.linalg.norm(b_vect)
            #b_vect = np.array([-vect1[1],vect1[0]])/norm_vect1           
        elif norm_vect2 != 0:
            c_vect = np.cross(ref_vect,vect2)
            b_vect = np.cross(vect2,c_vect)
            b_vect_norm = np.linalg.norm(b_vect)
            #b_vect = np.array([-vect2[1],vect2[0]])/norm_vect2
        else:
            b_vect = ref_vect
    
    b_vect = b_vect/b_vect_norm
    if np.dot(b_vect,ref_vect) < 0:
        b_vect = -b_vect

    return b_vect

def n2vector(vect,ref_v):
    """obtains the normalized vector perpendicular to a 2d vector vect, 
    along ref vect ref vect"""
    norm_v = np.linalg.norm(vect)
    if norm_v != 0:
        nv = np.array([vect[1],-vect[0]])/norm_v
        if np.dot(ref_v,nv) < 0:
            nv = -nv
    else:
        nv = vect
        
    return nv




def calc_normals_2_crack_asc(vcct_el,a_node_p_along_normal,vcct_el_dict,all_nodes_dict,option):
    """computes three tentative normals:
    one: bisection of crack faces
    two/three: normal to either crack face
    remaining: ref_v"""

    if option == 'carvalho' or option == 'carvalho2':
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        
        ref_v = all_nodes_dict[a_node_p_along_normal[0]]['coords'][:2] -  all_nodes_dict[vcct_el]['coords'][:2]
        el_type = vcct_el_type(vcct_el,vcct_el_dict)
        #If VCCT element is at an edge or corner:
        if (el_type == 'e' or el_type == 'c'): 
            normals_2_crack = np.zeros((2,2))
            normals_2_crack[0,:] = ref_v/np.linalg.norm(ref_v)
        else:
            if option == 'carvalho':
                normals_2_crack = np.zeros((3,2))
                normals_2_crack[0,:] = bisect_vectors_asc(cf_v[0,:],cf_v[1,:],ref_v)
                for i,cv in enumerate(cf_v):
                    normals_2_crack[1+i,:] = n2vector(cv,ref_v)
            elif option == 'carvalho2':
                normals_2_crack = np.array([bisect_vectors_asc(cf_v[0,:],cf_v[1,:],ref_v)])
            
            
    elif option == 'boeing':
        a_f_ns = a_nodes_cond(vcct_el,vcct_el_dict,all_nodes_dict,'deq1')    
        vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
        normals_2_crack = np.array([vcct_el_xy - all_nodes_dict[a_f]['coords'][:2] for a_f in a_f_ns])
        normals_2_crack = np.array([n/np.linalg.norm(n) for n in normals_2_crack])
    else:
        print('---option invalid in calc_normals_2_crack()---')
        normals_2_crack = []
        
        
    return normals_2_crack


def calc_normals_2_crack(vcct_el,a_node_p_along_normal,vcct_el_dict,all_nodes_dict,option):
    """computes three tentative normals:
    one: bisection of crack faces
    two/three: normal to either crack face
    remaining: ref_v"""

    if option == 'carvalho' or option == 'carvalho2':
        cf_v = all_nodes_dict[vcct_el]['cf_v']
      
        ref_v = all_nodes_dict[a_node_p_along_normal]['coords'] -  all_nodes_dict[vcct_el]['coords']
        el_type = vcct_el_type(vcct_el,vcct_el_dict)
        #If VCCT element is at an edge or corner:
        if (el_type == 'e' or el_type == 'c'): 
            
            normals_2_crack = ref_v/np.linalg.norm(ref_v)
        else:
            if option == 'carvalho':
                normals_2_crack = np.zeros((3,3))
                normals_2_crack[0,:] = bisect_vectors(cf_v[0,:],cf_v[1,:],ref_v)
                for i,cv in enumerate(cf_v):
                    normals_2_crack[1+i,:] = n2vector(cv,ref_v)
            elif option == 'carvalho2':
                normals_2_crack = bisect_vectors(cf_v[0,:],cf_v[1,:],ref_v)
            
            
    elif option == 'boeing':
        a_f_ns = a_nodes_cond(vcct_el,vcct_el_dict,all_nodes_dict,'deq1')    
        vcct_el_xy = all_nodes_dict[vcct_el]['coords']
        normals_2_crack = np.array([vcct_el_xy - all_nodes_dict[a_f]['coords'] for a_f in a_f_ns])
        normals_2_crack = np.array([n/np.linalg.norm(n) for n in normals_2_crack])
    else:
        print('---option invalid in calc_normals_2_crack()---')
        normals_2_crack = []
        
        
    return normals_2_crack    
    
def vcct_el_type(vcct_el,vcct_el_dict):
    
    no_nodes = np.sum(vcct_el_dict[vcct_el] > 0)
    el_type = None
    if no_nodes == 8:
        el_type = 'm' #middle
    elif no_nodes == 5:
        el_type = 'e' #edge
    elif no_nodes == 3:
        el_type = 'c' #corner
    
    return el_type

def index_roll(ind,indmax):
        
    return ind - int(ind/indmax)*indmax


def antenna_nodes_per_el(node_i,el_j):
    #given a set of nodes find the antenna nodes to a given node n based on the connectivity of the element
    no_nodes = 4
    node_ind = np.nonzero(el_j == node_i)[0]
    antena_nodes_per_el = np.zeros(2)

    antena_nodes_per_el[0] = el_j[index_roll(node_ind-1,no_nodes)]
    antena_nodes_per_el[1] = el_j[index_roll(node_ind+1,no_nodes)]

    return antena_nodes_per_el


def ref_cclock_el_2d(el_nodes,all_nodes_dict):
    """computes a sign associated with the counter clockwise normal to the element
    asssumes xy coordinates only"""
    
    v12 = all_nodes_dict[el_nodes[1]]['coords'] - all_nodes_dict[el_nodes[0]]['coords'] 
    v23 = all_nodes_dict[el_nodes[2]]['coords'] - all_nodes_dict[el_nodes[1]]['coords']
 
    return np.sign(np.cross(v12,v23))     



def sample_points_coords(vcct_el,vcct_el_coords,vcct_adj_els,node_ip1,a_node_wake,ind_a_node_wake,ns2crack,vcct_el_dict,
                     node_2vcct_adj_els_dict,all_els_dict,all_nodes_dict,option='boeing-mod'):
    
    tol = 1E-6
    
    if option == 'boeing':
        #four potential sample points, the maximum will be used
        
        no_antena_nodes= len(a_nodes_indexes)
        sample_point = np.zeros((no_antena_nodes,3))
        sample_point_f = np.zeros((no_antena_nodes,3))
        for i,a_ind in enumerate(a_nodes_indexes):
            
            sample_point[i,:] = 2*vcct_el_coords - all_nodes_dict[vcct_adj_els[a_ind]]['coords']
            sample_point_f[i,:] =all_nodes_dict[vcct_adj_els[a_ind]]['coords']
            
    elif option == 'boeing-mod':
        #uses the antena node identified as the next node pair based on the wake node
        #wake node identified based on whether the antena node is failed and the crack opening
        sample_point = np.array([2*vcct_el_coords - all_nodes_dict[node_ip1]['coords']])
        sample_point_f = np.array([all_nodes_dict[node_ip1]['coords']])
        
    elif option == 'normal' or option == 'normal_along_antena':
        el_type = vcct_el_type(vcct_el,vcct_el_dict)
        node_ip1_coords = all_nodes_dict[node_ip1]['coords'] 
        if el_type == 'e':
            sample_point = np.array([all_nodes_dict[a_node_wake]['coords']])
            sample_point_f = np.array([node_ip1_coords])
        else:
            inc_arr = np.array([1,-1])
            
            for inc in inc_arr:
                #find intersection point between normal and diagonal conecting next node 
                #to the two antena nodes on either side
                
                #antena node based on "inc"
                ind_a_node = next_a_node_ind_by_inc(a_nodes_indexes,ind_a_node_wake,inc) 
                a_node = vcct_el_dict[vcct_el][ind_a_node]
                a_node_coords = all_nodes_dict[a_node]['coords']
                
                #find element sharing antena node and node_ip1:
                nodes_arr = np.array([node_ip1,a_node])
                adj_el = adj_el_to_vcct_node_from2nodes(vcct_el,nodes_arr,node_2vcct_adj_els_dict,all_els_dict)
                
                #find coordinates of, a_node,node_ip1 and vcc_el in the natural space of the element
                a_node_nat_coords = inverse_iso_p(adj_el,a_node_coords,all_els_dict,all_nodes_dict)
                node_ip1_nat_coords = inverse_iso_p(adj_el,node_ip1_coords,all_els_dict,all_nodes_dict)
                vcct_el_nat_coords = inverse_iso_p(adj_el,vcct_el_coords,all_els_dict,all_nodes_dict)
                
                #find coordinates of a point along the normal (may need to 
                #project normal into the element plane):
                #print('vcct_el_coords,ns2crack',vcct_el_coords,ns2crack)
                p_along_normal_coords = vcct_el_coords + ns2crack
                p_along_normal_nat_coords = inverse_iso_p(adj_el,p_along_normal_coords,all_els_dict,all_nodes_dict)
                    
                #diagonal line in the natural space:
                line_diagonal_nat = line2D()
                line_diagonal_nat.line_two_points(a_node_nat_coords,node_ip1_nat_coords)
                
                #normal vector in the natural space:
                line_along_normal_vector = line2D()
                line_along_normal_vector.line_two_points(vcct_el_nat_coords,p_along_normal_nat_coords)
            
                intersection_nat_coords,t0, c_linear = line_diagonal_nat.intersect_w_line(line_along_normal_vector)
                if ((len(intersection_nat_coords) != 0) and (t0 >=0-tol) and (t0<=1.0+tol)):
                    sample_point_f = np.array([interpol_el_nodal_val_wnat_coords('coords',intersection_nat_coords,adj_el,3,
                                                                                 'quad_linear',all_els_dict,all_nodes_dict)])
                    sample_point = np.array([-sample_point_f[0]+2.0*vcct_el_coords])
            if option == 'normal_along_antena':
                sample_point_f = np.array([node_ip1_coords])
                    

    return sample_point, sample_point_f
                
    
    
def disp_sample_pos(vcct_el,node_ip1,a_node_ref,ind_a_node_ref,ns2crack,vcct_el_dict,
                    node_2vcct_adj_els_dict,all_els_dict,all_nodes_dict,option='carvalho'):
    """compute the location of the sampling points the element planar coordinates
    calculation is based on the vectors describing the crack front, node i the next node
    node_ip1"""
    
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
    no_ns2crack = len(ns2crack)
    sample_point = np.zeros((no_ns2crack,2))    
    if option == 'boeing':
        for i,n in enumerate(node_ip1):
            sample_point[i,:] = 2*vcct_el_xy - all_nodes_dict[node_ip1[i]]['coords'][:2]
            
    elif option == 'carvalho':
        
        el_type = vcct_el_type(vcct_el,vcct_el_dict)
        if el_type == 'e':
            sample_point = np.array([all_nodes_dict[a_node_ref]['coords'][:2]])
        else:
            cf_v = all_nodes_dict[vcct_el]['cf_v']
            c_face1 = line2D()
            c_face1.n = cf_v[0,:]
            c_face1.point = vcct_el_xy
            
            c_face2 = line2D()
            c_face2.n = cf_v[1,:]
            c_face2.point = vcct_el_xy
            
       
            #calculates the intersection between the normals to the crack front and the crack front faces
            #with that intersection if any defines the sampling vector which subtracted to the node_i position
            #gives the sampling point
            for i,n in enumerate(ns2crack):
                node_ip1_xy = all_nodes_dict[node_ip1[i]]['coords'][:2]
                n_2face = line2D() 
                n_2face.point = node_ip1_xy     
                n_2face.n = n
                
                int_l_1,t_1,col_1 = c_face1.intersect_w_line(n_2face)
                int_l_2,t_2,col_2 = c_face2.intersect_w_line(n_2face)
                if t_1 != None:
                    if t_1 >=0.0 and t_1 < 1.0:
                         sample_point[i,:] = int_l_1-node_ip1_xy + vcct_el_xy
                         
                if t_2 != None:
                    if t_2 >= 0.0 and t_2 < 1.0:
                        sample_point[i,:] = int_l_2-node_ip1_xy + vcct_el_xy
    elif option == 'carvalho2':
        el_type = vcct_el_type(vcct_el,vcct_el_dict)
        if el_type == 'e':
            sample_point = np.array([all_nodes_dict[a_node_ref]['coords'][:2]])
        else:
            cf_v = all_nodes_dict[vcct_el]['cf_v']
            a_nodes_indexes = [0,2,4,6]
          
            n2crack = ns2crack[0] #only one normal considered
            
            ind_a_node_s1 = a_nodes_indexes[index_roll(ind_a_node_ref+1,4)]
            a_node_s1 = vcct_el_dict[vcct_el][ind_a_node_s1]
            a_node_s1_xy = all_nodes_dict[a_node_s1]['coords'][:2]
            ind_a_node_s2 = a_nodes_indexes[index_roll(ind_a_node_ref-1,4)]
            a_node_s2 = vcct_el_dict[vcct_el][ind_a_node_s2]
            a_node_s2_xy =  all_nodes_dict[a_node_s2]['coords'][:2]
            node_ip1_xy =  all_nodes_dict[node_ip1[0]]['coords'][:2]

            l1 = np.linalg.norm(a_node_s1_xy - vcct_el_xy)
            l2 = np.linalg.norm(a_node_s2_xy - vcct_el_xy)
            lfwd = np.linalg.norm(node_ip1_xy - vcct_el_xy)

            lmin = min(l1,l2,lfwd)/10.0
            
            point_along_normal = n2crack*lmin+vcct_el_xy

            adj_els = node_2vcct_adj_els_dict[vcct_el]

            adj_el_with_point = []
            for a_el in adj_els:
                if point_in_element_asc(a_el,point_along_normal,all_els_dict,all_nodes_dict):
                    adj_el_with_point = a_el 
                    break

            point_along_normal_nat = inverse_iso_p_asc(adj_el_with_point,point_along_normal,all_els_dict,all_nodes_dict)
            vcct_el_nat = inverse_iso_p_asc(adj_el_with_point,vcct_el_xy,all_els_dict,all_nodes_dict)
            # a_node_ip1_nat = vi.inverse_iso_p(adj_el_with_point,a_node_ip1_xy,all_els_dict,all_nodes_dict)

            # v_fwd_nat = (a_node_ip1_nat - vcct_nat)
            # v_fwd_nat_norm = v_fwd_nat/np.linalg.norm(v_fwd_nat)

            v_norm_nat = (point_along_normal_nat-vcct_el_nat)
            v_norm_nat_norm = v_norm_nat/np.linalg.norm(v_norm_nat)
            sp_nat = v_norm_nat_norm*2.0 + vcct_el_nat
            sp_xy = interpol_el_nodal_val_wnat_coords('coords',sp_nat,adj_el_with_point,3,'quad_linear',all_els_dict,all_nodes_dict)
            sample_point = np.array([-sp_xy[:2]+2.0*vcct_el_xy])
        
                            
    return sample_point
    
def plot_point(point,indexes,axes, color='k'):
    
    axes.plot(point[indexes[0]],point[indexes[1]],'o'+color)
    

def plot_vect(point,vect,indexes,axes,color='b'):
    vect_for_plot = np.zeros((2,2))
    vect_for_plot[0,:] = point[indexes]
    vect_for_plot[1,:] = vect[indexes]+point[indexes]
    
    axes.plot(vect_for_plot[:,0],vect_for_plot[:,1],color)
    axes.plot(vect_for_plot[0,0],vect_for_plot[0,1],'o'+color)


def plot_crack_front_el(cf_v,vcct_el_xy,coord_indexes,axes):
    for vect in cf_v:
        plot_vect(vcct_el_xy,vect,coord_indexes,axes)
    

def build_all_els_dict(all_elements):
    all_els_dict = {}
    for el in all_elements:
          all_els_dict[int(el[0])] ={'nodes':el[1:5].astype(int)}
          
    return all_els_dict

def build_all_nodes_dict(all_nodes):

    all_nodes_dict = {}    
    for node in all_nodes:
          all_nodes_dict[int(node[0])] = {'coords':node[1:4],'d':0.0}
    return all_nodes_dict
        
def build_node_2vcct_adj_els_dict(all_nodes_dict,all_els_dict):
    
    node_2vcct_adj_els = {}

    for node in all_nodes_dict.keys():
        node_2vcct_adj_els[node] = []
        count_adj_els = 0
        for el_no, el_nodes in all_els_dict.items():
            if node in el_nodes['nodes']:
                node_2vcct_adj_els[node] += [el_no]
                count_adj_els += 1
                if count_adj_els == 4:
                    break
    return node_2vcct_adj_els


def build_node_2vcct_adj_els_dict_fast(all_els_dict):
    nodes_dict = {}
    for el,el_nodes in all_els_dict.items():
        for node in el_nodes['nodes']:
            if node in nodes_dict:
                nodes_dict[node] += [el]
            else:
                nodes_dict[node] = [el]
    
    return nodes_dict

# def all_vcct_el_nodes_failed(vcct_el_nodes,all_nodes_dict):
   
#     for node in vcct_el_nodes:
#         if ~('d' in all_nodes_dict[node]):
#             return False
#         else:
#             if all_nodes_dict[node]['d'] < 1.0:
#                 return False
        
        
#     return True


def self_similar_active_new(vcct_els_front_dict,vcct_el_dict,all_nodes_dict):
    
    
    #next_nodes = [el_info['a_node_next_p'] for el,el_info in vcct_els_front_dict.items()]
    vcct_els_front_dict_new = {}
    for el, el_info in vcct_els_front_dict.items():
        antena_stat = vcct_el_any_antena_status_bool(el,
                                        vcct_el_dict,all_nodes_dict,option='deq1')
        # el_is_next_node = False
        # if 'a_node_prev' in all_nodes_dict[el].keys():
        #     if all_nodes_dict[all_nodes_dict[el]['a_node_prev']]['d'] != 1.0:
        #         el_is_next_node = True
        
        vcct_els_front_dict[el]['ss_a'] =  antena_stat #and not el_is_next_node
        if vcct_els_front_dict[el]['ss_a']  == True:
            vcct_els_front_dict_new[el] = el_info
        
    return vcct_els_front_dict_new

def self_similar_active(vcct_els_front_dict,vcct_el_dict,all_nodes_dict):
    
    
    #next_nodes = [el_info['a_node_next_p'] for el,el_info in vcct_els_front_dict.items()]
    
    for el, el_info in vcct_els_front_dict.items():
        antena_stat = vcct_el_any_antena_status_bool(el,
                                        vcct_el_dict,all_nodes_dict,option='deq1')
        # el_is_next_node = False
        # if 'a_node_prev' in all_nodes_dict[el].keys():
        #     if all_nodes_dict[all_nodes_dict[el]['a_node_prev']]['d'] != 1.0:
        #         el_is_next_node = True
        
        vcct_els_front_dict[el]['ss_a'] =  antena_stat #and not el_is_next_node
        
    return vcct_els_front_dict


def vcct_el_any_antena_status_bool(vcct_el,vcct_el_dict,all_nodes_dict,option='dgt0'):
    
    a_nodes_indexes = [0,2,4,6]
    vcct_el_nodes = vcct_el_dict[vcct_el]
    status = False
    for a_node in vcct_el_nodes[a_nodes_indexes]:
        if a_node != 0:
            if 'd' in all_nodes_dict[a_node]:
                if option == 'deq1':
                    if all_nodes_dict[a_node]['d'] == 1.0:
                       status = True
                       break
                elif option == 'dgt0':
                    if all_nodes_dict[a_node]['d'] > 0.0:
                        status = True
                        break
    return status

def vcct_els_crackfront(vcct_el_dict,all_nodes_dict,option_node = 'dneq1',option_a_nodes='dgt0'):
    """determines VCCT elements at the crack front: one the antenna nodes
    is partially failed"""
        
    vcct_els_cfront = []
    if option_node == 'dneq1':
        aux_option_node = False
    else:
        aux_option_node = True
        
    for vcct_el in vcct_el_dict.keys():
        if all_nodes_dict[vcct_el]['d'] != 1.0 or aux_option_node:
            if vcct_el_any_antena_status_bool(vcct_el,vcct_el_dict,all_nodes_dict,option_a_nodes) == True:
                vcct_els_cfront += [vcct_el]
            
            
    return vcct_els_cfront

def sel_a_node_based_on_ns2crack_vector(vcct_el,ns2crack,vcct_el_dict,
                                        all_nodes_dict,option):

    a_nodes_indexes = [0,2,4,6]
    no_antena_nodes = 4
    no_of_normals = len(ns2crack[:,0])
    a_node_prist_along_normal = np.array([None]*no_of_normals)
    
    vcct_el_nodes = vcct_el_dict[vcct_el]
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
    
    #if progressive release has started, keep pristine node
    if all_nodes_dict[vcct_el]['d'] > 0.001:
        a_node_prist_along_normal[:] = all_nodes_dict[vcct_el]['a_node_next_p']

    else:
    
        for i,normal in enumerate(ns2crack):
            d_dist = np.zeros(no_antena_nodes)
            a_node_pristine = np.zeros(no_antena_nodes)
            #order the antena nodes closest to the normal considered
            for j,a_node in enumerate(vcct_el_nodes[a_nodes_indexes]):
                if a_node != 0:
                    if all_nodes_dict[a_node]['d'] == 0 or option == 'off':
                        a_node_xy = all_nodes_dict[a_node]['coords'][:2]
                        p_v = a_node_xy - vcct_el_xy
                        d_vect = p_v - np.dot(p_v,normal)*normal
                        d_dist[j] = np.linalg.norm(d_vect)
                        a_node_pristine[j] = a_node
                    
            a_node_pristine = a_node_pristine[np.argsort(d_dist)]

            if option == 'boeing':
                a_node_prist_along_normal[i] = a_node_pristine[0]

            elif option == 'carvalho':
                cf_v = all_nodes_dict[vcct_el]['cf_v']
                c_face1 = line2D()
                c_face1.n = cf_v[0,:]
                c_face1.point = vcct_el_xy

                c_face2 = line2D()
                c_face2.n = cf_v[1,:]
                c_face2.point = vcct_el_xy

                for a_node in a_node_pristine:
                    #check if the vector starting at a_node with direction normal
                    #intersects the crack front
                    if a_node!=0:
                        a_n_v = line2D() 
                        a_n_v.point = all_nodes_dict[a_node]['coords'][:2]     
                        a_n_v.n = normal

                        int_l_1,t_1,col_1 = c_face1.intersect_w_line(a_n_v)
                        int_l_2,t_2,col_2 = c_face2.intersect_w_line(a_n_v)
                        if t_1 != None:
                            if (t_1 >= 0 and t_1 < 1.0 and col_1 == False): 
                                a_node_prist_along_normal[i] = a_node
                                break
                        if t_2 != None:
                             if (t_2>= 0 and t_2 < 1.0 and col_2 == False): 
                                a_node_prist_along_normal[i] = a_node
                                break
    
    return a_node_prist_along_normal

def trim_based_on_a_node_prist(a_node_prist_along_normal,ns2crack):
    
    ns2crack = np.array([ns2crack[i] for i,a_node in enumerate(a_node_prist_along_normal) if a_node != None])
    a_node_prist_along_normal = np.array([a_node for a_node in a_node_prist_along_normal if a_node != None])
        
    return ns2crack,a_node_prist_along_normal


def antenna_node_max_d(vcct_el,vcct_el_dict,all_nodes_dict):
    """returns the antenna node no and its index in the element being considered
    with maximum d. 
    If multiple antenna are failed, returns the first detected"""

    dstat_ref = 0.0
    a_nodes_indexes = [0,2,4,6]
    ind_a_node = None
    a_node = None
    vcct_el_antena_nodes = vcct_el_dict[vcct_el][a_nodes_indexes]    
    for i,a_node_i in enumerate(vcct_el_antena_nodes):
        if a_node_i != 0:
            dstat_anode = all_nodes_dict[a_node_i]['d']
            #if there is a damaged antenna node
            #damaged nodes are considered in addition to failed to enable extrapolation of the propagation
            #i.e. the calculation of ERRs for nodes that are just about the be crack tip nodes (active)
            if dstat_anode > 0.0:
                if dstat_anode > dstat_ref:
                    #selects the reference antenna node to be the node with the largest d value
                    ind_a_node = a_nodes_indexes[i]  
                    a_node = a_node_i
                    dstat_ref = dstat_anode
                    if dstat_ref == 1.0:
                        break
    
    
    
    return a_node, ind_a_node


def a_nodes_cond(vcct_el,vcct_el_dict,all_nodes_dict,option='deq1'):
    
    """returns a failed antena nodes if any
    NOTE: can be used in multiple subroutines - cleanup"""    
    a_nodes_f = []
    vcct_el_nodes = vcct_el_dict[vcct_el]
    a_nodes_indexes = [0,2,4,6]
    for a_node in vcct_el_nodes[a_nodes_indexes]:
        if a_node != 0:
            if option == 'deq1':
                if  all_nodes_dict[a_node]['d'] == 1:
                    a_nodes_f += [a_node]
            elif option == 'dgt0':
                if  all_nodes_dict[a_node]['d'] > 0.0:
                    a_nodes_f += [a_node]
            elif option == 'deq0':
                if  all_nodes_dict[a_node]['d'] == 0.0:
                    a_nodes_f += [a_node]
                    
    return a_nodes_f

def vcct_el_nodes_pristine(vcct_el,vcct_el_dict,all_nodes_dict):
    """returns a failed antena nodes if any
    NOTE: can be used in multiple subroutines - cleanup"""    
    a_nodes_f = []
    vcct_el_nodes = vcct_el_dict[vcct_el]
    for a_node in vcct_el_nodes:
        if a_node != 0:
            if  all_nodes_dict[a_node]['d'] < 0.001:
                a_nodes_f += [a_node]
                
    return a_nodes_f

def cf_v_upd_dstat(cf_v,vcct_el,vcct_el_dict,all_nodes_dict):
    
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
    vcct_el_nodes = vcct_el_dict[vcct_el]
    for j,cv in enumerate(cf_v):
        c_face1 = line2D()
        c_face1.n = cv
        c_face1.point = vcct_el_xy
        
        for i in range(7):    
            node = vcct_el_nodes[index_roll(i,8)]
            if node != 0:
                node_xy = all_nodes_dict[node]['coords'][:2]
                
                next_node = vcct_el_nodes[index_roll(i+1,8)]
                next_node_xy = all_nodes_dict[next_node]['coords'][:2]
                
                boundary_line = line2D()
                boundary_line.point = node_xy
                boundary_line.n = next_node_xy - node_xy
                
                int_l_1,t_1,col_1 = c_face1.intersect_w_line(boundary_line)
                int_l_2,t_2,col_2 = boundary_line.intersect_w_line(c_face1)
                if t_1 != None:
                    if t_1 > 0.0 and t_2 >= 0.0 and t_2 <= 1.0:
                        cf_v[j,:] = int_l_1 - vcct_el_xy
                        break
                                   
    return cf_v


def expand_cfront_for_ctod_extraction(cfront_nodes,ref_vector,nsearches,vcct_el_dict,all_nodes_dict):

    expanded_cfront = []
    for i,f_node in enumerate(cfront_nodes):
        a_node_ref = antena_node_along_ctod_ref_vector(f_node,ref_vector,vcct_el_dict,all_nodes_dict)
        expanded_cfront += [a_node_ref]
        
        for j in range(1,nsearches+2): 
            node_ctod = search_m_line(a_node_ref,f_node,nsearches,vcct_el_dict)
            expanded_cfront +=[node_ctod]

    expanded_cfront = np.append(cfront_nodes,np.array(expanded_cfront))
    
    expanded_cfront = np.unique(expanded_cfront)
    return expanded_cfront

def ctod_calc(cfront_nodes,ref_vect,nsearches_ctod,vcct_el_dict,all_nodes_dict):
    
    nodes_ctod = np.zeros((len(cfront_nodes),11))
    for i,f_node in enumerate(cfront_nodes):
        a_node_ref = antena_node_along_ctod_ref_vector(f_node,ref_vect,vcct_el_dict,all_nodes_dict)
        
        node_ctod = int(search_m_line(a_node_ref,f_node,nsearches_ctod,vcct_el_dict))
        
        if 'delta' in all_nodes_dict[node_ctod].keys():
            ctod = all_nodes_dict[node_ctod]['delta']
        else:
            ctod = np.array([0,0,0])
            print('no data available for the requested number of ctod searches at node:',node_ctod)
        coords = all_nodes_dict[node_ctod]['coords'] 
        

        nodes_ctod[i,0] = f_node
        nodes_ctod[i,1:4] = all_nodes_dict[f_node]['coords']
        nodes_ctod[i,4] = node_ctod
        nodes_ctod[i,5:8] = coords
        nodes_ctod[i,8:] = ctod
        
    return nodes_ctod


def antena_node_along_ctod_ref_vector(node_ref,ctod_ref_vect,vcct_el_dict,all_nodes_dict):
    
    
    vcct_adj_el = vcct_el_dict[node_ref]
    #antena nodes:
    a_nodes = vcct_adj_el[a_nodes_indexes]
    #coordinates of reference node:
    node_ref_coords = all_nodes_dict[node_ref]['coords']
    v_dot_node_ref_a_nodes = []
    for a_i in a_nodes:
        if a_i != 0:
            coords_i = all_nodes_dict[a_i]['coords']
            v = coords_i-node_ref_coords
            v_dot_node_ref = np.dot(v,ctod_ref_vect)
            v_dot_node_ref_a_nodes +=[v_dot_node_ref]
        
    a_node_index = np.argmax(v_dot_node_ref_a_nodes)    
   
    return  a_nodes[a_node_index]



def search_m_line(node_ref,node_wake,nsearches,vcct_el_dict):
    
    node_front = np.copy(node_ref)
    node_back = np.copy(node_wake)
    for i in range(nsearches):
        next_node = search_next_node_mline(node_front,node_back,vcct_el_dict)
        if next_node:
            node_back= node_front
            node_front = next_node
        else:
            print('next_node',next_node)
            # print('search incomplete in "search_m_line"')
            print('target nsearches',nsearches)
            print('current nsearches',i)
            print('node_ref',node_ref)
            print('node_wake',node_wake)
            break
        
       
    return node_front



def search_m_line_cond(node_ref,node_wake,d_val,nsearches,vcct_el_dict,all_nodes_dict):

    for i in range(nsearches):
        next_node = search_next_node_mline(node_ref,node_wake,vcct_el_dict)
        if next_node :
            if all_nodes_dict[next_node]['d'] != d_val:
                break
            else:
                node_wake = node_ref
                node_ref = next_node
        else:
            break
            
            
    if next_node:        
        if d_val == 1.0:
            node_d_0 = next_node
            node_d_1 =  search_next_node_mline(node_ref,node_wake,vcct_el_dict)
        elif d_val == 0.0:
            node_d_0 = node_ref            
            node_d_1 = next_node
    else:
        if d_val == 1.0:
            node_d_0 = node_ref
            node_d_1 = node_wake
        elif d_val == 0.0:
            node_d_0 = node_ref          
            node_d_1 = node_wake
       
    return node_d_0,node_d_1



def search_next_node_mline(node_ref,node_wake,vcct_el_dict):
    
    a_nodes_indexes = [0,2,4,6]
    vcct_adj_el = vcct_el_dict[int(node_ref)]
    a_nodes = vcct_adj_el[a_nodes_indexes]
    next_node = None
    for i,node in enumerate(a_nodes):
        if node == node_wake:
            next_node = vcct_adj_el[a_nodes_indexes[index_roll(i+2,4)]] 
            if next_node == 0:
                next_node = None
            
    return next_node



def normal_vcct_el(vcct_el,vcct_el_dict,all_nodes_dict):
    
    tol = 1.0E-9
    a_nodes_indexes = [0,2,4,6]
    vcct_adj_els = vcct_el_dict[vcct_el]
    
    vcct_el_coords = all_nodes_dict[vcct_el]['coords']
        
    coords_a_nodes = np.zeros((4,3))
    ind_coords = 0
    for a_node_ind in a_nodes_indexes:
        if vcct_adj_els[a_node_ind] != 0:
            node_ref = vcct_adj_els[a_node_ind] 
            coords_a_nodes[ind_coords,:] = all_nodes_dict[node_ref]['coords']
            ind_coords += 1
        
    coords_a_nodes = coords_a_nodes[:ind_coords]
    no_a_nodes = len(coords_a_nodes[:,0])
    normal_vectors = np.zeros_like(coords_a_nodes)
    count_normal_vectors = 0
    for i,cs_a_n in enumerate(coords_a_nodes):
        n_vector = np.cross(cs_a_n-vcct_el_coords,coords_a_nodes[index_roll(i+1,no_a_nodes),:]-vcct_el_coords)
        n_vector_norm = np.linalg.norm(n_vector)  
        if n_vector_norm > tol:
            normal_vectors[count_normal_vectors,:] = n_vector/n_vector_norm
            count_normal_vectors +=1
    normal_vectors = normal_vectors[:count_normal_vectors,:]
    
    normal_vect_final = np.sum(normal_vectors,0)
    
    return normal_vect_final/np.linalg.norm(normal_vect_final) 


def calc_area_VCCT_el_pristine(vcct_el_coords,vcct_adj_els,all_nodes_dict):
    
    area_vcct_el = 1E-6
    for a_node_ind in a_nodes_indexes:
        a_node_no = vcct_adj_els[a_node_ind]
        if a_node_no != 0:
            a_node_coords = all_nodes_dict[a_node_no]['coords']
            v_a = a_node_coords - vcct_el_coords
            
            next_node_ind =  next_a_node_ind_by_inc(a_nodes_indexes,a_node_ind,1)
            next_node_no = vcct_adj_els[next_node_ind]
            if next_node_no != 0:    
                next_a_node_coords = all_nodes_dict[next_node_no]['coords']
                v_anext = next_a_node_coords - vcct_el_coords
                a_inc = np.linalg.norm(np.cross(v_a,v_anext))/2.0
                area_vcct_el += a_inc
                
    
    return area_vcct_el


def calc_area_VCCT_el_pristine_v2(vcct_el,vcct_el_dict,all_nodes_dict):

    #finds the vectors between the central node in the points half way between the center node and the antenna nodes
    #computes the area of the squares defined by those vectors
    #one square for a corner node, two squares for an edge node and four squares for a 'middle' node
    el_type = vcct_el_type(vcct_el,vcct_el_dict)
    vcct_el_coords = all_nodes_dict[vcct_el]['coords']
    
    a_node_1 = vcct_el_dict[vcct_el][0]
    a_node_2 = vcct_el_dict[vcct_el][2]

    a_node_1_coords = all_nodes_dict[a_node_1]['coords']
    a_node_2_coords = all_nodes_dict[a_node_2]['coords']

    a_node_1_mid = (a_node_1_coords-vcct_el_coords)*0.5
    a_node_2_mid = (a_node_2_coords-vcct_el_coords)*0.5

    #area of a corner node
    area_vcct_el = np.linalg.norm(np.cross(a_node_1_mid,a_node_2_mid))
  
    if ((el_type == 'e') or  (el_type == 'm')):
        a_node_3 = vcct_el_dict[vcct_el][4]
        a_node_3_coords = all_nodes_dict[a_node_3]['coords']
        a_node_3_mid = (a_node_3_coords-vcct_el_coords)*0.5
        area_vcct_el += np.linalg.norm(np.cross(a_node_2_mid,a_node_3_mid))
        
        if el_type == 'm':
              a_node_4 = vcct_el_dict[vcct_el][6]
              a_node_4_coords = all_nodes_dict[a_node_4]['coords']
              a_node_4_mid = (a_node_4_coords-vcct_el_coords)*0.5
              
              area_vcct_el += np.linalg.norm(np.cross(a_node_3_mid,a_node_4_mid))
              
              area_vcct_el += np.linalg.norm(np.cross(a_node_4_mid,a_node_1_mid))

    
    return area_vcct_el


def calc_failed_area(vcct_el,vcct_el_coords,vcct_adj_els,a_node_wake_index,a_node_prist_along_normal,cf_vs,
                                  node_2vcct_adj_els_dict,all_nodes_dict,all_els_dict):
    
    f_area_inc = 0.0
    tol = 1E-4 #based on CELS-skewed mesh
    no_antena_nodes = len(a_nodes_indexes)
    inc_arr = np.array([1,-1])
    for inc in inc_arr:
        a_node_1_ind = a_node_wake_index
        #selects the crack face_clockwise/conter_clockwise
        if inc == 1:
            cf_v_i = cf_vs[0,:]
        else:
            cf_v_i = cf_vs[1,:]
        bool_intersect = False
        crack_pos = vcct_el_coords + cf_v_i
        for i in range(2):
            a_node_1 = vcct_adj_els[a_node_1_ind]
            a_node_2_ind = next_a_node_ind_by_inc(a_nodes_indexes,a_node_1_ind, inc)
            a_node_2 = vcct_adj_els[a_node_2_ind]
            if a_node_2 == 0:
                break
            else:
                a_node_1_coords = all_nodes_dict[a_node_1]['coords']
                a_node_2_coords = all_nodes_dict[a_node_2]['coords']
                
                nodes_arr = np.array([a_node_1,a_node_2])
                adj_el = adj_el_to_vcct_node_from2nodes(vcct_el,nodes_arr,node_2vcct_adj_els_dict,all_els_dict)
                           
                vcct_el_nat_coords = inverse_iso_p(adj_el,vcct_el_coords,all_els_dict,all_nodes_dict)

                normal_2adj_el = calc_normal2el(adj_el,all_els_dict,all_nodes_dict) 
                crack_pos_proj = project_point_on_plane(normal_2adj_el,vcct_el_coords,crack_pos)
                crack_pos_nat_coords = inverse_iso_p(adj_el,crack_pos_proj,all_els_dict,all_nodes_dict)
                v1 = a_node_1_coords - vcct_el_coords
                
                #finds diagonal line between the two antena nodes within element adj_el
                a_node_1_nat_coords = inverse_iso_p(adj_el,a_node_1_coords,all_els_dict,all_nodes_dict)
                a_node_2_nat_coords = inverse_iso_p(adj_el,a_node_2_coords,all_els_dict,all_nodes_dict)
                
                line_diagonal_nat = line2D()
                line_diagonal_nat.line_two_points(a_node_1_nat_coords,a_node_2_nat_coords)
                
                line_cfront_nat = line2D()
                line_cfront_nat.line_two_points(vcct_el_nat_coords,crack_pos_nat_coords)
                
                intersection_nat_coords,t0, c_linear = line_diagonal_nat.intersect_w_line(line_cfront_nat)
                if ((len(intersection_nat_coords) != 0) and (t0 >=0-tol) and (t0<=1.0+tol)):
                    intersection_coords = interpol_el_nodal_val_wnat_coords('coords',intersection_nat_coords,adj_el,3,'quad_linear',all_els_dict,all_nodes_dict)
                    v2 = intersection_coords - vcct_el_coords
                    bool_intersect = True
                else:
                    v2 = a_node_2_coords - vcct_el_coords
                    a_node_1_ind = a_node_2_ind
                    
                f_area_inc += abs(np.linalg.norm(np.cross(v1,v2)))/2.0 
                if bool_intersect == True:
                    break
    
    return f_area_inc

def calc_area_vcct_el_cfront(vcct_el,vcct_el_coords,vcct_adj_els,a_node_wake_index,
                             a_node_prist_along_normal,cf_vs,
                                 node_2vcct_adj_els_dict,all_nodes_dict,all_els_dict):
    
    area_vcct_el = calc_area_VCCT_el_pristine(vcct_el_coords,vcct_adj_els,all_nodes_dict)
    
    area_vcct_failed = calc_failed_area(vcct_el,vcct_el_coords,vcct_adj_els,a_node_wake_index,a_node_prist_along_normal,cf_vs,
                                     node_2vcct_adj_els_dict,all_nodes_dict,all_els_dict)
    
    
    return np.array([area_vcct_el - area_vcct_failed])
       


def cfront_faces(vcct_el,vcct_el_coords,vcct_adj_els,vcct_el_dict,
                          all_nodes_dict,a_node_wake_index,a_front_node,option = 'corner_nodes'):
    
    #counter and clockwise increment [counter_clock,clock] if including corner nodes in the crack position determination
    cf_vs = np.zeros((2,3))
    
    if option == 'corner_nodes':
        #cycles one node at the time
        inc_arr =  np.array([1,-1])
        
    # elif option == 'no_corner_nodes':
    #     #cycles two nodes at the time, ignoring the corner nodes
    #     inc_arr = np.array([2,-2])
    no_cf_vs = 0

    front_node_coords = all_nodes_dict[a_front_node]['coords']
    vcct_el_d = all_nodes_dict[vcct_el]['d']
    vcct_el_coords_actual = vcct_el_coords + (front_node_coords - vcct_el_coords)*vcct_el_d
    for j,inc in enumerate(inc_arr):
        ind_a_node = a_node_wake_index
        for i in range(len(a_nodes_indexes_w_corners)):
            
            #connects vcct_el to the trailing wake node for failed area calculation
            ind_a_node = next_a_node_ind_by_inc(a_nodes_indexes_w_corners,ind_a_node, inc)
            a_node = vcct_adj_els[ind_a_node]
            if a_node == 0:
                break
            else: 
                next_a_node_ind = next_a_node_ind_by_inc(a_nodes_indexes_w_corners,ind_a_node, inc)
                next_a_node = vcct_adj_els[next_a_node_ind]
                d_a_node = all_nodes_dict[a_node]['d']
                if next_a_node == a_front_node:
                    if option == 'corner_nodes':
                        if d_a_node > 0.0:
                            cf_vs[j,:] = all_nodes_dict[a_node]['coords'] - vcct_el_coords_actual
                            break
                            
               
                a_node_coords = all_nodes_dict[a_node]['coords']
                next_a_node_coords = all_nodes_dict[next_a_node]['coords'] 
                if d_a_node <= 0.0:
                    all_nodes_dict[a_node]['d'] = 0.0
                    crack_pos = a_node_coords  
                    cf_vs[j,:] = crack_pos - vcct_el_coords_actual
                    no_cf_vs += 1
                    break
                elif d_a_node > 0.0 and d_a_node <1.0:
                    crack_pos = a_node_coords+(next_a_node_coords-a_node_coords)*d_a_node  
                    cf_vs[j,:]  = crack_pos - vcct_el_coords_actual
                    no_cf_vs += 1
                    break

    return cf_vs


# def area_and_cfront_faces(vcct_el,vcct_el_coords,vcct_adj_els,
#                                  node_2vcct_adj_els,vcct_el_dict,
#                           all_nodes_dict,all_els_dict,a_node_wake_index,a_front_node,option = 'corner_nodes'):
    
#     area_vcct_el_pristine = calc_area_VCCT_el_pristine(vcct_el_coords,vcct_adj_els,all_nodes_dict)
    
#     #counter and clockwise increment [counter_clock,clock] if including corner nodes in the crack position determination
#     cf_v = np.zeros((2,3))
    
#     if option == 'corner_nodes':
#         inc_arr =  np.array([1,-1])
        
#     elif option == 'no_corner_nodes':
#         inc_arr = np.array([2,-2])
#     no_cf_vs = 0
    
#     f_area = 0.0
   
#     front_node_coords = all_nodes_dict[a_front_node]['coords']
#     vcct_el_d = all_nodes_dict[vcct_el]['d']
#     vcct_el_coords_actual = vcct_el_coords + (front_node_coords - vcct_el_coords)*vcct_el_d
#     for inc in inc_arr:
#         ind_a_node = a_node_wake_index
#         a_w_node = vcct_adj_els[ind_a_node]
#         cw_v = all_nodes_dict[a_w_node]['coords'] - vcct_el_coords 
#         for i in range(7):
            
#             #connects vcct_el to the trailing wake node for failed area calculation
#             ind_a_node = next_a_node_ind_by_inc(a_nodes_indexes_w_corners, ind_a_node, inc)
#             a_node = vcct_adj_els[ind_a_node]
#             if a_node == 0:
#                 break
#             else: 
#                 f_area_inc =0.0 
#                 next_a_node_ind = next_a_node_ind_by_inc(a_nodes_indexes_w_corners, ind_a_node, inc)
#                 next_a_node = vcct_adj_els[next_a_node_ind]
#                 d_a_node = all_nodes_dict[a_node]['d']
#                 if next_a_node == a_front_node:
#                     if d_a_node == 1.0:
#                         d_a_node = 0.99999
               
#                 a_node_coords = all_nodes_dict[a_node]['coords']
               
#                 next_a_node_coords = all_nodes_dict[next_a_node]['coords'] 
#                 nodes_arr = np.array([a_node,next_a_node])
#                 if d_a_node <= 0.0:
#                     all_nodes_dict[a_node]['d'] = 0.0
#                     crack_pos = a_node_coords
#                     cf_v_i =  crack_pos - vcct_el_coords_actual
#                     cf_v[no_cf_vs,:] = cf_v_i
#                     no_cf_vs += 1
#                     if d_a_node < 0.0:
#                         print('warning: d registered negative, reset to zero',a_node,d_a_node)
#                     f_area_inc,cw_v = failed_area_inc(vcct_el,vcct_el_coords,cw_v,nodes_arr,cf_v_i,
#                                                   all_nodes_dict,
#                                         node_2vcct_adj_els,all_els_dict)
#                     break
#                 elif d_a_node > 0.0 and d_a_node <1.0:
#                     crack_pos = a_node_coords+(next_a_node_coords-a_node_coords)*d_a_node
#                     cf_v_i =  crack_pos - vcct_el_coords_actual
#                     cf_v[no_cf_vs,:]  = cf_v_i
#                     no_cf_vs += 1
#                     f_area_inc,cw_v = failed_area_inc(vcct_el,vcct_el_coords,cw_v,nodes_arr,cf_v_i,all_nodes_dict,
#                                         node_2vcct_adj_els,all_els_dict)
#                     break
#                 elif d_a_node == 1.0:
#                     crack_pos = a_node_coords
#                     cf_v_i = crack_pos - vcct_el_coords
#                     f_area_inc,cw_v = failed_area_inc(vcct_el,vcct_el_coords,cw_v,nodes_arr,cf_v_i,all_nodes_dict,
#                                         node_2vcct_adj_els,all_els_dict)
#                 if vcct_el == 492:    
#                     print('f_area_inc',f_area_inc)
#                     print('cw_v',cw_v)
#                 f_area += f_area_inc
#     area = area_vcct_el_pristine - f_area

#     return area,cf_v[:no_cf_vs,:]
    




# def cfront_faces(vcct_el,vcct_el_dict,all_nodes_dict,ind_a_f,option = 'carvalho'):
    
#     a_nodes_indexes = [0,2,4,6]
#     ind_side_a_nodes = np.array([a_nodes_indexes[index_roll(ind_a_f + 1,4)],
#                               a_nodes_indexes[index_roll(ind_a_f - 1,4)]])
    
#     cf_v = np.empty((0,3))
    
#     vcct_adj_els = vcct_el_dict[vcct_el]
#     inc_ori = np.array([1,-1])
#     vcct_el_xy = all_nodes_dict[vcct_el]['coords']   
#     vcct_el_dstat = all_nodes_dict[vcct_el]['d']
#     for i,ind in enumerate(ind_side_a_nodes):
#         node_ref = vcct_adj_els[ind]
#         if node_ref != 0:       
#             d_ref = all_nodes_dict[node_ref]['d'] #damage state of the antenna node
#             if d_ref < 0.0:
#                 d_ref = 0.0
#                 all_nodes_dict[node_ref]['d'] = 0.0
#                 print('warning: d registered negative(el,d) - reset to zero:',node_ref,d_ref)
            
#             if d_ref == 0.00:    
#                 if option == 'carvalho':
#                     node_d_0 = vcct_adj_els[index_roll(ind - inc_ori[i],8)]   
#                 elif option == 'carvalho_no_corner':
#                     node_d_0 = node_ref
#                 node_d_1 = node_ref
#             elif d_ref > 0.0:
#                 node_d_0 = node_ref
#                 node_d_1 = vcct_adj_els[index_roll(ind + inc_ori[i],8)]

#             node_d_0_coords = all_nodes_dict[node_d_0]['coords']
#             node_d_1_coords = all_nodes_dict[node_d_1]['coords']
#             cf_v_tip = node_d_0_coords+(node_d_1_coords-node_d_0_coords)*all_nodes_dict[node_d_0]['d']
            
#             vcct_el_next_node = vcct_adj_els[a_nodes_indexes[index_roll(ind_a_f + 2,4)]]
#             if vcct_el_next_node != 0:
#                 vcct_el_next_node_xy = all_nodes_dict[vcct_el_next_node]['coords']
#             else:
#                 vcct_el_next_node_xy = vcct_el_xy
#             cf_v = np.append(cf_v,np.array([cf_v_tip - (vcct_el_xy+(vcct_el_next_node_xy-vcct_el_xy)*vcct_el_dstat)]),axis=0)
            
#     all_nodes_dict[vcct_el]['cf_v'] = cf_v
    
#     return all_nodes_dict




# def cfront_faces(vcct_el,vcct_el_dict,all_nodes_dict,ind_a_f,option = 'carvalho2'):
    
#     nsearches = 10
    
#     el_type = vcct_el_type(vcct_el,vcct_el_dict)
#     a_nodes_indexes = [0,2,4,6]
#     adj_nodes_indexes = [0,1,2,3,4,5,6,7]
#     ind_side_a_nodes = np.array([a_nodes_indexes[index_roll(ind_a_f + 1,4)],
#                              a_nodes_indexes[index_roll(ind_a_f - 1,4)]])
    
#     cf_v = np.empty((0,2))

    
#     vcct_adj_els = vcct_el_dict[vcct_el]
#     inc_ori = np.array([1,-1])
#     vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]   
#     vcct_el_dstat = all_nodes_dict[vcct_el]['d']
#     for i,ind in enumerate(ind_side_a_nodes):
#         node_ref = vcct_adj_els[ind]
#         if node_ref != 0:       
#             d_ref = all_nodes_dict[node_ref]['d'] 
#             if d_ref == 0.00:
#                 #search "back" -> towards crack wake
#                 node_wake = vcct_adj_els[index_roll(ind + inc_ori[i],8)]
#                 node_d_1,node_d_0 = search_m_line_cond(node_ref,node_wake,0.0,
#                                                        nsearches,vcct_el_dict,all_nodes_dict)
         
#             elif d_ref == 1.0:
#                 node_wake = vcct_adj_els[index_roll(ind - inc_ori[i],8)]
#                 #search "forward" -> towards crack front
#                 node_d_0,node_d_1 = search_m_line_cond(node_ref,node_wake,1.0,
#                                                        nsearches,vcct_el_dict,all_nodes_dict)
#             else:
#                 node_d_0 = node_ref
#                 node_d_1 = vcct_adj_els[ind+inc_ori[i]]
      
#             node_d_0_coords = all_nodes_dict[node_d_0]['coords'][:2]
#             node_d_1_coords = all_nodes_dict[node_d_1]['coords'][:2]
#             cf_v_tip = node_d_0_coords+(node_d_1_coords-node_d_0_coords)*all_nodes_dict[node_d_0]['d']
#             if option == 'carvalho':
#                 if vcct_el_dstat > 0.001:
#                     vcct_el_next_node = all_nodes_dict[vcct_el]['a_node_next_p']
#                     vcct_el_next_node_xy = all_nodes_dict[vcct_el_next_node]['coords'][:2]
#                 else:
#                     vcct_el_next_node_xy = vcct_el_xy
#             else:
#                 vcct_el_next_node = vcct_adj_els[a_nodes_indexes[index_roll(ind_a_f + 2,4)]]
#                 if vcct_el_next_node != 0:
#                     vcct_el_next_node_xy = all_nodes_dict[vcct_el_next_node]['coords'][:2]
#                 else:
#                     vcct_el_next_node_xy = vcct_el_xy
#             cf_v = np.append(cf_v,np.array([cf_v_tip - (vcct_el_xy+(vcct_el_next_node_xy-vcct_el_xy)*vcct_el_dstat)]),axis=0)
            
#     all_nodes_dict[vcct_el]['cf_v'] = cf_v
    
#     return all_nodes_dict


# def a_node_cf_v_ind():
  
#         cf_v0_n = cf_v[0,:]/np.linalg.norm(cf_v[0,:])
#         cf_v1_n = cf_v[1,:]/np.linalg.norm(cf_v[1,:])
        
#         ind_a_node_s1 = a_nodes_indexes[index_roll(ind_a_node_ref+1,4)]
#         a_node_s1 = vcct_el_dict[vcct_el][ind_a_node_s1]
#         a_node_s1_xy = all_nodes_dict[a_node_s1]['coords'][:2]
#         ind_a_node_s2 = a_nodes_indexes[index_roll(ind_a_node_ref-1,4)]
#         a_node_s2 = vcct_el_dict[vcct_el][ind_a_node_s2]
#         a_node_s2_xy =  all_nodes_dict[a_node_s2]['coords'][:2]
#         node_ip1_xy =  all_nodes_dict[node_ip1[0]]['coords'][:2]
        
#         vx = node_ip1_xy - vcct_el_xy
#         vx_l =  np.linalg.norm(vx)
#         vx_n = vx/vx_l
#         l1 = a_node_s1_xy - vcct_el_xy
#         if np.dot(cf_v[0,:],l1) > 0:
#             a_node_cf_v0_xy = a_node_s1_xy
#             a_node_cf_v1_xy = a_node_s2_xy
#         else:
#             a_node_cf_v0_xy = a_node_s2_xy
#             a_node_cf_v1_xy = a_node_s1_xy
            
#         return a_node_cf_v_ind

# def cfront_faces(vcct_el,vcct_el_dict,all_nodes_dict,option = 'carvalho'):
#     """Inputs:
#         calculates the crack front faces for each vcct element
#         admits an option parameter
#         option == carvalho, calculates the crack faces
#         option != returns an empty list"""
    
#     if option == 'carvalho':
        
#         #finds a reference node for the cf calculation
#         #if node is a next node the reference node is the previous node
#         if 'a_node_prev' in all_nodes_dict[vcct_el].keys():
#             a_nodes_indexes = [0,2,4,6]
#             vcct_el_antena_nodes = vcct_el_dict[vcct_el][a_nodes_indexes]    
#             for i,node in enumerate(vcct_el_antena_nodes):   
#                 if node == all_nodes_dict[vcct_el]['a_node_prev']:
#                       ind_a_f = a_nodes_indexes[i] 
#         else:
#             a_node_dummy,ind_a_f = antenna_node_max_d(vcct_el,vcct_el_dict,all_nodes_dict)
            
#         vcct_el_dstat = all_nodes_dict[vcct_el]['d']
#         vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]   
        
#         vcct_el_nodes = vcct_el_dict[vcct_el]
        
        
#         el_type = vcct_el_type(vcct_el,vcct_el_dict)
#         nodes_pristine_count = 0
#         if el_type == 'm':
#             for ind,node in enumerate(vcct_el_nodes):
#                 if node != 0:
#                     if all_nodes_dict[node]['d'] != 1.0:
#                         nodes_pristine_count += 1
#                         node_pristine_ind = ind
        
#         if nodes_pristine_count == 1:
#             inc_arr = [-1,1]
#             cf_v = np.zeros((2,2))
#             node = vcct_el_nodes[node_pristine_ind]
#             node_xy = all_nodes_dict[node]['coords'][:2]
#             for j,inc in enumerate(inc_arr):
#                 next_node = vcct_el_nodes[index_roll(node_pristine_ind+inc,8)]
#                 next_node_xy = all_nodes_dict[next_node]['coords'][:2]
#                 cf_v_tip =  node_xy+(next_node_xy-node_xy)*0.5
#                 if vcct_el_dstat > 0.001:
#                     vcct_el_next_node = all_nodes_dict[vcct_el]['a_node_next_p']
#                     vcct_el_next_node_xy = all_nodes_dict[vcct_el_next_node]['coords'][:2]
#                 else:
#                     vcct_el_next_node_xy = vcct_el_xy
                
#                 cf_v[j,:]  = cf_v_tip - (vcct_el_xy+(vcct_el_next_node_xy-vcct_el_xy)*vcct_el_dstat)
            
#         else:   
#             #if a reference antenna node was found, calculates the crack front vectors:          
#             if ind_a_f != None:
#                 inc_arr = [1,-1]
#                 cf_v = np.zeros((2,2))
#                 for j,inc in enumerate(inc_arr):
#                     ind = ind_a_f
#                     for i in range(6):
#                         node = vcct_el_nodes[index_roll(ind+inc,8)]
#                         if node != 0:
#                             node_dstat = all_nodes_dict[node]['d'] 
#                             if node_dstat !=1:
#                                 node_xy = all_nodes_dict[node]['coords'][:2]
#                                 next_node = vcct_el_nodes[index_roll(ind+inc*2,8)]
#                                 if next_node != 0:
#                                     next_node_xy = all_nodes_dict[next_node]['coords'][:2]
#                                     cf_v_tip = node_xy+(next_node_xy-node_xy)*node_dstat
                                
#                                     if vcct_el_dstat > 0.001:
#                                         vcct_el_next_node = all_nodes_dict[vcct_el]['a_node_next_p']
#                                         vcct_el_next_node_xy = all_nodes_dict[vcct_el_next_node]['coords'][:2]
#                                     else:
#                                         vcct_el_next_node_xy = vcct_el_xy
#                                     cf_v[j,:] = cf_v_tip - (vcct_el_xy+(vcct_el_next_node_xy-vcct_el_xy)*vcct_el_dstat)
#                                     break
#                             ind += inc
#                         else:
#                             break
#             else:
#                 print('WARNING: vcct element(node) '+str(vcct_el) +
#                       ' should not have been selected as a front element; wrong solution likely...')

#     else:
#         cf_v = []
    
   
    
#     if vcct_el_dstat > 0.001: 
#         el_type = vcct_el_type(vcct_el,vcct_el_dict)
#         if el_type == 'm':
#             cf_v = cf_v_upd_dstat(cf_v,vcct_el,vcct_el_dict,all_nodes_dict)
          
#     all_nodes_dict[vcct_el]['cf_v'] = cf_v
    
#     return  all_nodes_dict

def build_vcct_el_dict(node_2vcct_adj_els,all_nodes_dict,all_els_dict):
    vcct_el_dict = {}
    for node_i,adj_els2node_i in node_2vcct_adj_els.items():
        no_adj_els= len(adj_els2node_i)
        if no_adj_els == 1: #one adjacent element
            #finds nodes of adjacent element
            nodes_of_adj_el = all_els_dict[adj_els2node_i[0]]['nodes']
            #determines the antenna nodes of the adjacent element (given node_i)
            an_nodes = antenna_nodes_per_el(node_i,nodes_of_adj_el)
            #extracts the coordinates of node_i and the antenna nodes
            node_i_xy = all_nodes_dict[node_i]['coords']
            an_node1_xy =  all_nodes_dict[an_nodes[0]]['coords']
            an_node2_xy =  all_nodes_dict[an_nodes[1]]['coords']
            #determines vectors from the antenna nodes to the 
            va1 = an_node1_xy - node_i_xy
            va2 = an_node2_xy - node_i_xy
            
            ref_n = ref_cclock_el_2d(nodes_of_adj_el,all_nodes_dict)
            vcct_el_dict[node_i] = np.zeros(8,dtype=int)
            if np.dot(np.cross(va1,va2),ref_n) > 0:
                vcct_el_dict[node_i][0] = an_nodes[0] 
            else:
                vcct_el_dict[node_i][0] = an_nodes[1]
            ind = np.nonzero(nodes_of_adj_el == vcct_el_dict[node_i][0])[0]
            vcct_el_dict[node_i][1] = nodes_of_adj_el[index_roll(ind+1,4)]
            vcct_el_dict[node_i][2] = nodes_of_adj_el[index_roll(ind+2,4)]
        else:
                
            #extracts all nodes and antenna nodes adjacent elements
            nodes_of_adj_el = np.zeros((no_adj_els,4))
            an_nodes_el = np.zeros((no_adj_els,2))
            for i in range(no_adj_els):
                nodes_of_adj_el[i,:] = all_els_dict[adj_els2node_i[i]]['nodes']
                an_nodes_el[i,:] = antenna_nodes_per_el(node_i,nodes_of_adj_el[i,:] )
            
            #finds whether the two first elements on the the list have a common antenna node
            common_antena = np.intersect1d(an_nodes_el[0,:],an_nodes_el[1,:])
            if len(common_antena) != 0:
                common_antena = common_antena[0]
                #check next if next element needs to be swapped:
                if no_adj_els == 4:
                    ca = np.intersect1d(an_nodes_el[1,:],an_nodes_el[2,:])
                    if len(ca) ==0:
                        an_nodes_el[[2,3],:] = an_nodes_el[[3,2],:] 
                        nodes_of_adj_el[[2,3],:] =  nodes_of_adj_el[[3,2],:]             
            else:
                common_antena = np.intersect1d(an_nodes_el[0,:],an_nodes_el[2,:])[0]
                an_nodes_el[[1,2],:] = an_nodes_el[[2,1],:] 
                nodes_of_adj_el[[1,2],:] =  nodes_of_adj_el[[2,1],:] 
           
            common_antena_xy = all_nodes_dict[common_antena]['coords']
            
            el_unique_an = np.zeros(no_adj_els)    
            for i in range(no_adj_els):
                el_unique_an[i] = an_nodes_el[i,an_nodes_el[i,:] != common_antena][0]

            el_0_unique_an_xy = all_nodes_dict[el_unique_an[0]]['coords']
            ref_n = ref_cclock_el_2d(nodes_of_adj_el[0,:],all_nodes_dict)
            
            node_i_xy = all_nodes_dict[node_i]['coords']
            va_el0 = el_0_unique_an_xy - node_i_xy
            va_common = common_antena_xy - node_i_xy
            
            vcct_el_dict[node_i] = np.zeros(9,dtype=int)
            if np.dot(np.cross(va_el0,va_common),ref_n) < 0:
                nodes_of_adj_el[[0,1],:] = nodes_of_adj_el[[1,0],:]
                el_unique_an[[0,1]] = el_unique_an[[1,0]]
                an_nodes_el[[0,1],:] = an_nodes_el[[1,0],:]
                if no_adj_els == 4:
                    common_antena = np.intersect1d(an_nodes_el[1,:],an_nodes_el[2,:])
                    if len(common_antena) == 0:
                        nodes_of_adj_el[[2,3],:] =  nodes_of_adj_el[[3,2],:] 
                
            vcct_el_dict[node_i][0] = el_unique_an[0]
            nodes_added_per_el = 2
            for i in range(no_adj_els):    
                ind = np.nonzero(nodes_of_adj_el[i,:] == vcct_el_dict[node_i][i*nodes_added_per_el])[0][0]
                vcct_el_dict[node_i][i*nodes_added_per_el+1] = nodes_of_adj_el[i,index_roll(ind+1,4)]
                vcct_el_dict[node_i][i*nodes_added_per_el+2] = nodes_of_adj_el[i,index_roll(ind+2,4)]
            vcct_el_dict[node_i] = vcct_el_dict[node_i][:8]
    
    
    return vcct_el_dict

def adj_els_with_sample_points(vcct_el,sample_points_xy,node_2vcct_adj_els_dict,all_nodes_dict,all_els_dict):  
    """determines the adjacent elements that contain sample points
    Inputs: sample points; 
    Outputs: indexes of the sample points inside an adjacent element;
    adjacent els with a sample poing
    """

    adj_els = node_2vcct_adj_els_dict[vcct_el]
    no_sample_points_xy = len(sample_points_xy)
    a_els_with_s_points = np.array([None]*no_sample_points_xy)
    for i,point_xy in enumerate(sample_points_xy):
        for a_el in adj_els:
            if point_in_element_asc(a_el,point_xy,all_els_dict,all_nodes_dict):
                a_els_with_s_points[i] = a_el
                break
    
    ind_inside_el = a_els_with_s_points != None
    
    sample_points_xy = sample_points_xy[ind_inside_el]
    a_els_with_s_points = a_els_with_s_points[ind_inside_el]
    
    return ind_inside_el,a_els_with_s_points
    
    
def load_el_I_2nodes(el_force_I,el_force_onset_I,el_delta_I,all_nodes_dict,all_els_dict):
    for j,el in enumerate(el_force_I):
        for i,node in enumerate(all_els_dict[int(el[0])]['nodes']):
            #shear forces always assumed positive (prevents issues with local in-plane rotaitons of the the elements not being tracked)
            force = el[1+i*3:i*3+4]
            force[0:2] = abs(force[0:2])
            if 'force' in  all_nodes_dict[node]:
                all_nodes_dict[node]['force'] = all_nodes_dict[node]['force'] + force
            else:
               all_nodes_dict[node]['force'] = force
                
            all_nodes_dict[node]['delta'] = el_delta_I[j,1+i*3:i*3+4]
            
    onset_node_list  = []   
    for j,el in enumerate(el_force_onset_I):
        for i,node in enumerate(all_els_dict[int(el[0])]['nodes']):
            #shear forces always assumed positive (prevents issues with local in-plane rotations of the the elements not being tracked)
            force = el[1+i*3:i*3+4]
            force[0:2] = abs(force[0:2])
            onset_node_list += [node]
            if 'force_onset' in  all_nodes_dict[node]:
                all_nodes_dict[node]['force_onset'] = all_nodes_dict[node]['force_onset'] + force
            else:
                all_nodes_dict[node]['force_onset'] = force
    onset_node_list = np.unique(onset_node_list)      
    return all_nodes_dict, onset_node_list
    

def min_distance_to_vertices_nat_space(sp_point):
    
    nat_space_corners = np.array([[-1,-1],
                                  [-1,1],
                                  [1,-1],
                                  [1,1]])
    dist = np.zeros(len(nat_space_corners[:,0]))
    for i,corner in enumerate(nat_space_corners):
        dist[i] = np.linalg.norm(sp_point-corner)
        
    return min(dist)


def eval_sample_point(sp_point,quantity,reg_adj_els,all_els_dict,
                              all_nodes_dict):
    #options for quantity = ['delta','force']
    quant_sp = np.zeros_like(sp_point)
    no_adj_reg_els = len(reg_adj_els)
    for i,sp_i in enumerate(sp_point):
        tol = 1E-6
        sp_point_nat = np.zeros((no_adj_reg_els,2))
        sp_dist_arr = np.zeros(no_adj_reg_els)
        sp_in_el = False
        for j,el in enumerate(reg_adj_els):
            sp_point_nat[j,:] = inverse_iso_p(el,sp_i,all_els_dict,all_nodes_dict)
            if ~(np.any(sp_point_nat[j,:] > 1.0 + tol) or np.any(sp_point_nat[j,:] < -1.0 - tol)):
                sp_in_el = True
                adj_el_ind = j
                break
            else:
                sp_dist_arr[j] = min_distance_to_vertices_nat_space(sp_point_nat[j,:])
        if sp_in_el == False:
            adj_el_ind = np.argmin(sp_dist_arr)
    
        quant_sp[i,:] = interpol_el_nodal_val_wnat_coords(quantity,
                                                    sp_point_nat[adj_el_ind,:],reg_adj_els[adj_el_ind],3,'quad_linear',
                                                    all_els_dict,all_nodes_dict)

    return quant_sp

# def calc_delta_sample_points(sample_points_xy,a_els_with_s_points,all_els_dict,
#                               all_nodes_dict):
    
#     no_sample_points = len(sample_points_xy[:,0])
#     delta_sample_points = np.zeros((no_sample_points,3))
#     for i,sp in enumerate(sample_points_xy):
#         nat = inverse_iso_p(a_els_with_s_points[i],sp,all_els_dict,all_nodes_dict)
#         delta_sample_points[i,:] = interpol_el_nodal_val_wnat_coords('delta',
#                                                     nat,a_els_with_s_points[i],3,'quad_linear',
#                                                     all_els_dict,all_nodes_dict)
 
#         # delta_nodes = np.zeros((4,3))
#         # for j,node in enumerate(all_els_dict[a_els_with_s_points[i]]['nodes']):
#         #     if 'delta' in all_nodes_dict[node]:
#         #         delta_nodes[j,:] = all_nodes_dict[node]['delta']
#         # n_shape = shape_f(nat,'quad_linear')
#         # delta_sample_points[i,:] = np.array([np.dot(n_shape,delta_nodes[:,0]),
#         #                                      np.dot(n_shape,delta_nodes[:,1]),
#         #                                      np.dot(n_shape,delta_nodes[:,2])])

#     return delta_sample_points


# def calc_delta_sample_points(sample_points_xy,a_els_with_s_points,all_els_dict,
#                              all_nodes_dict):
    
#     no_sample_points = len(sample_points_xy[:,0])
#     delta_sample_points = np.zeros((no_sample_points,3))
#     for i,sp in enumerate(sample_points_xy):
#         nat = inverse_iso_p(a_els_with_s_points[i],sp,all_els_dict,all_nodes_dict)
#         delta_sample_points[i,:] = interpol_el_nodal_val_wnat_coords('delta',
#                                                     nat,a_els_with_s_points[i],3,'quad_linear',
#                                                     all_els_dict,all_nodes_dict)
 
#         # delta_nodes = np.zeros((4,3))
#         # for j,node in enumerate(all_els_dict[a_els_with_s_points[i]]['nodes']):
#         #     if 'delta' in all_nodes_dict[node]:
#         #         delta_nodes[j,:] = all_nodes_dict[node]['delta']
#         # n_shape = shape_f(nat,'quad_linear')
#         # delta_sample_points[i,:] = np.array([np.dot(n_shape,delta_nodes[:,0]),
#         #                                      np.dot(n_shape,delta_nodes[:,1]),
#         #                                      np.dot(n_shape,delta_nodes[:,2])])

#     return delta_sample_points


def calc_delta_sample_points_asc(sample_points_xy,a_els_with_s_points,all_els_dict,
                             all_nodes_dict):
    
    no_sample_points = len(sample_points_xy[:,0])
    delta_sample_points = np.zeros((no_sample_points,3))
    for i,sp in enumerate(sample_points_xy):
        nat = inverse_iso_p_asc(a_els_with_s_points[i],sp,all_els_dict,all_nodes_dict)
        delta_sample_points[i,:] = interpol_el_nodal_val_wnat_coords('delta',
                                                    nat,a_els_with_s_points[i],3,'quad_linear',
                                                    all_els_dict,all_nodes_dict)
 
        # delta_nodes = np.zeros((4,3))
        # for j,node in enumerate(all_els_dict[a_els_with_s_points[i]]['nodes']):
        #     if 'delta' in all_nodes_dict[node]:
        #         delta_nodes[j,:] = all_nodes_dict[node]['delta']
        # n_shape = shape_f(nat,'quad_linear')
        # delta_sample_points[i,:] = np.array([np.dot(n_shape,delta_nodes[:,0]),
        #                                      np.dot(n_shape,delta_nodes[:,1]),
        #                                      np.dot(n_shape,delta_nodes[:,2])])

    return delta_sample_points

def calc_force_sample_points_asc(sample_points_xy,a_node_prist_along_normal,vcct_el,
                             node_2vcct_adj_els_dict,all_els_dict,
                             all_nodes_dict,option='next_a_node'):
    
    no_sample_points = len(sample_points_xy[:,0])
    force_sample_points = np.zeros((no_sample_points,3))
    
    if option == 'next_a_node':
        for i,node in enumerate(a_node_prist_along_normal):
            if node == 0:
                force_sample_points[i,:]  = 0
            else:
                if 'force' in all_nodes_dict[node].keys():
                    force_sample_points[i,:] = all_nodes_dict[node]['force']
    
    elif option == 'sample':
        vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
        sample_points_xy_f = np.array([2.0*vcct_el_xy - sp for sp in sample_points_xy])
        ind_inside_el,a_els_with_s_points_f = adj_els_with_sample_points(vcct_el,
                                            sample_points_xy_f,node_2vcct_adj_els_dict,all_nodes_dict,all_els_dict)
        for i,sp in enumerate(sample_points_xy_f):

            nat = inverse_iso_p_asc(a_els_with_s_points_f[i],sp,all_els_dict,all_nodes_dict)
            
            force_sample_points[i,:] = interpol_el_nodal_val_wnat_coords('force',
                                                        nat,a_els_with_s_points_f[i],3,'quad_linear',
                                                        all_els_dict,all_nodes_dict)
    return force_sample_points



def calc_area_carvalho7(vcct_el,ind_ref_adj_el,vcct_el_dict,
                            all_nodes_dict,all_els_dict):
    
    
    # if vcct_el == 37396:
    #     print('here')
    #     fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9, 9))
    #     axes.set_xlim([49.0,50.5])
    #     axes.set_ylim([9.3,10.5])
    #     plot_elements_dict(axes,all_els_dict,all_nodes_dict,'g')
        
    el_type = vcct_el_type(vcct_el,vcct_el_dict)
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
    
   
    
    if el_type == 'e':
        area_p = np.zeros((3,2))
        vcct_el_nodes = vcct_el_dict[vcct_el]
        inc_dir = 1
        if vcct_el_nodes[index_roll(ind_ref_adj_el+1,8)] == 0:
            inc_dir = -1

        area_p[0] = all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p[1] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+inc_dir*2,5)]]['coords'][:2]
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4*inc_dir,8)]
        area_p[2] = all_nodes_dict[next_node]['coords'][:2]
        
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        
        c_face1 = line2D()
        c_face1.n = cf_v[0]
        c_face1.point = vcct_el_xy
        c_faces_lines = [c_face1]
        
        #plot_vect(vcct_el_xy, cf_v[0], axes)
    else:
            
        area_p = np.zeros((4,2))
        vcct_el_nodes = vcct_el_dict[vcct_el]

        #points defining the area1 associated with the node
        area_p[0] =  all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p[1] =  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+2,8)]]['coords'][:2]
        area_p[2] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]]['coords'][:2]
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]
        area_p[3] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+6,8)]]['coords'][:2]
    
       
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        c_face1 = line2D()
        c_face1.n = cf_v[0,:]
        c_face1.point = vcct_el_xy
        
        c_face2 = line2D()
        c_face2.n = cf_v[1,:]
        c_face2.point = vcct_el_xy
        c_faces_lines = [c_face1,c_face2]
        

    start_area_calc = False
    end_area_calc = False
    area = 0
    no_points = len(area_p)
    for i,pi in enumerate(area_p):
        #line from pi to pi+1
       
        line_pi_pip1 = line2D()
        line_pi_pip1.n =  area_p[index_roll(i+1,no_points)] -  pi
        line_pi_pip1.point = pi
        area_inc = abs(np.cross(pi-vcct_el_xy,area_p[index_roll(i+1,no_points)]-vcct_el_xy))/2.0
      
        for j,c_line in enumerate(c_faces_lines):
            int_l,t,col = line_pi_pip1.intersect_w_line(c_line)
        
                
            if len(int_l) > 0:
                int_l,t_aux,col = c_line.intersect_w_line(line_pi_pip1)
                
                if ((t >= 0 and t <= 1.0) and (t_aux >= 0 and t_aux <= 1.0)):  
                    
                    if (((i==0 or i==2) and t < 0.5) or ((i == 1 or i==3) and t > 0.5)):
                            int_l = pi + 0.5*line_pi_pip1.n
                    #print('t',t)
                    #plot_point(int_l,axes,'b')
                    if start_area_calc == False:
                        start_area_calc = True
                        area_inc = abs(np.cross(int_l-vcct_el_xy,area_p[index_roll(i+1,no_points)]-vcct_el_xy))/2.0
                        c_faces_lines.pop(j)
                        break
                    if start_area_calc == True:
                        end_area_calc = True
                        area_inc = abs(np.cross(area_p[i]-vcct_el_xy,int_l-vcct_el_xy))/2.0
                    break
        if start_area_calc == True:
            area += area_inc
        if end_area_calc == True:
            break   
  
    area = (1.0-all_nodes_dict[next_node]['d'])*area + 1.0E-6
    
    return area



def calc_area_carvalho6(vcct_el,ind_ref_adj_el,vcct_el_dict,
                            all_nodes_dict,all_els_dict):

        
    el_type = vcct_el_type(vcct_el,vcct_el_dict)
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
    
    
    if el_type == 'e':
        area_p = np.zeros((3,2))
        vcct_el_nodes = vcct_el_dict[vcct_el]
        inc_dir = 1
        if vcct_el_nodes[index_roll(ind_ref_adj_el+1,8)] == 0:
            inc_dir = -1
        #print('inc_dir',inc_dir)
        area_p[0] = all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p[1] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+inc_dir*2,5)]]['coords'][:2]
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4*inc_dir,8)]
        area_p[2] = all_nodes_dict[next_node]['coords'][:2]
       
        cf_v = all_nodes_dict[vcct_el]['cf_v']
       
        c_face1 = line2D()
        c_face1.n = cf_v[0]
        c_face1.point = vcct_el_xy
        c_faces_lines = [c_face1]
        #plot_vect(vcct_el_xy, cf_v[0], axes)
        #plot_vect(vcct_el_xy, -cf_v[0], axes)
    else:
            
        area_p = np.zeros((4,2))
        vcct_el_nodes = vcct_el_dict[vcct_el]

        #points defining the area1 associated with the node
        area_p[0] =  all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p[1] =  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+2,8)]]['coords'][:2]
        area_p[2] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]]['coords'][:2]
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]
        area_p[3] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+6,8)]]['coords'][:2]
      
        vcct_el_xy = vcct_el_xy
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        c_face1 = line2D()
        c_face1.n = cf_v[0,:]
        c_face1.point = vcct_el_xy
        
        c_face2 = line2D()
        c_face2.n = cf_v[1,:]
        c_face2.point = vcct_el_xy
        c_faces_lines = [c_face1,c_face2]
        
        #plot_vect(vcct_el_xy, cf_v[0], axes)
        #plot_vect(vcct_el_xy, cf_v[1], axes)
        
        # if vcct_el == 37396:
        
        #plot_crack_front_el(cf_v, vcct_el_xy, axes)
   
   
    #for p in area_p:
    #    plot_point(p,axes,'r')
    start_area_calc = False
    end_area_calc = False
    area = 0
    no_points = len(area_p)
    for i,pi in enumerate(area_p):
        #line from pi to pi+1
       
        line_pi_pip1 = line2D()
        line_pi_pip1.n =  area_p[index_roll(i+1,no_points)] -  pi
        line_pi_pip1.point = pi
        area_inc = abs(np.cross(pi-vcct_el_xy,area_p[index_roll(i+1,no_points)]-vcct_el_xy))/2.0
      
        for j,c_line in enumerate(c_faces_lines):
            int_l,t,col = line_pi_pip1.intersect_w_line(c_line)
            
            if len(int_l) > 0:
                int_l,t_aux,col = c_line.intersect_w_line(line_pi_pip1)
                
                if ((t >= 0 and t <= 1.0) and (t_aux >= 0 and t_aux <= 1.0)):  
                    #plot_point(int_l,axes,'c')
                    if start_area_calc == False:
                        start_area_calc = True
                        area_inc = abs(np.cross(int_l-vcct_el_xy,area_p[index_roll(i+1,no_points)]-vcct_el_xy))/2.0
                        c_faces_lines.pop(j)
                        break
                    if start_area_calc == True:
                        end_area_calc = True
                        area_inc = abs(np.cross(area_p[i]-vcct_el_xy,int_l-vcct_el_xy))/2.0
                    break
        if start_area_calc == True:
            area += area_inc
        if end_area_calc == True:
            break
   # print('ratio',area_i[0]/area_i[1])
    #print(all_nodes_dict[next_node]['d'])
    area = (1.0-all_nodes_dict[next_node]['d'])*area + 1.0E-6
    
    return area

def calc_area_carvalho5(vcct_el,ind_ref_adj_el,vcct_el_dict,
                            all_nodes_dict,all_els_dict):

        
    el_type = vcct_el_type(vcct_el,vcct_el_dict)
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
    
    area_p = np.zeros((4,2))
    if el_type == 'e':
        
        vcct_el_nodes = vcct_el_dict[vcct_el]
        inc_dir = 1
        if vcct_el_nodes[index_roll(ind_ref_adj_el+1,8)] == 0:
            inc_dir = -1
        #print('inc_dir',inc_dir)
        area_p[0] = all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p[1] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+inc_dir*2,5)]]['coords'][:2]
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4*inc_dir,8)]
        area_p[2] = all_nodes_dict[next_node]['coords'][:2]
        area_p[3] = -area_p[1]+2.0*vcct_el_xy
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        
        c_face1 = line2D()
        c_face1.n = cf_v[0]
        c_face1.point = vcct_el_xy
        c_face2 = line2D()
        c_face2.n = -cf_v[0]
        c_face2.point = vcct_el_xy
        c_faces_lines = [c_face1,c_face2]
        #plot_vect(vcct_el_xy, cf_v[0], axes)
        #plot_vect(vcct_el_xy, -cf_v[0], axes)
    else:
            
        area_p = np.zeros((4,2))
        vcct_el_nodes = vcct_el_dict[vcct_el]

        #points defining the area1 associated with the node
        area_p[0] =  all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p[1] =  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+2,8)]]['coords'][:2]
        area_p[2] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]]['coords'][:2]
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]
        area_p[3] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+6,8)]]['coords'][:2]
      
        vcct_el_xy = vcct_el_xy
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        c_face1 = line2D()
        c_face1.n = cf_v[0,:]
        c_face1.point = vcct_el_xy
        
        c_face2 = line2D()
        c_face2.n = cf_v[1,:]
        c_face2.point = vcct_el_xy
        c_faces_lines = [c_face1,c_face2]
        
        #plot_vect(vcct_el_xy, cf_v[0], axes)
        #plot_vect(vcct_el_xy, cf_v[1], axes)
        
        # if vcct_el == 37396:
        
        #plot_crack_front_el(cf_v, vcct_el_xy, axes)
   
   
    #for p in area_p:
    #    plot_point(p,axes,'r')
    start_area_calc = False
    end_area_calc = False
    area = 0
    for i,pi in enumerate(area_p):
        #line from pi to pi+1
       
        line_pi_pip1 = line2D()
        line_pi_pip1.n =  area_p[index_roll(i+1,4)] -  pi
        line_pi_pip1.point = pi
        area_inc = abs(np.cross(pi-vcct_el_xy,area_p[index_roll(i+1,4)]-vcct_el_xy))/2.0
      
        for j,c_line in enumerate(c_faces_lines):
            int_l,t,col = line_pi_pip1.intersect_w_line(c_line)
            
            if len(int_l) > 0:
                int_l,t_aux,col = c_line.intersect_w_line(line_pi_pip1)
                
                if ((t >= 0 and t <= 1.0) and (t_aux >= 0 and t_aux <= 1.0)):  
                    #plot_point(int_l,axes,'c')
                    if start_area_calc == False:
                        start_area_calc = True
                        area_inc = abs(np.cross(int_l-vcct_el_xy,area_p[index_roll(i+1,4)]-vcct_el_xy))/2.0
                        c_faces_lines.pop(j)
                        break
                    if start_area_calc == True:
                        end_area_calc = True
                        area_inc = abs(np.cross(area_p[i]-vcct_el_xy,int_l-vcct_el_xy))/2.0
                    break
        if start_area_calc == True:
            area += area_inc
        if end_area_calc == True:
            break
    if el_type == 'e':
        area = area/2.0
   # print('ratio',area_i[0]/area_i[1])
    #print(all_nodes_dict[next_node]['d'])
    #area = (1.0-all_nodes_dict[next_node]['d'])*area + 1.0E-6
    
    return area


def calc_area_carvalho4(vcct_el,ind_ref_adj_el,vcct_el_dict,
                            all_nodes_dict,all_els_dict):
    
    
    # if vcct_el == 37396:
    #     print('here')
    #     fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9, 9))
    #     axes.set_xlim([49.0,50.5])
    #     axes.set_ylim([9.3,10.5])
    #     plot_elements_dict(axes,all_els_dict,all_nodes_dict,'g')
        
    el_type = vcct_el_type(vcct_el,vcct_el_dict)
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
    
   
    
    
    if el_type == 'e':
        area_p = np.zeros((3,2))
        vcct_el_nodes = vcct_el_dict[vcct_el]
        inc_dir = 1
        if vcct_el_nodes[index_roll(ind_ref_adj_el+1,8)] == 0:
            inc_dir = -1

        area_p[0] = all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p[1] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+inc_dir*2,5)]]['coords'][:2]
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4*inc_dir,8)]
        area_p[2] = all_nodes_dict[next_node]['coords'][:2]
        
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        
        c_face1 = line2D()
        c_face1.n = cf_v[0]
        c_face1.point = vcct_el_xy
        c_faces_lines = [c_face1]
    else:
            
        area_p = np.zeros((4,2))
        vcct_el_nodes = vcct_el_dict[vcct_el]

        #points defining the area1 associated with the node
        area_p[0] =  all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p[1] =  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+2,8)]]['coords'][:2]
        area_p[2] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]]['coords'][:2]
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]
        area_p[3] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+6,8)]]['coords'][:2]
    
       
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        c_face1 = line2D()
        c_face1.n = cf_v[0,:]
        c_face1.point = vcct_el_xy
        
        c_face2 = line2D()
        c_face2.n = cf_v[1,:]
        c_face2.point = vcct_el_xy
        c_faces_lines = [c_face1,c_face2]
        
        # if vcct_el == 37396:
        
        #plot_crack_front_el(cf_v, vcct_el_xy, axes)
   
   
    #for p in area_p:
    #    plot_point(p,axes,'r')
    start_area_calc = False
    end_area_calc = False
    area = 0
    for i,pi in enumerate(area_p[:-1]):
        line_pi_pip1 = line2D()
        line_pi_pip1.n =  area_p[i+1] - pi
        line_pi_pip1.point = pi
        area_inc = abs(np.cross(pi-vcct_el_xy,area_p[i+1]-vcct_el_xy))/2.0
        
        for j,c_line in enumerate(c_faces_lines):
            int_l,t,col = line_pi_pip1.intersect_w_line(c_line)
            
            if len(int_l) > 0:
                int_l,t_aux,col = c_line.intersect_w_line(line_pi_pip1)
                
                if ((t >= 0 and t <= 1.0) and (t_aux >= 0 and t_aux <= 1.0)):  
                    #plot_point(int_l,axes,'c')
                    if start_area_calc == False:
                        start_area_calc = True
                        area_inc = abs(np.cross(int_l-vcct_el_xy,area_p[i+1]-vcct_el_xy))/2.0
                        c_faces_lines.pop(j)
                        break
                    if start_area_calc == True:
                        end_area_calc = True
                        area_inc = abs(np.cross(area_p[i]-vcct_el_xy,int_l-vcct_el_xy))/2.0
                    break
        if start_area_calc == True:
            area += area_inc
        if end_area_calc == True:
            break

   # print('ratio',area_i[0]/area_i[1])
    #print(all_nodes_dict[next_node]['d'])
    #area = (1.0-all_nodes_dict[next_node]['d'])*area + 1.0E-6
    
    return area



def calc_area_carvalho3(vcct_el,ind_ref_adj_el,vcct_el_dict,
                            all_nodes_dict,all_els_dict):
    
    
    # if vcct_el == 37396:
    #     print('here')
    #     fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9, 9))
    #     axes.set_xlim([49.0,50.5])
    #     axes.set_ylim([9.3,10.5])
    #     plot_elements_dict(axes,all_els_dict,all_nodes_dict,'g')
        
    el_type = vcct_el_type(vcct_el,vcct_el_dict)
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
    
   
    
    
    if el_type == 'e':
        area_p = np.zeros((5,2))
        vcct_el_nodes = vcct_el_dict[vcct_el]
        inc_dir = 1
        if vcct_el_nodes[index_roll(ind_ref_adj_el+1,8)] == 0:
            inc_dir = -1

        area_p[0] = all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p[1] = all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]*0.5+all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+inc_dir,8)]]['coords'][:2]*0.5
        area_p[2] = vcct_el_xy*0.5 + all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+inc_dir*2,5)]]['coords'][:2]*0.5
        area_p[3] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+3*inc_dir,8)]]['coords'][:2]*0.5 +  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+4*inc_dir,8)]]['coords'][:2]*0.5
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4*inc_dir,8)]
        area_p[4] = all_nodes_dict[next_node]['coords'][:2]
        
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        
        c_face1 = line2D()
        c_face1.n = cf_v[0]
        c_face1.point = vcct_el_xy
        c_faces_lines = [c_face1]
        area_p_list = [area_p]
    else:
            
        area_p1 = np.zeros((8,2))
        area_p2 = np.zeros((8,2))
        vcct_el_nodes = vcct_el_dict[vcct_el]

        #points defining the area1 associated with the node
        area_p1[0] =  all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p1[1] = all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]*0.5+all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+1,8)]]['coords'][:2]*0.5
        area_p1[2] = vcct_el_xy*0.5 + all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+2,8)]]['coords'][:2]*0.5
        area_p1[3] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+3,8)]]['coords'][:2]*0.5 +  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]]['coords'][:2]*0.5
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]
        print('next_node',next_node)
        area_p1[4] = all_nodes_dict[next_node]['coords'][:2]
        area_p1[5] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]]['coords'][:2]*0.5 +  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+5,8)]]['coords'][:2]*0.5
        area_p1[6] = vcct_el_xy*0.5 + all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+6,8)]]['coords'][:2]*0.5
        area_p1[7] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+7,8)]]['coords'][:2]*0.5 + all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]*0.5
    
        #points defining the area2 associated with the node
        area_p2[0] = (all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2] + vcct_el_xy)*0.5
        area_p2[1] = (all_nodes_dict[vcct_el_nodes[ind_ref_adj_el+1]]['coords'][:2]+all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+2,8)]]['coords'][:2])*0.5
        area_p2[2] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+2,8)]]['coords'][:2]
        area_p2[3] = (all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+2,8)]]['coords'][:2] + all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+3,8)]]['coords'][:2])*0.5
        area_p2[4] = (all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]]['coords'][:2] + vcct_el_xy)*0.5
        area_p2[5] = (all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+5,8)]]['coords'][:2] +  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+6,8)]]['coords'][:2])*0.5
        area_p2[6] =  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+6,8)]]['coords'][:2]
        area_p2[7] = (all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+6,8)]]['coords'][:2] + all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+7,8)]]['coords'][:2])*0.5
      
        area_p_list = [area_p1,area_p2]
       
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        c_face1 = line2D()
        c_face1.n = cf_v[0,:]
        c_face1.point = vcct_el_xy
        
        c_face2 = line2D()
        c_face2.n = cf_v[1,:]
        c_face2.point = vcct_el_xy
        c_faces_lines = [c_face1,c_face2]
        # if vcct_el == 37396:
        
        #plot_crack_front_el(cf_v, vcct_el_xy, axes)
    area_i = []
    for area_p in area_p_list: 
        
        #for p in area_p:
        #    plot_point(p,axes,'r')
        if el_type == 'e':
            c_faces_lines = [c_face1]
        else:
            c_faces_lines = [c_face1,c_face2]
        start_area_calc = False
        end_area_calc = False
        area = 0
        for i,pi in enumerate(area_p[:-1]):
            line_pi_pip1 = line2D()
            line_pi_pip1.n =  area_p[i+1] - pi
            line_pi_pip1.point = pi
            area_inc = abs(np.cross(pi-vcct_el_xy,area_p[i+1]-vcct_el_xy))/2.0
            
            for j,c_line in enumerate(c_faces_lines):
                int_l,t,col = line_pi_pip1.intersect_w_line(c_line)
                
                if len(int_l) > 0:
                    int_l,t_aux,col = c_line.intersect_w_line(line_pi_pip1)
                    
                    if ((t >= 0 and t <= 1.0) and (t_aux >= 0 and t_aux <= 1.0)):  
                        #plot_point(int_l,axes,'c')
                        if start_area_calc == False:
                            start_area_calc = True
                            area_inc = abs(np.cross(int_l-vcct_el_xy,area_p[i+1]-vcct_el_xy))/2.0
                            c_faces_lines.pop(j)
                            break
                        if start_area_calc == True:
                            end_area_calc = True
                            area_inc = abs(np.cross(area_p[i]-vcct_el_xy,int_l-vcct_el_xy))/2.0
                        break
            if start_area_calc == True:
                area += area_inc
            if end_area_calc == True:
                break
        area_i += [area]  
    
    
    #print(all_nodes_dict[next_node]['d'])
    #area = (1.0-all_nodes_dict[next_node]['d'])*area + 1.0E-6
    
    return min(area_i)

def calc_area_carvalho2(vcct_el,ind_ref_adj_el,vcct_el_dict,
                            all_nodes_dict,all_els_dict):
    
    
    # if vcct_el == 37396:
    #     print('here')
    #     fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9, 9))
    #     axes.set_xlim([49.0,50.5])
    #     axes.set_ylim([9.3,10.5])
    #     plot_elements_dict(axes,all_els_dict,all_nodes_dict,'g')
        
    el_type = vcct_el_type(vcct_el,vcct_el_dict)
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
    
   
    
    
    if el_type == 'e':
        area_p = np.zeros((5,2))
        vcct_el_nodes = vcct_el_dict[vcct_el]
        inc_dir = 1
        if vcct_el_nodes[index_roll(ind_ref_adj_el+1,8)] == 0:
            inc_dir = -1

        area_p[0] = all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p[1] = all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]*0.5+all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+inc_dir,8)]]['coords'][:2]*0.5
        area_p[2] = vcct_el_xy*0.5 + all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+inc_dir*2,5)]]['coords'][:2]*0.5
        area_p[3] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+3*inc_dir,8)]]['coords'][:2]*0.5 +  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+4*inc_dir,8)]]['coords'][:2]*0.5
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4*inc_dir,8)]
        area_p[4] = all_nodes_dict[next_node]['coords'][:2]
        
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        
        c_face1 = line2D()
        c_face1.n = cf_v[0]
        c_face1.point = vcct_el_xy
        c_faces_lines = [c_face1]
    else:
            
        area_p = np.zeros((8,2))
        vcct_el_nodes = vcct_el_dict[vcct_el]

        #points defining the maximum area associated with the node
        area_p[0] =  all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]
        area_p[1] = all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]*0.5+all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+1,8)]]['coords'][:2]*0.5
        area_p[2] = vcct_el_xy*0.5 + all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+2,8)]]['coords'][:2]*0.5
        area_p[3] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+3,8)]]['coords'][:2]*0.5 +  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]]['coords'][:2]*0.5
        next_node = vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]
        area_p[4] = all_nodes_dict[next_node]['coords'][:2]
        area_p[5] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+4,8)]]['coords'][:2]*0.5 +  all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+5,8)]]['coords'][:2]*0.5
        area_p[6] = vcct_el_xy*0.5 + all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+6,8)]]['coords'][:2]*0.5
        area_p[7] = all_nodes_dict[vcct_el_nodes[index_roll(ind_ref_adj_el+7,8)]]['coords'][:2]*0.5 + all_nodes_dict[vcct_el_nodes[ind_ref_adj_el]]['coords'][:2]*0.5
        
        
       
        cf_v = all_nodes_dict[vcct_el]['cf_v']
        c_face1 = line2D()
        c_face1.n = cf_v[0,:]
        c_face1.point = vcct_el_xy
        
        c_face2 = line2D()
        c_face2.n = cf_v[1,:]
        c_face2.point = vcct_el_xy
        c_faces_lines = [c_face1,c_face2]
        
        # if vcct_el == 37396:
        #     for p in area_p:
        #         plot_point(p,axes,'r')
        #     plot_crack_front_el(cf_v, vcct_el_xy, axes)
            
    start_area_calc = False
    end_area_calc = False
    area = 0
    for i,pi in enumerate(area_p[:-1]):
        line_pi_pip1 = line2D()
        line_pi_pip1.n =  area_p[i+1] - pi
        line_pi_pip1.point = pi
        area_inc = abs(np.cross(pi-vcct_el_xy,area_p[i+1]-vcct_el_xy))/2.0
        
        for j,c_line in enumerate(c_faces_lines):
            int_l,t,col = line_pi_pip1.intersect_w_line(c_line)
            
            if len(int_l) > 0:
                int_l,t_aux,col = c_line.intersect_w_line(line_pi_pip1)
                if ((t >= 0 and t <= 1.0) and (t_aux >= 0 and t_aux <= 1.0)):  
                    if start_area_calc == False:
                        start_area_calc = True
                        area_inc = abs(np.cross(int_l-vcct_el_xy,area_p[i+1]-vcct_el_xy))/2.0
                        c_faces_lines.pop(j)
                        break
                    if start_area_calc == True:
                        end_area_calc = True
                        area_inc = abs(np.cross(area_p[i]-vcct_el_xy,int_l-vcct_el_xy))/2.0
                        break
        if start_area_calc == True:
            area += area_inc
        if end_area_calc == True:
            break
        
    area = (1.0-all_nodes_dict[next_node]['d'])*area + 1.0E-6

    return area



def interpol_el_nodal_val_wnat_coords(quantity,nat,el,val_dim,shape_f_type,
                                      all_els_dict,all_nodes_dict): 
    if 'quad_linear':
        val_nodes = np.zeros((4,val_dim))
        for j,node in enumerate(all_els_dict[el]['nodes']):
            if quantity in all_nodes_dict[node]:
                val_nodes[j,:] = all_nodes_dict[node][quantity]
        n_shape = shape_f(nat,'quad_linear')
        val = np.array([np.dot(n_shape,val_nodes[:,i]) for i in range(val_dim)])
    else:
        print("shape function not defined")
    
    
    return val


def calc_area_sample_points_asc(vcct_el,ind_ref,sample_points_xy,a_node_prist_along_normal,vcct_el_dict,
                            all_nodes_dict,all_els_dict,option='carvalho'):
    
    no_sample_points = len(sample_points_xy[:,0])
    area_sample_points = np.zeros(no_sample_points)
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
    cf_v = all_nodes_dict[vcct_el]['cf_v']
    
    if option == 'boeing':
        nodes_vcct_el = vcct_el_dict[vcct_el]
        area = 0
        for i,node in enumerate(nodes_vcct_el):
            if nodes_vcct_el[index_roll(i+1,8)] == 0:
                break
            else:
                coords_next = all_nodes_dict[nodes_vcct_el[index_roll(i+1,8)]]['coords'][:2]
                coords_curr = all_nodes_dict[nodes_vcct_el[i]]['coords'][:2]
                area += np.cross((coords_curr-vcct_el_xy),(coords_next - vcct_el_xy))/4.0
        area_sample_points[:] = area
    elif option == 'carvalho':
        for i,sp in enumerate(sample_points_xy):
            c_dir = vcct_el_xy[:2] - sp
            area_sample_points[i] = 1.0/2.0*(abs(np.cross(cf_v[0,:2],c_dir) + np.cross(c_dir,cf_v[1,:2])))
            if area_sample_points[i] == 0.0:
                area_sample_points[i] = 1e-6
    elif option == 'carvalho2':
        a_nodes_indexes = [0,2,4,6]
        ind_ref_adj_el = a_nodes_indexes[ind_ref]
        area_sample_points[0] = calc_area_carvalho7(vcct_el,ind_ref_adj_el,vcct_el_dict,
                                    all_nodes_dict,all_els_dict)
    return area_sample_points

def calc_G_sample_points(vcct_el,delta_sample_points,force_sample_points,area,
                    all_nodes_dict):
    #if option == 'nonlocal':
    no_sample_points = len(delta_sample_points[:,0]) 
    G_sample_points = np.zeros((no_sample_points,3))
    
    f_vcct_el = all_nodes_dict[vcct_el]['force']
    delta_vcct_el = all_nodes_dict[vcct_el]['delta']
    for i,delta_sp in enumerate(delta_sample_points):
        G_sample_points[i,:] = (f_vcct_el*delta_sp + force_sample_points[i,:]*delta_vcct_el)/(2.0*area[i]) 
            
    return G_sample_points


# def calc_G_sample_points(vcct_el,sample_points_xy,delta_sample_points,force_sample_points,area_sample_points,
#                     all_nodes_dict,option = 'nonlocal'):
    
#     no_sample_points = len(sample_points_xy[:,0]) 
#     G_sample_points = np.zeros((no_sample_points,3))
#     if option == 'local' or option == 'local_upd_area':
#         if 'fmax' in all_nodes_dict[vcct_el]:
#             fmax = all_nodes_dict[vcct_el]['fmax'] 
#         else:
#             fmax = all_nodes_dict[vcct_el]['force']
            
#         if 'delta_max' in all_nodes_dict[vcct_el]:
#             delta_sample_points = np.array([all_nodes_dict[vcct_el]['delta_max']])
            
#         if option == 'local':
#             if 'area_max' in all_nodes_dict[vcct_el]:
#                 area_sample_points = np.array([all_nodes_dict[vcct_el]['area_max']]*no_sample_points)
        
#         for i,delta_sp in enumerate(delta_sample_points):
#             G_sample_points[i,:] = fmax*delta_sp/(2.0*area_sample_points[i])
#     elif option == 'nonlocal':
#         f_vcct_el = all_nodes_dict[vcct_el]['force']
#         delta_vcct_el = all_nodes_dict[vcct_el]['delta']
#         for i,delta_sp in enumerate(delta_sample_points):
#             G_sample_points[i,:] = (f_vcct_el*delta_sp + force_sample_points[i,:]*delta_vcct_el)/(2.0*area_sample_points[i]) 
            
#     return G_sample_points


def write_pnewdt(filename,pnewdt):
    with open(filename,'wb') as f_handle:
        np.savetxt(f_handle,pnewdt,fmt='%f')


def adj_el_to_vcct_node_from2nodes(vcct_el,nodes_arr,node_2vcct_adj_els,all_els_dict):

    #given an array with two nodes, finds the adjacent element to 
    #a vcct element that contains t
    
    adj_els = node_2vcct_adj_els[vcct_el]
    for el in adj_els:
        nodes_of_adj_el = all_els_dict[el]
        common_nodes = np.intersect1d(nodes_of_adj_el['nodes'],nodes_arr)
        if len(common_nodes) == 2:
            return el
    return None
    


# def failed_area_inc(vcct_el,vcct_el_coords,cw_v,nodes_arr,cf_v_i,all_nodes_dict,
#                     node_2vcct_adj_els_dict,all_els_dict):
    
    
#     f_area_inc = 0.0
#     adj_el = adj_el_to_vcct_node_from2nodes(vcct_el,nodes_arr,node_2vcct_adj_els_dict,all_els_dict)
#     if adj_el != None:
#         adj_el_nodes = all_els_dict[adj_el]['nodes']
        
#         vcct_node_ind = np.nonzero(adj_el_nodes == vcct_el)[0][0]
#         a_node_1_ind = index_roll(vcct_node_ind+1,4)
#         a_node_2_ind = index_roll(vcct_node_ind-1,4)
#         a_node_1_coords = all_nodes_dict[adj_el_nodes[a_node_1_ind]]['coords']
#         a_node_2_coords = all_nodes_dict[adj_el_nodes[a_node_2_ind]]['coords']
        
#         a_node_1_nat_coords = inverse_iso_p(adj_el,a_node_1_coords,all_els_dict,all_nodes_dict)
#         a_node_2_nat_coords = inverse_iso_p(adj_el,a_node_2_coords,all_els_dict,all_nodes_dict)
        
#         line_diagonal_nat = line2D()
#         line_diagonal_nat.line_two_points(a_node_1_nat_coords,a_node_2_nat_coords)
#         #based on a_node1/a_node_2 define the element diagonal connecting the two antena nodes
           
        
#         #find crack face line in natural coordinates
#         crack_pos = vcct_el_coords + cf_v_i
#         vcct_el_nat_coords = inverse_iso_p(adj_el,vcct_el_coords,all_els_dict,all_nodes_dict)
#         crack_pos_nat_coords = inverse_iso_p(adj_el,crack_pos,all_els_dict,all_nodes_dict)

#         line_cfront_nat = line2D()
#         line_cfront_nat.line_two_points(vcct_el_nat_coords,crack_pos_nat_coords)
        
#         intersection_nat_coords,t0, c_linear = line_diagonal_nat.intersect_w_line(line_cfront_nat)
#         cf_intersection = len(intersection_nat_coords) == 1
        
#         if len(intersection_nat_coords):
#             intersection_coords = interpol_el_nodal_val_wnat_coords('coords',intersection_nat_coords,adj_el,3,'quad_linear',all_els_dict,all_nodes_dict)
        
#         v1 = cw_v
#         v2 = intersection_coords - vcct_el_coords
#         f_area_inc = abs(np.linalg.norm(np.cross(v1,v2)))/2.0    
        

#     return f_area_inc,cf_intersection



# def failed_area_inc(vcct_el,vcct_el_coords,cw_v,nodes_arr,cf_v_i,all_nodes_dict,
#                     node_2vcct_adj_els_dict,all_els_dict):
    
    
#     f_area_inc = 0.0
#     adj_el = adj_el_to_vcct_node_from2nodes(vcct_el,nodes_arr,node_2vcct_adj_els_dict,all_els_dict)
#     if adj_el != None:
#         adj_el_nodes = all_els_dict[adj_el]['nodes']
        
#         vcct_node_ind = np.nonzero(adj_el_nodes == vcct_el)[0][0]
#         a_node_1_ind = index_roll(vcct_node_ind+1,4)
#         a_node_2_ind = index_roll(vcct_node_ind-1,4)
#         a_node_1_coords = all_nodes_dict[adj_el_nodes[a_node_1_ind]]['coords']
#         a_node_2_coords = all_nodes_dict[adj_el_nodes[a_node_2_ind]]['coords']
        
#         a_node_1_nat_coords = inverse_iso_p(adj_el,a_node_1_coords,all_els_dict,all_nodes_dict)
#         a_node_2_nat_coords = inverse_iso_p(adj_el,a_node_2_coords,all_els_dict,all_nodes_dict)
        
#         line_diagonal_nat = line2D()
#         line_diagonal_nat.line_two_points(a_node_1_nat_coords,a_node_2_nat_coords)
#         #based on a_node1/a_node_2 define the element diagonal connecting the two antena nodes
           
        
#         #find crack face line in natural coordinates
#         crack_pos = vcct_el_coords + cf_v_i
#         vcct_el_nat_coords = inverse_iso_p(adj_el,vcct_el_coords,all_els_dict,all_nodes_dict)
#         crack_pos_nat_coords = inverse_iso_p(adj_el,crack_pos,all_els_dict,all_nodes_dict)

#         line_cfront_nat = line2D()
#         line_cfront_nat.line_two_points(vcct_el_nat_coords,crack_pos_nat_coords)
        
#         intersection_nat_coords,t0, c_linear = line_diagonal_nat.intersect_w_line(line_cfront_nat)
        
#         intersection_coords = interpol_el_nodal_val_wnat_coords('coords',intersection_nat_coords,adj_el,3,'quad_linear',all_els_dict,all_nodes_dict)
        
#         v1 = cw_v
#         v2 = intersection_coords - vcct_el_coords
#         f_area_inc = abs(np.linalg.norm(np.cross(v1,v2)))/2.0    
        

#     return f_area_inc,v2


def nodes_vcct(all_els_dict):

    nodes = np.array([el_nodes['nodes'] for el_no,el_nodes in all_els_dict.items()])
    nodes = nodes.flatten()
    
    return np.unique(nodes.flatten())


def assign_node_prop_from_el_set(el_set,all_els_dict,all_nodes_dict,prop_val,prop_key):

    for el in el_set:
        nodes_el = all_els_dict[el]['nodes']
        for node in nodes_el:
            for p_val,p_key in zip(prop_val,prop_key):
                if p_key in all_nodes_dict[node].keys():
                    if all_nodes_dict[node][p_key] < p_val:
                        all_nodes_dict[node][p_key] = p_val
                else:
                    all_nodes_dict[node][p_key] = p_val
                    
    return all_nodes_dict


def assign_fract_params(mat_list,set_list,all_els_dict,all_nodes_dict):
    """based on the element material specifications made in the input file, assign the properties to the nodes
    for vcct calculations"""
                        
    no_mats = len(mat_list)
    for i in range(no_mats):
        
        d =  mat_list[i][ini.ind_ini_dstat]
        
        theta_T = mat_list[i][ini.ind_ini_theta_T]
        theta_B = mat_list[i][ini.ind_ini_theta_B]
        GIC_B  = mat_list[i][ini.ind_ini_GIC_B]
        GIIC_B = mat_list[i][ini.ind_ini_GIIC_B]
        BK_B   = mat_list[i][ini.ind_ini_BK_B]
        
        GIC_T  = mat_list[i][ini.ind_ini_GIC_T]
        GIIC_T = mat_list[i][ini.ind_ini_GIIC_T]
        BK_T   = mat_list[i][ini.ind_ini_BK_T]
        
        Z_T = mat_list[i][ini.ind_ini_Z_T]
        Z_S = mat_list[i][ini.ind_ini_Z_S]
        
        dadn_B = mat_list[i][ini.ind_ini_dadn_B]
        dadn_T = mat_list[i][ini.ind_ini_dadn_T]
        
        r_curve_B = mat_list[i][ini.ind_ini_r_curve_B]
        r_curve_T = mat_list[i][ini.ind_ini_r_curve_T]
        
        props_vals = np.array([d,theta_T,theta_B,GIC_B,GIIC_B,BK_B,GIC_T,GIIC_T,BK_T,Z_T,Z_S,dadn_B,dadn_T,r_curve_B,r_curve_T])
        props_keys = ['d','theta_T','theta_B','GIC_B','GIIC_B','BK_B','GIC_T','GIIC_T','BK_T','Z_T','Z_S','dadn_B','dadn_T','r_curve_B','r_curve_T']
        all_nodes_dict = assign_node_prop_from_el_set(set_list[i],all_els_dict,all_nodes_dict,props_vals,props_keys)
       

    return all_nodes_dict


def save_vcct_data(current_dir,list_of_dicts,list_of_names):

    for vcct_dict,name in zip(list_of_dicts,list_of_names):  
        with open(current_dir+name+'.pickle', 'wb') as handle:
            pickle.dump(vcct_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    return


def save_dadn_r_curve_data(current_dir,dadn_mm_table,rcurve_mm_table):
    with open(current_dir+'dadn_mm_table.pickle', 'wb') as handle:
        pickle.dump(dadn_mm_table, handle,protocol=pickle.HIGHEST_PROTOCOL)
    with open(current_dir+'rcurve_mm_table.pickle', 'wb') as handle:
        pickle.dump(rcurve_mm_table, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return

def load_vcct_data(current_dir,list_of_names):
    list_of_dicts = []
    for name in list_of_names:
        with open(current_dir+name+'.pickle', 'rb') as handle:
            if sys.version_info.major > 2:
                list_of_dicts += [pickle.load(handle,encoding='latin1')]
            else:
                list_of_dicts += [pickle.load(handle)]
       
    return list_of_dicts


def load_dadn_r_curve_data(current_dir):
    
    with open(current_dir+'dadn_mm_table.pickle', 'rb') as handle:
        if sys.version_info.major > 2:
            dadn_mm_table = pickle.load(handle,encoding='latin1')
        else:
            dadn_mm_table = pickle.load(handle)
        
    with open(current_dir+'rcurve_mm_table.pickle', 'rb') as handle:
        if sys.version_info.major > 2:
            rcurve_mm_table = pickle.load(handle,encoding='latin1')
        else:
            rcurve_mm_table = pickle.load(handle)

    return dadn_mm_table,rcurve_mm_table


def all_nodes_failed(all_nodes_dict):
    
    nodes_d = [node for node, node_info in all_nodes_dict.items() if node_info['d'] == 1.0 ]

    return nodes_d


def assign_nodal_info(el,ind_node,tip_els,node,all_nodes_dict):
    
    #tip_els info:
    #[el_no,m_x_node_1,m_x_node_2,m_x_node_3,m_x_node_4,
    #       m_y_node_1,m_y_node_2,m_y_node_3,m_y_node_4,
    #       m_z_node_1,m_z_node_2,m_z_node_3,m_z_node_4,
    #       pnlty_f_node_1,penalty_f_node_2,penalty_f_node_3,penalty_f_node_4,
    #       d_node_1,d_node_2,d_node_3,d_node_4]
    no_nodes = 4
    no_dofs = 3

    if el not in tip_els.keys():
        tip_els[el] = np.array([0.0]*21)
        tip_els[el][0] = el
        
    #write slopes: 'm' (vector)  
    tip_arr_offset = 1
    if 'm' in all_nodes_dict[node].keys():
        for i,m_i in enumerate(all_nodes_dict[node]['m']):
            tip_els[el][tip_arr_offset+ind_node+i*no_nodes] = m_i
    tip_arr_offset = tip_arr_offset + no_nodes*no_dofs
    
    #writes 'pnlty_f'
    if 'pnlty_f' not in all_nodes_dict[node].keys():
        all_nodes_dict[node]['pnlty_f'] = 2 
    tip_els[el][tip_arr_offset+ind_node]  = all_nodes_dict[node]['pnlty_f']    
    tip_arr_offset = tip_arr_offset + no_nodes    
    
    #writes 'd'
    tip_els[el][tip_arr_offset+ind_node] = all_nodes_dict[node]['d']
   
    return tip_els


def reset_force_delta(vcct_els_cfront,onset_node_list,all_nodes_dict,node_2vcct_adj_els_dict,all_els_dict):

    #may be optimized to avoid redundancy (nodes will be reset several times)
    #elements are used to ensure nodes reset are consistent with tip elements provided to abaqus
    for vcct_el in vcct_els_cfront:
        adj_els = node_2vcct_adj_els_dict[vcct_el]
        for el in adj_els:
            nodes = all_els_dict[el]['nodes']
            for node in nodes:
                all_nodes_dict[node]['force'] = np.array([0,0,0])
                all_nodes_dict[node]['delta'] = np.array([0,0,0])

    for node in onset_node_list:
        all_nodes_dict[node]['force_onset'] = np.array([0,0,0])
        
    return all_nodes_dict

def map_vcctnodes2els(nodes,node_2vcct_adj_els,all_nodes_dict,all_els_dict):
    
    tip_els = {}
    for node in nodes:
        #finds adjacent elements
        adj_els = node_2vcct_adj_els[node]
        #for each element updates the nodal information
        for el in adj_els:
            #finds the index to the node to be updated in the element being considered
            nodes = all_els_dict[el]['nodes']
            ind_node = np.nonzero(nodes == node)[0]
            tip_els = assign_nodal_info(el,ind_node,tip_els,node,all_nodes_dict)
           
    tip_els_arr = np.array([info for t_el, info in tip_els.items()])
            
    return tip_els_arr
    

def write_to_file_els(filename,els):
    if len(els) != 0:
        lenels = np.array([len(els[:,0])])
    else:
        lenels = np.array([0])
    
    with open(filename,'wb') as f_handle:
        np.savetxt(f_handle,lenels,fmt='%d')
        if lenels != 0:
            l_format = ['%d'] #element number
            no_variables = len(els[0,:])
            if no_variables > 1:
                l_format.extend(['%1.17E']*(no_variables-1))
            np.savetxt(f_handle,els,fmt=l_format)
        

def read_nodal_damage_map_info(current_dir,all_nodes_dict):
    
    filename = current_dir+'node_map_info.ndmg'
    if os.path.isfile(filename):
        print('reading...')
        #assumes nodal map is given as an abaqus comma separated set (no set same)
        with open(filename) as f:
            for line in f:
                line = line.strip() #added to remove \r\n
                s = np.fromstring(line, dtype=int, count=-1, sep=',')
                for node in s:
                    all_nodes_dict[node]['d'] = 1.0
    else:
        print('no damage map info found')
        
    return all_nodes_dict


def read_lead_nodes(current_dir,all_nodes_dict,analysis_dict):

    analysis_dict['lead_nodes_bool'] = False
    filename = current_dir+'lead_nodes.nlead'
    if os.path.isfile(filename):
        with open(filename) as f:
            for line in f:
                line = line.strip() #added to remove \r\n
                s = np.fromstring(line, dtype=int, count=-1, sep=',')
                for node in s:
                    all_nodes_dict[node]['lead_node_bool'] = True
                    analysis_dict['lead_nodes_bool'] = True    
        
    return all_nodes_dict,analysis_dict

def minmax_interval(interval,coord):
    
    if interval[0] > coord:
        interval[0] = coord
    elif interval[1] < coord:
        interval[1] = coord

    return interval


def coords_corner_of_bucket(bucket):
    
    no_of_coners = 8
    dofs = 3
    coords_corners = np.zeros((no_of_coners,dofs))
    coords_corners[0,:] = np.array([bucket[0],bucket[1],bucket[2]])
    coords_corners[1,:] = np.array([bucket[0]+bucket[3],bucket[1],bucket[2]])
    coords_corners[2,:] = np.array([bucket[0]+bucket[3],bucket[1]+bucket[4],bucket[2]])
    coords_corners[3,:] = np.array([bucket[0],bucket[1]+bucket[4],bucket[2]])
    coords_corners[4:,:] = coords_corners[:4,:] + np.array([[0,0,bucket[5]]]*4)
    
    return coords_corners
    
def create_bounding_cube(center,edge_length):
    
    half_edge_length = edge_length/2.0
    cube_coords = np.zeros((8,3))
    cube_coords[0,:] = center + np.array([-half_edge_length,-half_edge_length,-half_edge_length])
    cube_coords[1,:] = center + np.array([half_edge_length,-half_edge_length,-half_edge_length])
    cube_coords[2,:] = center + np.array([half_edge_length,half_edge_length,-half_edge_length])
    cube_coords[3,:] = center + np.array([-half_edge_length,half_edge_length,-half_edge_length])
    cube_coords[4:,:] = cube_coords[:4,:] + np.array([[0,0,edge_length]]*4)
    
    return cube_coords

def point_in_cube(point,cube):
    
    cond_x = point[0] >=  cube[0,0] and  point[0] <= cube[1,0]
    cond_y = point[1] >=  cube[0,1] and  point[1] <= cube[2,1]
    cond_z = point[2] >=  cube[0,2] and  point[2] <= cube[4,2]
   
    return cond_x*cond_y*cond_z

    

def reset_no_fail_zones(all_nodes_dict):
    
    for node,info in all_nodes_dict.items():
        if 'nodes_in_cube' in all_nodes_dict[node].keys():
            del all_nodes_dict[node]['nodes_in_cube']
            
    return all_nodes_dict
    
    

def setup_no_fail_zones(all_nodes_dict):
    
    if hasattr(inp_vars,'no_fail_zone_char_length'):
        nfz_length = inp_vars.no_fail_zone_char_length

        #t1 = time.time()
        
        ##first key
        node_no = next(iter(all_nodes_dict))
        
        #selects min_max_coords
        x_interval = [all_nodes_dict[node_no]['coords'][0],
                      all_nodes_dict[node_no]['coords'][0]]
        y_interval =[all_nodes_dict[node_no]['coords'][1],
                     all_nodes_dict[node_no]['coords'][1]]
        z_interval =[all_nodes_dict[node_no]['coords'][2],
                     all_nodes_dict[node_no]['coords'][2]]
        
        if hasattr(inp_vars,'no_of_buckets'):
            no_of_buckets = inp_vars.no_of_buckets
        else:
            no_of_buckets = 4
            
            
        
        for node, node_info in all_nodes_dict.items():
            node_coords = node_info['coords']
            x_interval = minmax_interval(x_interval,node_coords[0])
            y_interval = minmax_interval(y_interval,node_coords[1])
            z_interval = minmax_interval(z_interval,node_coords[2])
        
        #creates buckets
        x_arr_bucket = np.linspace(x_interval[0],x_interval[1], num=no_of_buckets+1)
        y_arr_bucket = np.linspace(y_interval[0],y_interval[1], num=no_of_buckets+1)
        z_arr_bucket = np.linspace(z_interval[0],z_interval[1], num=no_of_buckets+1)
        
        bucket_counter = 0
        all_buckets_dict = {}
        for l,z in enumerate(z_arr_bucket[:-1]):
            for m,y in enumerate(y_arr_bucket[:-1]):
                for j,x in enumerate(x_arr_bucket[:-1]):
                    bucket_xyz_delta_xyz = [x,y,z,
                                            x_arr_bucket[j+1]-x,
                                            y_arr_bucket[m+1]-y,
                                            z_arr_bucket[l+1]-z]
                    all_buckets_dict[bucket_counter] = {}
                    all_buckets_dict[bucket_counter]['cube'] = coords_corner_of_bucket(bucket_xyz_delta_xyz)
                    all_buckets_dict[bucket_counter]['nodes'] = []
                    bucket_counter += 1
                    all_buckets = list(all_buckets_dict.keys())
                    
        for node,node_info in all_nodes_dict.items():
            node_coords = node_info['coords']
            for bucket_no,bucket_info in all_buckets_dict.items():
                if point_in_cube(node_coords,bucket_info['cube']):
                    all_buckets_dict[bucket_no]['nodes'] += [node] 
                    all_nodes_dict[node]['bucket'] = bucket_no
                    break
                
        for node_i,node_info in all_nodes_dict.items():
            
            all_nodes_dict[node_i]['nodes_in_cube'] = []
            #creates a cube associated with the node i:
            cube_node_i = create_bounding_cube(all_nodes_dict[node_i]['coords'],nfz_length)
            
            #obtains bucket of node i
            bucket_of_node_i = all_nodes_dict[node_i]['bucket']
            bucket_info = all_buckets_dict[bucket_of_node_i]
            bucket_search_list = [bucket_of_node_i]
            
            #check if volume is fully whithin the node's bucket
            #check if volume is whithin the node's bucket
            vol_whithin = True
            for corner in cube_node_i:
                if point_in_cube(corner,bucket_info['cube']) == False:
                    #vol_whithin == False
                    vol_whithin = False
                    break
                #vol_whithin = False 
            
            #check how many buckets have at least a corner whithin the node cube
            if vol_whithin == False:
                all_other_buckets = [bucket for bucket in all_buckets if bucket != bucket_of_node_i]
                for bucket in all_other_buckets:
                    for corner_cube in all_buckets_dict[bucket]['cube']:
                        if point_in_cube(corner_cube,cube_node_i):
                            bucket_search_list += [bucket]
                            break
        
            # if vol_whithin == False:
            #     print('here')
            #     all_other_buckets = [bucket for bucket in all_buckets if bucket != bucket_of_node_i]
            #     print('all_buckets',all_buckets)
            #     for bucket in all_other_buckets:
            #         for corner in cube_node_i:
            #             if point_in_cube(corner,all_buckets_dict[bucket]['cube']):
            #                 bucket_search_list += [bucket]
            #                 print('here-2')
            #                 break
                        
            #retrieves nodes in bucket_search_list
            for bucket in bucket_search_list:
                nodes_in_bucket= all_buckets_dict[bucket]['nodes']
                for node_ib in nodes_in_bucket:
                    node_ib_coords = all_nodes_dict[node_ib]['coords']
                    if point_in_cube(node_ib_coords,cube_node_i):
                        all_nodes_dict[node_i]['nodes_in_cube'] += [node_ib]
    else:
        print('no fail zone not defined')
                        
    return all_nodes_dict


def initialize_analysis_options(analysis_dict):
        
    PR_default = True
    max_normalized_inc_default = 0.2
    cut_back_factor_default = 0.1
    static_scale_f_default = 10.0
    stress_crit_tol_default = 0.05
    unst_gth_tol_default = 0.05
    
    analysis_dict['PR'] = PR_default
    analysis_dict['max_normalized_inc'] = max_normalized_inc_default
    analysis_dict['cut_back_factor'] = cut_back_factor_default
    analysis_dict['static_scale_f'] = static_scale_f_default
    analysis_dict['stress_crit_tol'] = stress_crit_tol_default
    analysis_dict['unst_gth_tol'] = unst_gth_tol_default
    
    
    if hasattr(inp_vars,'PR'):
        analysis_dict['PR'] = inp_vars.PR
    
    if hasattr(inp_vars,'max_normalized_inc'): 
        analysis_dict['max_normalized_inc'] = inp_vars.max_normalized_inc
        
    if hasattr(inp_vars,'cut_back_factor'):
        analysis_dict['cut_back_factor']  = inp_vars.cut_back_factor
        
    if hasattr(inp_vars,'static_scale_f'):
        analysis_dict['static_scale_f'] = inp_vars.static_scale_f

        
    if hasattr(inp_vars,'stress_crit_tol'): 
        analysis_dict['stress_crit_tol'] = inp_vars.stress_crit_tol
        
    if hasattr(inp_vars,'unst_gth_tol'): 
        analysis_dict['unst_gth_tol']  = inp_vars.unst_gth_tol
    
    #print outs:
    not_default = False
    out_string = 'PR (progressive release)'
    if analysis_dict['PR'] != PR_default:
        out_string = out_string + '**'
        not_default = True
    print(out_string,analysis_dict['PR'])
        
    out_string = 'normalized inc'
    if analysis_dict['max_normalized_inc'] != max_normalized_inc_default:
       out_string = out_string + '**'
       not_default = True
    print(out_string,analysis_dict['max_normalized_inc'])

    out_string = 'unst gth tol'
    if analysis_dict['unst_gth_tol'] != unst_gth_tol_default:
       out_string = out_string + '**'
       not_default = True
    print(out_string,analysis_dict['unst_gth_tol'])
    
    out_string = 'cut back factor'
    if analysis_dict['cut_back_factor']  != cut_back_factor_default:
        out_string = out_string + '**'
        not_default = True
    print(out_string,analysis_dict['cut_back_factor'])
    
    out_string = 'static scale f'
    if  analysis_dict['static_scale_f'] != static_scale_f_default:
        out_string = out_string + '**'
        not_default = True
    print(out_string,analysis_dict['static_scale_f'] )
        
    out_string = 'stress crit tol'
    if analysis_dict['stress_crit_tol'] != stress_crit_tol_default:
        out_string = out_string + '**'
        not_default = True
    print(out_string,analysis_dict['stress_crit_tol'])

    if 'lead_nodes_bool' in analysis_dict.keys():
        if analysis_dict['lead_nodes_bool'] == True:
            out_string = 'propagation based on lead nodes**'
            not_default = True
         
    if not_default:
        print('** - not default')
    return analysis_dict


def ctod_reporting_option(analysis_dict):
    
    ctod_bool_default = False
    nsearches_ctod_default = 3
    ref_vect_ctod_default = np.array([-1.0,0.0,0.0])
    
    analysis_dict['ctod_bool'] = ctod_bool_default
    analysis_dict['nsearches_ctod'] = nsearches_ctod_default
    analysis_dict['ref_vect_ctod'] = ref_vect_ctod_default
    
    if hasattr(inp_vars, 'ctod_bool'):
        analysis_dict['ctod_bool'] = inp_vars.ctod_bool
        
        if inp_vars.ctod_bool == True:
            if hasattr(inp_vars, 'nsearches_ctod'):
                analysis_dict['nsearches_ctod'] = inp_vars.nsearches_ctod
                
            if hasattr(inp_vars, 'ref_vect_ctod'):
                analysis_dict['ref_vect_ctod'] = inp_vars.ref_vect_ctod
                
            print('ctod reporting option activated')
            print('nsearches_ctod',analysis_dict['nsearches_ctod'])
            print('ref_vect_ctod',analysis_dict['ref_vect_ctod'])
            
            # nodes_coords_ctod_arr = ctod_calc(el_number,ref_vect,nsearches_ctod,vcct_el_dict,all_nodes_dict)
            # filename = current_dir + 'ctod'+str(int(step))+'_'+str(int(inc))+'.txt'
            # with open(filename,'ab') as f_handle:
            #     np.savetxt(f_handle,nodes_coords_ctod_arr,fmt='%f') 
               
    return analysis_dict

def VCCT_init():
    
    print('Initilization...')
    if test_init == True:
        print('Load elements/nodes test:')
        all_els_dict,all_nodes_dict = tn.test_el_node_data()
        
    else:
        print('Loading VCCT elements...')
        ti = time.time()
        all_els_dict = ini.read_elements(current_dir + input_file + '.inp')
        tf =  time.time()
        print('done, time elapsed:',tf-ti)
        
    
        print('Loading VCCT nodes...')
        # #obtains a list with all nodes belonging to the user elements
        ns_vcct = nodes_vcct(all_els_dict)
        
        # #reads nodes in the '_VCCT...' part in the *.inp file
        # #that belowng to the VCCT elements
        all_nodes_dict = ini.read_nodes(current_dir +input_file+'.inp',ns_vcct)
        ti = tf
        tf =  time.time()
        print('done, time elapsed:',tf-ti)
    
    print('Building VCCT elements...')
    #for each VCCT node determines the adjacent elements
    node_2vcct_adj_els_dict = build_node_2vcct_adj_els_dict_fast(all_els_dict)

    #builds a 'vcct element' corresponding to each node
    vcct_el_dict = build_vcct_el_dict(node_2vcct_adj_els_dict,all_nodes_dict,all_els_dict)
    ti = tf
    tf =  time.time()
    print('done, time elapsed:',tf-ti)

    #global variables helpful to keep track of
    analysis_dict = {}
     
    print('Defining no fail zones...')
    all_nodes_dict = setup_no_fail_zones(all_nodes_dict)
    ti = tf
    tf =  time.time()
    print('done, time elapsed:',tf-ti)
    
    
    if test_init == False:
        print('Reading material definition...')
        mat_list,set_list,dadn_mm_table,rcurve_mm_table = ini.read_mat_list_elset(current_dir +input_file + '.inp')
        
        all_nodes_dict = assign_fract_params(mat_list,set_list,all_els_dict,all_nodes_dict)
        ti = tf
        tf =  time.time()
        print('done, time elapsed:',tf-ti)
    
    print('reading nodal damage mapping info...')
    all_nodes_dict = read_nodal_damage_map_info(current_dir,all_nodes_dict)
    ti = tf
    tf =  time.time()
    print('done, time elapsed:',tf-ti)

    print('reading leading nodes info...')
    all_nodes_dict,analysis_dict = read_lead_nodes(current_dir,all_nodes_dict,analysis_dict)
    ti = tf
    tf =  time.time()
    print('done, time elapsed:',tf-ti)
    
    analysis_dict = ctod_reporting_option(analysis_dict)
    
    print('Preparing damage initialization files...')
    #identifies all failed nodes
    nodes_deq1 = all_nodes_failed(all_nodes_dict)
    #identifies nodes at the vcct crack front
    vcct_els_cfront = vcct_els_crackfront(vcct_el_dict,all_nodes_dict,option_node = 'dneq1',option_a_nodes = 'dgt0')
    if analysis_dict['ctod_bool'] == True:
        vcct_els_cfront = expand_cfront_for_ctod_extraction(vcct_els_cfront,analysis_dict['ref_vect_ctod'],analysis_dict['nsearches_ctod'],vcct_el_dict,all_nodes_dict)
        
    #merges the two arrays and finds a unique set of nodes:
    vcct_els_cfront = np.unique(np.append(vcct_els_cfront,nodes_deq1))

    tip_els_I = map_vcctnodes2els(vcct_els_cfront,node_2vcct_adj_els_dict,all_nodes_dict,all_els_dict)
    write_to_file_els(current_dir+'TIP_ELS.txt',tip_els_I)
    #map_to_elements
    
    print('Initializes analysis options')
    analysis_dict = initialize_analysis_options(analysis_dict)
    
    #saves all vcct data
    list_of_dicts = [analysis_dict,all_els_dict,all_nodes_dict,node_2vcct_adj_els_dict,vcct_el_dict]
    list_of_names = ['analysis_dict','all_els_dict','all_nodes_dict','node_2vcct_adj_els_dict','vcct_el_dict']
    save_vcct_data(current_dir+input_file,list_of_dicts,list_of_names)
    if test_init == False:
        save_dadn_r_curve_data(current_dir,dadn_mm_table,rcurve_mm_table)
    #initialize file that controls cut backs
    pnewdt = np.array([1.0E6])
    write_pnewdt(current_dir+'pnewdt.txt',pnewdt)
    ti = tf
    tf =  time.time()
    print('done, time elapsed:',tf-ti)
    
    print('Initialization completed')
    
    
    return


def max_normalized_inc_PR(PR,max_normalized_inc = 0.2):
    
    if (PR == False and max_normalized_inc != 1.0):
        print("max_normalized_inc updated to 1.0 since using instantaneous release")
        max_normalized_inc = 1.0
    
    return max_normalized_inc


def f_cr_static_cfront(rcurve_mm_table,vcct_els_front_dict,all_nodes_dict,analysis_dict):
    
    bool_static_f = False
    if analysis_dict['lead_nodes_bool'] == True:
        for vcct_el,vcct_info in vcct_els_front_dict.items():
            if 'lead_node_bool' in all_nodes_dict[vcct_el].keys():
                if all_nodes_dict[vcct_el]['lead_node_bool'] == True:
                    G_s_cfront = vcct_info['G_s']
                    vcct_els_front_dict[vcct_el]['f_cr_s'] = np.array([f_cr_static(G,vcct_el,rcurve_mm_table,all_nodes_dict) for G in G_s_cfront])
                    if vcct_els_front_dict[vcct_el]['ss_a'] == True:
                        f_cr_s = np.max(vcct_els_front_dict[vcct_el]['f_cr_s'])
                        if f_cr_s > 1.0:
                            bool_static_f = True
                            break
        
        for vcct_el,vcct_info in vcct_els_front_dict.items():
            G_s_cfront = vcct_info['G_s']
            if bool_static_f == True:
                vcct_els_front_dict[vcct_el]['f_cr_s'] = np.array([f_cr_s]*len(G_s_cfront))
            else:
                vcct_els_front_dict[vcct_el]['f_cr_s'] = np.array([[0.0]*len(G_s_cfront)])
            
                
    else:      
        for vcct_el,vcct_info in vcct_els_front_dict.items():
            G_s_cfront = vcct_info['G_s']
            vcct_els_front_dict[vcct_el]['f_cr_s'] = np.array([f_cr_static(G,vcct_el,rcurve_mm_table,all_nodes_dict) for G in G_s_cfront])
            
            #if  bool_static_f == False:
            if  bool_static_f == False and vcct_els_front_dict[vcct_el]['ss_a'] == True:   
                bool_static_f =  np.any(vcct_els_front_dict[vcct_el]['f_cr_s'] > 1.0)
        
    return vcct_els_front_dict,bool_static_f


def GTmax_GSmax(G):
    #prevents spurious negative GS_max
    GS_max = np.sum(abs(G[:2]))
    
    GTmax = GS_max + G[2]
        
    return GTmax,GS_max

def mode_mixity(GT_max,GS_max):
    
    if GT_max > 0:
        mmixity = GS_max/GT_max        
    else:
        mmixity = 0.0
        
    return mmixity


def r_curve_mm_interp(r_curve,a_acc,mmixity):
    
    unique_mm = np.sort(np.unique(r_curve[:,2]))
    no_of_unique_mms = len(unique_mm)
    mm_vs_f = np.zeros((no_of_unique_mms,2))
 
    #obtains a vector mm_vs_f with mode-mixity vs f_val for a given crack length
    for i,mm in enumerate(unique_mm):
        ind_mm = r_curve[:,2] == mm
        f_r = r_curve[ind_mm,0]
        a_length = r_curve[ind_mm,1]
        ind_sort = np.argsort(a_length)
        f_r_sorted = f_r[ind_sort]
        a_length_sorted =  a_length[ind_sort]
        f_val = np.interp(a_acc,a_length_sorted,f_r_sorted)
        mm_vs_f[i,0] = mm
        mm_vs_f[i,1] = f_val

    #interpolates the vector mm_vs_f given the mode-mixity mmixity 
    f_val = np.interp(mmixity,mm_vs_f[:,0],mm_vs_f[:,1])
    
    return f_val

def f_cr_static(G,vcct_el,rcurve_mm_table,all_nodes_dict,option='BK'):
        
    GIC = all_nodes_dict[vcct_el]['GIC_B']
    GIIC = all_nodes_dict[vcct_el]['GIIC_B']
    BK = all_nodes_dict[vcct_el]['BK_B']

    f_r_curve = 1.0
    if 'a_acc' in all_nodes_dict[vcct_el].keys():
        a_acc =  all_nodes_dict[vcct_el]['a_acc']
    else:
        a_acc = 0.0
    
    #prevents near zero negative values to be inlcuded in the calculation
    GT_max,GS_max = GTmax_GSmax(G)
   
    
    mmixity = mode_mixity(GT_max,GS_max)

    #selects r_curve factor given mode-mixity and accumulated crack growth
    r_curve_mm = rcurve_mm_table[all_nodes_dict[vcct_el]['r_curve_B'].astype(int)-1]     
    f_r_curve = r_curve_mm_interp(r_curve_mm,a_acc,mmixity)
    
    if option == 'BK':        
        GC = (GIC + (GIIC-GIC)*mmixity**BK)*f_r_curve
    else:
        print('fracture criterion not recognized - failure prevented, VCCT_el:',vcct_el)
        GC = 1.0E52 #prevents failure

    if GT_max/GC > 100:
        sys.exit("check case, f_cr >>>> 1.0")
        
    return GT_max/GC    
 
def dadn_static_sample(vcct_els_front_dict,stat_scale_f):
            
    
    for el,el_info in vcct_els_front_dict.items():
        vcct_els_front_dict[el]['dadn_s'] = el_info['f_cr_s']**stat_scale_f

    return vcct_els_front_dict

def select_bounds(x_mmixity,mmixity):
    
    index_gt = x_mmixity[:,0] >= mmixity
    upper_b = x_mmixity[index_gt,0] 
    index_lt = x_mmixity[:,0] <= mmixity
    lower_b = x_mmixity[index_lt,0]
    
    if len(upper_b) == 0:
        lower_b = x_mmixity[-1,0]
        upper_b = x_mmixity[-1,0]
    elif len(lower_b) == 0:
        lower_b = x_mmixity[0,0]
        upper_b = x_mmixity[0,0]
    else:
        upper_b = upper_b[0]
        lower_b = lower_b[-1]
    
    index_lower_b = np.nonzero(x_mmixity[:,0] == lower_b)[0][0]
    index_upper_b = np.nonzero(x_mmixity[:,0] == upper_b)[0][0]
    
    return np.array([index_lower_b,index_upper_b])    


def dadn_fatigue_sample(vcct_els_front_dict,dadn_mm_table,rcurve_mm_table,all_nodes_dict):

    for vcct_el,vcct_info in vcct_els_front_dict.items():
        G_s_cfront = vcct_info['G_s']
        vcct_els_front_dict[vcct_el]['dadn_s'] = np.array([dadn_mm_fatigue(G,vcct_el,dadn_mm_table,rcurve_mm_table,all_nodes_dict) for G in G_s_cfront])
    
    return vcct_els_front_dict

def read_from_dadn_table(table_no,dadn_mm_table):
    
    #nomenclature used to designate the coeficients
    c_mm_ind = 0 #index to coeficient 'c'
    beta_mm_ind = 1 #index to coeficient 'beta'
    mu_mm_ind = 2   #index to coeficient 'mu'
    gamma_mm_ind = 3 #index to coeficient 'gamma'
    rho_mm_ind = 4 #index to coeficient 'rho'
    
    mm_ind = -1 #index to mode-mixity 'mm = [0,1]
    c_mm = np.array([dadn_mm_table[table_no-1][:,mm_ind],dadn_mm_table[table_no-1][:,c_mm_ind]]).T    
    beta_mm = np.array([dadn_mm_table[table_no-1][:,mm_ind],dadn_mm_table[table_no-1][:,beta_mm_ind]]).T
    mu_mm = np.array([dadn_mm_table[table_no-1][:,mm_ind],dadn_mm_table[table_no-1][:,mu_mm_ind]]).T
    gamma_mm = np.array([dadn_mm_table[table_no-1][:,mm_ind],dadn_mm_table[table_no-1][:,gamma_mm_ind]]).T
    rho_mm = np.array([dadn_mm_table[table_no-1][:,mm_ind],dadn_mm_table[table_no-1][:,rho_mm_ind]]).T
    
    return c_mm,beta_mm,mu_mm,gamma_mm,rho_mm

def paris_law_C(c_mm,gamma_mm,mu_mm,beta_mm,rho_mm,R_ratio,f_GR):
    
    #C = c_mm*(1-R^mu_mm)^gamma_mm
    C = np.zeros((2,2))
    C[:,0] = c_mm[:,0] #mode-mixity
    C[:,1] = c_mm[:,1]*((1.0-R_ratio**mu_mm[:,1])**gamma_mm[:,1]*1.0/f_GR)**(beta_mm[:,1]*(1-R_ratio)**rho_mm[:,1])
            
    
    return C

def paris_law_n(beta_mm,rho_mm,R_ratio):
    
    #n = beta_mm*(1-R)^rho_mm
    n = np.zeros((2,2))
    n[:,0] = beta_mm[:,0] #mode-mixity
    n[:,1] = beta_mm[:,1]*(1.0-R_ratio)**rho_mm[:,1]

    return n


def fatigue_mm_int(C,n,mmixity_b,mmixity_int):
    
    #two points, value not relevant
    GIc = 0.21 
    GIIc = 0.71    

    log_Gc_p1 = np.log10(np.array([GIc,GIc]))
    log_dadn_Gc_p1 = n[:,1]*log_Gc_p1[:]+np.log10(C[:,1])
    log_Gc_p2 = np.log10(np.array([GIIc,GIIc]))
    log_dadn_Gc_p2 = n[:,1]*log_Gc_p2[:]+np.log10(C[:,1])

    line0 = line2D()
    line1 = line2D()
    
    p1_l0 = np.array([log_Gc_p1[0],log_dadn_Gc_p1[0]])
    p2_l0 = np.array([log_Gc_p2[0],log_dadn_Gc_p2[0]])

    p1_l1 = np.array([log_Gc_p1[1],log_dadn_Gc_p1[1]])
    p2_l1 = np.array([log_Gc_p2[1],log_dadn_Gc_p2[1]])
    
    line0.line_two_points(p1_l0,p2_l0)
    line1.line_two_points(p1_l1,p2_l1)
    
    n0 = line0.n
    n1 = line1.n
   
    l_intersect,t0, c_linear = line0.intersect_w_line(line1)
    
    n0_norm = n0/np.linalg.norm(n0)
    n1_norm = n1/np.linalg.norm(n1)

    dtheta = np.arccos(np.clip(np.dot(n1_norm,n0_norm),-1.0,1.0))
    
    if  mmixity_b[1] == mmixity_b[0]:        
        dtheta_int = 0.0
        mm_factor = 0.0
    else:
        mm_factor = (mmixity_int-mmixity_b[0])/(mmixity_b[1]-mmixity_b[0])
        dtheta_int = dtheta*mm_factor
        cross_p_n = np.cross(n1_norm,n0_norm)
        if cross_p_n > 0:
            dtheta_int = - dtheta_int

    rot_a = np.array([[np.cos(dtheta_int),-np.sin(dtheta_int)],
                       [np.sin(dtheta_int),np.cos(dtheta_int)]])
    
    nint = np.dot(rot_a,n0_norm)
    line_int = line2D()
    line_int.n = nint.flatten()
    if len(l_intersect) != 0:
        #line_int.point = l_intersect
        line_int.point = l_intersect
    elif c_linear:
        line_int.point = p1_l0
        
    n_int = line_int.slope()
    log_C_int = line_int.y_axis_int()
    C_int =10**log_C_int 
    
    #BK interpolation

#    print(GC_int)
#    #linear interpolation
#    Gth_bounds = Gth[b_inds,1]
#    Gth_int = Gth_bounds[0] + (Gth_bounds[1]-Gth_bounds[0])*mm_factor
#    
#    
#    log_Gc_int = np.log10(GC_int)
#    log_dadnc_int = log_C_int+n_int*log_Gc_int
#    log_Gth_int =  np.log10(Gth_int)
#    log_dadngth_int = log_C_int+n_int*log_Gth_int
#    
#    point_2 = np.array([log_Gc_int,log_dadnc_int]).flatten()
#    point_1 = np.array([log_Gth_int,log_dadngth_int]).flatten()
#    
#    line_int.line_two_points(point_1,point_2)
#

    
    return C_int,n_int


def paris_law(GT_i,pl_C,pl_n):
    
    dadn = pl_C*(GT_i)**pl_n 

    return dadn

def dadn_mm_fatigue(G,vcct_el,dadn_mm_table,rcurve_mm_table,all_nodes_dict):

    #hardcoded for now
    table_no = int(all_nodes_dict[vcct_el]['dadn_T'])
    c_mm,beta_mm,mu_mm,gamma_mm,rho_mm = read_from_dadn_table(table_no,dadn_mm_table)
    R_ratio = inp_vars.R_ratio
    
    #calculate total ERR and shear ERR
    GT_max,GS_max = GTmax_GSmax(G)
    #calculate mode mixity
    mmixity = mode_mixity(GT_max,GS_max)
    #calculate indexes to the dadn_table that bound the mode-mixity
    mm_b_inds  = select_bounds(c_mm,mmixity)
    
    #replace by an R curve factor - assumed 1 for now, needs to be coded
    f_GR = 1
    
    #extracts Paris law coeficient and exponent bounding mode-mixity
    pl_C = paris_law_C(c_mm[mm_b_inds,:],gamma_mm[mm_b_inds,:],mu_mm[mm_b_inds,:],beta_mm[mm_b_inds,:],rho_mm[mm_b_inds,:],R_ratio,f_GR)
    pl_n = paris_law_n(beta_mm[mm_b_inds,:],rho_mm[mm_b_inds,:],R_ratio)

    #interpolates mode-mixity
    mmixity_b = c_mm[mm_b_inds,0] #values of mode_mixity corresponding to the bounds
    pl_C_int,pl_n_int = fatigue_mm_int(pl_C,pl_n,mmixity_b,mmixity)
    dadn = paris_law(GT_max,pl_C_int,pl_n_int) 
    
    return dadn
    
    
def update_pnewdt(pnewdt,vcct_els_front_dict,cut_back_factor,unst_gth_tol):
    
    for el, el_info in vcct_els_front_dict.items():
        if el_info['f_cr_s'] > 1.0+unst_gth_tol and el_info['ss_a'] == True:
            pnewdt = np.array([cut_back_factor])
            break
            
    return pnewdt      


def dl0_sample_points_asc(sample_points_xy,vcct_el,all_nodes_dict):
    
    vcct_el_coords = all_nodes_dict[vcct_el]['coords'][:2]
    dl0 = vcct_el_coords - sample_points_xy
    
    return np.linalg.norm(dl0,axis=1)

def dl0_sample_points(sample_points_xy,vcct_el,all_nodes_dict):
    
    vcct_el_coords = all_nodes_dict[vcct_el]['coords']
    dl0 = vcct_el_coords - sample_points_xy
    
    return np.linalg.norm(dl0,axis=1)

def select_sample(vcct_els_front_dict,option='max_d'):

    el_info_0 = list(vcct_els_front_dict.items())[0][1]
    list_vcct_els_front_keys = list(el_info_0.keys())
    list_vcct_els_front_keys.remove('ss_a')
    for el,el_info in vcct_els_front_dict.items():
        #print('el',el)
        if option == 'max_d':

            ind = np.argmax(np.divide(el_info['dl0_s'], el_info['dadn_s'], 
                                      out=np.zeros_like(el_info['dl0_s']), 
                                      where=el_info['dadn_s']!=0))
        elif option == 'max_dadn':
            ind = np.argmax(el_info['dadn_s'])          
        elif option == 'max_gmax':
            ind = np.argmax(np.sum(el_info['G_s'],axis=1))   
        
        for vcct_els_front_key in list_vcct_els_front_keys:
            #print('vcct_els_front_key',vcct_els_front_key)
            vcct_els_front_dict[el][vcct_els_front_key] = el_info[vcct_els_front_key][ind]
        
    
    return vcct_els_front_dict
# def select_sample_ind(dadn_s_cfront,G_s_cfront,dl0_s_cfront,option='max_d'):

    
#     if option == 'max_d':
#         s_max_index =[np.argmax(dl0_el/dadn_el) 
#                                 for dadn_el,dl0_el in zip(dadn_s_cfront,dl0_s_cfront)]
#     elif option == 'max_dadn':
#         s_max_index = [np.argmax(dadn_el) for dadn_el in dadn_s_cfront]
#     elif option == 'max_gmax':
#         s_max_index = [np.argmax(np.sum(G_s,axis=1)) for G_s in G_s_cfront]
#     else:
#         print('----error----')
            
#     return  s_max_index  

def read_step_inc(filename):
    step,inc = np.loadtxt(filename)
    return int(step),int(inc)


def steps_input():

    if hasattr(inp_vars,'steps_cycles'):
       steps_cycles = inp_vars.steps_cycles        
    else:
#       print('WARNING: number of cycles in each step not specified - A quasi-static analysis assumed')
       steps_cycles = None
       
    return steps_cycles

def fatigue_load_step(steps_cycles,step):
    
    try: 
        return steps_cycles[step-1] != 0
    except:
        return False 



def min_cycles_dN(vcct_els_front_dict,max_normalized_inc,all_nodes_dict,option='extrap'):
    
    dN_c = np.array([el_info['dl0_s']*max_normalized_inc/el_info['dadn_s'] if (el_info['dadn_s'] > 0 and el_info['ss_a'] == True) else 1.0E50
                     for el,el_info in vcct_els_front_dict.items()])
    
    dN = np.min(dN_c)
    
    if option == 'exact':
        no_els_front = len(vcct_els_front_dict)
        dN_c_inc_to_fail = np.array([]*no_els_front)
        for i,(el,el_info) in enumerate(vcct_els_front_dict.items()):
            d_to_fail = 1.0-all_nodes_dict[el]['d']
            c_inc_to_fail = el_info['dl0_s']*d_to_fail
            if el_info['dadn_s'] > 0.0:
                dN_c_inc_to_fail[i] = c_inc_to_fail/el_info['dadn_s'] 
            else:
                dN_c_inc_to_fail[i] = 1.0E50
                
        dN = min(dN,np.min(dN_c_inc_to_fail))

    return dN

def unloading_slope(vcct_el,all_nodes_dict):
    
    fmaxi = all_nodes_dict[vcct_el]['fmax']
    delta_maxi = all_nodes_dict[vcct_el]['delta_max']
    
    m_max = -10000.0
    m_arr = np.array([m_max]*3)
    
    for i,(f,delta) in enumerate(zip(fmaxi,delta_maxi)):
        if delta != 0.0:
            m_aux = -f/delta 
            if m_aux > m_max and (m_aux < 0):
                m_arr[i] = m_aux

    all_nodes_dict[vcct_el]['m'] = m_arr

      
    return all_nodes_dict

def cycles_vs_crack_length_to_file(filename,ncycles,xy_node_d,glob_dir='x'):
    
    if glob_dir == 'x':
        ind = 0
    elif glob_dir == 'y':
        ind = 1
    elif glob_dir == 'z':
        print('extraction along global z not admissible')
    
    pos_max = np.max(xy_node_d[:,ind])
    
    with open(filename,'a') as f_handle:
        f_handle.write(str(ncycles)+','+str(pos_max)+'\n')    

def cycles_to_file(filename,cycles):
    "outputs the total cycles at the end of each iteration to a file"
    with open(filename,'a') as f_handle:
        f_handle.write(str(cycles)+'\n')    

def read_cycles(filename):
    
    if os.path.isfile(filename) == True:
        cycles = np.loadtxt(filename,ndmin=1)
    else:
        with open(filename,'w') as f_handle:
            f_handle.write(str(0)+'\n')
        cycles = np.array([0])
        
    return cycles[-1]


def calc_softening_EXP(vcct_els_front_dict,vcct_els_edge_of_model_dict,n_min_dN,ntotal,PR,all_nodes_dict,option='extrap'):
	
    #fatigue:
    d_tolerance_split = 0.001 #tolerance for splitting the element
    if option == 'exact':
        exact_tol = 0.03
    else:
        exact_tol = 1E-4
     
    for el,el_info in vcct_els_front_dict.items():
        if el_info['ss_a'] == True:
           
            #calculates crack increment
            dc = el_info['dadn_s']*n_min_dN
            #updates next node being considered
            all_nodes_dict[el]['a_node_next_p'] = el_info['a_node_next_p']
            
            if el_info['dl0_s'] == 0:
                d_inc = 1.0 
                print('Warning: no crack front identified - node released!')
            else:
                d_inc = dc/el_info['dl0_s'] #computes increment in 'd'
    
            #element just started softening, record fmax and gmax
            if ((all_nodes_dict[el]['d'] < d_tolerance_split) and 
                (all_nodes_dict[el]['d'] +d_inc > d_tolerance_split)):
                    
                all_nodes_dict[el]['fmax'] = all_nodes_dict[el]['force']
                all_nodes_dict[el]['delta_max'] = el_info['delta_s']
                all_nodes_dict[el]['Gmax'] = np.sum(el_info['G_s'])
                
               # all_nodes_dict[el_info['a_node_next_p']]['a_node_prev'] = el
               # all_nodes_dict[el]['a_node_ref'] = el_info['a_node_ref']
                
                all_nodes_dict = unloading_slope(el,all_nodes_dict)
                #all_nodes_dict[el]['a_node_ref'] = el_info['a_node_ref']
                
            #increment damage
            if d_inc < 0.0:
                print('el',el)
                print('el_info_dadn',el_info['dadn_s'])
                
            all_nodes_dict[el]['d'] = all_nodes_dict[el]['d'] + d_inc
            
            delta_d =  all_nodes_dict[el]['d']  - 1.0
            if  delta_d  >=  - exact_tol:               
                if option == 'exact':
                    all_nodes_dict[el]['d'] = 1.0
                    all_nodes_dict[el]['ncycles_failed'] = ntotal+n_min_dN
                    
                #caps damage variable at one
                if option == 'extrap':
                    n_min_delta_d = 0.0
                    if delta_d > exact_tol:
                        n_min_delta_d =  delta_d*el_info['dl0_s']/el_info['dadn_s']
                    
                    #caps damage at 1.0                  
                    all_nodes_dict[el]['d'] = 1.0
                    #corrects for the fact that n_min_dN is calclulated based on the total d_inc
                    all_nodes_dict[el]['ncycles_failed'] = ntotal+(n_min_dN - n_min_delta_d)
                              
                    a_next_node_p = all_nodes_dict[el]['a_node_next_p']
                    #finds whether next node is already split 
                    if a_next_node_p not in vcct_els_edge_of_model_dict:
                        if all_nodes_dict[a_next_node_p]['d'] > d_tolerance_split:
                            print('WARNING: self similar crack growth may have been compromised')
                            print('all_nodes_dict[a_next_node_p]-d',all_nodes_dict[a_next_node_p]['d'])
                            all_nodes_dict[a_next_node_p]['d'] = all_nodes_dict[a_next_node_p]['d'] + delta_d
                        else:
                            #adds delta_d to next node since if larger than exact_tol, i.e. significant
                            if delta_d > exact_tol:
                                all_nodes_dict[a_next_node_p]['d'] += delta_d
                            
                            #checks next node's d variable and if larger than the split tolerance, save fmax delta_max and calculate unloading slope
                            if (all_nodes_dict[a_next_node_p]['d'] > d_tolerance_split):                       
                                all_nodes_dict[a_next_node_p]['fmax'] = all_nodes_dict[a_next_node_p]['force']
                                
                                all_nodes_dict[a_next_node_p]['delta_max'] = vcct_els_front_dict[a_next_node_p]['delta_s']
                                all_nodes_dict[a_next_node_p]['Gmax'] = np.sum(vcct_els_front_dict[a_next_node_p]['G_s'])
                                all_nodes_dict[a_next_node_p]['a_node_next_p'] = vcct_els_front_dict[a_next_node_p]['a_node_next_p']
                            
                                all_nodes_dict = unloading_slope(a_next_node_p,all_nodes_dict)
                            
                                if 'a_acc' in all_nodes_dict[el].keys():   
                                     all_nodes_dict[a_next_node_p]['a_acc'] = all_nodes_dict[el]['a_acc'] + el_info['dl0_s']
                                else:
                                     all_nodes_dict[a_next_node_p]['a_acc'] = el_info['dl0_s']
                
                # if option == 'extrap':                 
                #     a_next_node_p = all_nodes_dict[el]['a_node_next_p']
                #     print('el',el)
                #     if all_nodes_dict[a_next_node_p]['d'] != 0:
                #         print('WARNING: self similar crack growth may have been compromised')
                #         print('all_nodes_dict[a_next_node_p]-d',all_nodes_dict[a_next_node_p]['d'])
                #         all_nodes_dict[a_next_node_p]['d'] = all_nodes_dict[a_next_node_p]['d'] + (all_nodes_dict[el]['d']  - 1.0)
                #     else:
                #         n_min_delta_d = 0.0
                #         if delta_d > exact_tol:
                #             all_nodes_dict[a_next_node_p]['d'] = delta_d
                #             n_min_delta_d =  delta_d*el_info['dl0_s']/el_info['dadn_s']
                            
                #     all_nodes_dict[el]['d'] = 1.0
                #     all_nodes_dict[el]['ncycles_failed'] = ntotal+(n_min_dN - n_min_delta_d)
                #     if 'a_acc' in all_nodes_dict[el].keys():   
                #         all_nodes_dict[a_next_node_p]['a_acc'] = all_nodes_dict[el]['a_acc'] + el_info['dl0_s']
                #     else:
                #         all_nodes_dict[a_next_node_p]['a_acc'] = el_info['dl0_s']
                        
                #     if (all_nodes_dict[a_next_node_p]['d'] > d_tolerance_split):                       
                #         all_nodes_dict[a_next_node_p]['fmax'] = all_nodes_dict[a_next_node_p]['force']
                        
                #         all_nodes_dict[a_next_node_p]['delta_max'] = vcct_els_front_dict[a_next_node_p]['delta_s']
                #         all_nodes_dict[a_next_node_p]['Gmax'] = np.sum(vcct_els_front_dict[a_next_node_p]['G_s'])
                #         all_nodes_dict[a_next_node_p]['a_node_next_p'] = vcct_els_front_dict[a_next_node_p]['a_node_next_p']
                        
                        
                #         #all_nodes_dict[all_nodes_dict[a_next_node_p]['a_node_next_p']]['a_node_prev'] =  a_next_node_p
                #         #all_nodes_dict[a_next_node_p]['a_node_ref'] = el
                #         all_nodes_dict = unloading_slope(a_next_node_p,all_nodes_dict)
                        
                    
    return all_nodes_dict
    
     
def calc_softening_EXP_new(vcct_els_front_dict,vcct_els_edge_of_model_dict,n_min_dN,ntotal,PR,all_nodes_dict,option='extrap'):
	
    #fatigue:
    d_tolerance_split = 0.001 #tolerance for splitting the element
    if option == 'exact':
        exact_tol = 0.03
    else:
        exact_tol = 1E-4
            
    for el,el_info in vcct_els_front_dict.items():
        if el_info['ss_a'] == True:
           
            #calculates crack increment
            dc = el_info['dadn_s']*n_min_dN
            #updates next node being considered
            all_nodes_dict[el]['a_node_next_p'] = el_info['a_node_next_p']
            
            if el_info['dl0_s'] == 0:
                d_inc = 1.0 
                print('Warning: no crack front identified - node released!')
            else:
                d_inc = dc/el_info['dl0_s'] #computes increment in 'd'
    
            #element just started softening, record fmax and gmax
            if ((all_nodes_dict[el]['d'] < d_tolerance_split) and 
                (all_nodes_dict[el]['d'] +d_inc > d_tolerance_split)):
                    
                all_nodes_dict[el]['fmax'] = all_nodes_dict[el]['force']
                all_nodes_dict[el]['delta_max'] = el_info['delta_s']
                all_nodes_dict[el]['Gmax'] = np.sum(el_info['G_s'])
                
               # all_nodes_dict[el_info['a_node_next_p']]['a_node_prev'] = el
               # all_nodes_dict[el]['a_node_ref'] = el_info['a_node_ref']
                
                all_nodes_dict = unloading_slope(el,all_nodes_dict)
                #all_nodes_dict[el]['a_node_ref'] = el_info['a_node_ref']
                
            #increment damage
            if d_inc < 0.0:
                print('el',el)
                print('el_info_dadn',el_info['dadn_s'])
                
            all_nodes_dict[el]['d'] = all_nodes_dict[el]['d'] + d_inc
            
            delta_d =  all_nodes_dict[el]['d']  - 1.0
            if  delta_d  >=  - exact_tol:               
                if option == 'exact':
                    all_nodes_dict[el]['d'] = 1.0
                    all_nodes_dict[el]['ncycles_failed'] = ntotal+n_min_dN
                    
                #caps damage variable at one
                if option == 'extrap':
                    n_min_delta_d = 0.0
                    if delta_d > exact_tol:
                        n_min_delta_d =  delta_d*el_info['dl0_s']/el_info['dadn_s']
                    
                    #caps damage at 1.0                  
                    all_nodes_dict[el]['d'] = 1.0
                    #corrects for the fact that n_min_dN is calclulated based on the total d_inc
                    all_nodes_dict[el]['ncycles_failed'] = ntotal+(n_min_dN - n_min_delta_d)
                              
                    a_next_node_p = all_nodes_dict[el]['a_node_next_p']

                    if 'a_acc' in all_nodes_dict[el].keys():   
                        all_nodes_dict[a_next_node_p]['a_acc'] = all_nodes_dict[el]['a_acc'] + el_info['dl0_s']
                    else:
                        all_nodes_dict[a_next_node_p]['a_acc'] = el_info['dl0_s']
                                     
                    #finds whether next node is already split 
                    if a_next_node_p not in vcct_els_edge_of_model_dict:
                        if all_nodes_dict[a_next_node_p]['d'] > d_tolerance_split:
                            print('WARNING: self similar crack growth may have been compromised')
                            print('all_nodes_dict[a_next_node_p]-d',all_nodes_dict[a_next_node_p]['d'])
                            all_nodes_dict[a_next_node_p]['d'] = all_nodes_dict[a_next_node_p]['d'] + delta_d
                        else:
                            #adds delta_d to next node since if larger than exact_tol, i.e. significant
                            if delta_d > exact_tol:
                                all_nodes_dict[a_next_node_p]['d'] += delta_d
                            
                            #checks next node's d variable and if larger than the split tolerance, save fmax delta_max and calculate unloading slope
                            if (all_nodes_dict[a_next_node_p]['d'] > d_tolerance_split):                       
                                all_nodes_dict[a_next_node_p]['fmax'] = all_nodes_dict[a_next_node_p]['force']
                                
                                all_nodes_dict[a_next_node_p]['delta_max'] = vcct_els_front_dict[a_next_node_p]['delta_s']
                                all_nodes_dict[a_next_node_p]['Gmax'] = np.sum(vcct_els_front_dict[a_next_node_p]['G_s'])
                                all_nodes_dict[a_next_node_p]['a_node_next_p'] = vcct_els_front_dict[a_next_node_p]['a_node_next_p']
                            
                                all_nodes_dict = unloading_slope(a_next_node_p,all_nodes_dict)
                
                # if option == 'extrap':                 
                #     a_next_node_p = all_nodes_dict[el]['a_node_next_p']
                #     print('el',el)
                #     if all_nodes_dict[a_next_node_p]['d'] != 0:
                #         print('WARNING: self similar crack growth may have been compromised')
                #         print('all_nodes_dict[a_next_node_p]-d',all_nodes_dict[a_next_node_p]['d'])
                #         all_nodes_dict[a_next_node_p]['d'] = all_nodes_dict[a_next_node_p]['d'] + (all_nodes_dict[el]['d']  - 1.0)
                #     else:
                #         n_min_delta_d = 0.0
                #         if delta_d > exact_tol:
                #             all_nodes_dict[a_next_node_p]['d'] = delta_d
                #             n_min_delta_d =  delta_d*el_info['dl0_s']/el_info['dadn_s']
                            
                #     all_nodes_dict[el]['d'] = 1.0
                #     all_nodes_dict[el]['ncycles_failed'] = ntotal+(n_min_dN - n_min_delta_d)
                #     if 'a_acc' in all_nodes_dict[el].keys():   
                #         all_nodes_dict[a_next_node_p]['a_acc'] = all_nodes_dict[el]['a_acc'] + el_info['dl0_s']
                #     else:
                #         all_nodes_dict[a_next_node_p]['a_acc'] = el_info['dl0_s']
                        
                #     if (all_nodes_dict[a_next_node_p]['d'] > d_tolerance_split):                       
                #         all_nodes_dict[a_next_node_p]['fmax'] = all_nodes_dict[a_next_node_p]['force']
                        
                #         all_nodes_dict[a_next_node_p]['delta_max'] = vcct_els_front_dict[a_next_node_p]['delta_s']
                #         all_nodes_dict[a_next_node_p]['Gmax'] = np.sum(vcct_els_front_dict[a_next_node_p]['G_s'])
                #         all_nodes_dict[a_next_node_p]['a_node_next_p'] = vcct_els_front_dict[a_next_node_p]['a_node_next_p']
                        
                        
                #         #all_nodes_dict[all_nodes_dict[a_next_node_p]['a_node_next_p']]['a_node_prev'] =  a_next_node_p
                #         #all_nodes_dict[a_next_node_p]['a_node_ref'] = el
                #         all_nodes_dict = unloading_slope(a_next_node_p,all_nodes_dict)
                               
                    
    return all_nodes_dict


def select_based_on_index(list_of_val_lists,s_max_index):
    
    list_arrs = []
    for list_of_val in list_of_val_lists:
        list_arrs += [np.array([val_s[ind] for ind,val_s in zip(s_max_index,list_of_val)])]
        
    return list_arrs
    

def onset_crit(vcct_el_dict,all_els_dict,all_nodes_dict,onset_node_list, inc,crit = 'quads'):

    for i,node in enumerate(onset_node_list):
        
        # area_onset = calc_area_carvalho7(node,ind_ref_adj_el_dummy,vcct_el_dict,
        #                             all_nodes_dict,all_els_dict)
        #vcct_el_coords = all_nodes_dict[node]['coords']
        #vcct_adj_els = vcct_el_dict[node]
        
#        area_onset_old = calc_area_VCCT_el_pristine(vcct_el_coords,vcct_adj_els,all_nodes_dict)

        area_onset = calc_area_VCCT_el_pristine_v2(node,vcct_el_dict,all_nodes_dict)
        
        stress_v = all_nodes_dict[node]['force_onset']/area_onset
        
        Z_S = all_nodes_dict[node]['Z_S']
        Z_T = all_nodes_dict[node]['Z_T']
        
        #print('Z_S,Z_S',Z_S,Z_T)
        if crit == 'quads':
            #print('stress_v',stress_v)
            f_crt_s = (stress_v[0]/Z_S)**2.0+(stress_v[1]/Z_S)**2.0+(stress_v[2]/Z_T)**2.0

#            if ((node == 25915) or (node == 25412)):
#                print('node',node)
#                print('area_onset_old',area_onset_old)
#                print('area_onset',area_onset)
#                print('vcct_el_dict',vcct_adj_els)
#                print('stress_v',stress_v)
#                print('area_onset',area_onset)
#                print('z_s,z_t',Z_S,Z_T)
        
        all_nodes_dict[node]['f_crt_stress'] = f_crt_s
        
              
    return all_nodes_dict


def hot_spot(vcct_el_dict,all_nodes_dict,node):
    
    nodes_a = vcct_el_dict[node]
    el_type = vcct_el_type(node,vcct_el_dict)
    
    f_crt_ref = all_nodes_dict[node]['f_crt_stress']
    
    tol = 1.0E-3
 
    if el_type == 'c':
        antena_nodes_ind = np.array([0,2])
    elif el_type == 'e':
        antena_nodes_ind = np.array([0,2,4])
    elif el_type == 'm':
        antena_nodes_ind = np.array([0,2,4,6])
    
    f_crt_a = np.zeros(len(antena_nodes_ind))
    for i,ind in enumerate(antena_nodes_ind):
        if 'f_crt_stress' in all_nodes_dict[nodes_a[ind]].keys():
            #print('nodes_a[ind]',nodes_a[ind])
            f_crt_a[i] = all_nodes_dict[nodes_a[ind]]['f_crt_stress']

    #equal sign prevents erroneous results whene stress is uniform
    equal_bool = np.all(abs(f_crt_ref - f_crt_a) < tol) #plateau
    peak_bool = np.all(f_crt_ref > f_crt_a) #peak
    
    bool_hot_spot =  (equal_bool or peak_bool) and (f_crt_ref > 1.0)
    
    #print('node,hot_spot',node,bool_hot_spot)
    
    return bool_hot_spot

def onset_cracks(vcct_el_dict,all_nodes_dict,onset_node_list,stress_crit_tol):
    f_node_list = []
    
    #first loop identifies nodes that should be failed:
    for node in onset_node_list:
        
        #first criterion - no failed/damaged antena node
        bool_pristine_a_nodes = vcct_el_any_antena_status_bool(node,vcct_el_dict,all_nodes_dict,'dgt0')
        if (bool_pristine_a_nodes == False):  
            #second criterion - highest stress of all antena nodes = hot_spot
            bool_hot_spot = hot_spot(vcct_el_dict,all_nodes_dict,node)
            if bool_hot_spot:
                f_node_list += [node]
    
    #second loop fails nodes:
    #Failing nodes progressively affects the failure pattern due to the criterion 
    #of not failing nodes that have failed antena nodes - code potentially be coded differently
    if len(f_node_list) != 0:
        f_node_list_fcrt_norm = np.array([all_nodes_dict[node]['f_crt_stress'] for node in f_node_list])
        fcrt_max = np.max(f_node_list_fcrt_norm)
        f_node_list_fcrt_norm = f_node_list_fcrt_norm/fcrt_max
        for i,f_node in enumerate(f_node_list):
            if (f_node_list_fcrt_norm[i] > 1.0 - stress_crit_tol):
                all_nodes_dict[f_node]['d'] = 1.0
                print('node onset',f_node)
                for a_node in vcct_el_dict[f_node]:
                    if a_node != 0:
                        all_nodes_dict[a_node]['d'] = 1.0
                    #        print('a_node_onset',a_node)
    
    return all_nodes_dict

def update_edge_of_model(all_nodes_dict,vcct_els_edge_of_model_dict):
    
    for edge_node,ref_node_info in vcct_els_edge_of_model_dict.items():
        ref_node = ref_node_info['ref_node']     
        if all_nodes_dict[ref_node]['d'] == 1.0:
           all_nodes_dict[edge_node]['d'] = 1.0
    
    return all_nodes_dict


def island_remover(vcct_el_dict,all_nodes_dict):
    
    island_detected = True
    max_iterations = 10
    for i in range(max_iterations):
        island_detected = False
        for el,adj_nodes in vcct_el_dict.items():
            if all_nodes_dict[el]['d'] != 1.0:
                
                el_type = vcct_el_type(el,vcct_el_dict)
                
                if el_type == 'c':
                    antena_nodes_ind = np.array([0,2])
                elif el_type == 'e':
                    antena_nodes_ind = np.array([0,2,4])
                elif el_type == 'm':
                    antena_nodes_ind = np.array([0,2,4,6])
                
                d_antena_nodes = np.zeros(len(antena_nodes_ind))
                for i,ind in enumerate(antena_nodes_ind):
                    d_antena_nodes[i] = all_nodes_dict[adj_nodes[ind]]['d']
                if np.all(d_antena_nodes > 0.0):
                    all_nodes_dict[el]['d'] = 1.0
                    print('island resolved',el)
                    island_detected = True
        if island_detected == False:
            break
    if i == max_iterations:
        print('WARNING: model is unzipping - check solution')
            
    return all_nodes_dict

def no_fail_zone(vcct_els_front_dict,all_nodes_dict):
    no_fail_node_list = []
    
    for node,vcct_el_info in vcct_els_front_dict.items():
        no_fail_node_list += all_nodes_dict[node]['nodes_in_cube']
        
    no_fail_node_arr = np.array(no_fail_node_list)
    no_fail_node_arr = np.unique(no_fail_node_arr)
    
    return no_fail_node_arr    

def remove_nodes_in_no_fail_zone(onset_nodes_arr,no_fail_nodes_arr):
    
    return np.setdiff1d(onset_nodes_arr,no_fail_nodes_arr)


def VCCT_main():

    #determines if cycles have been specified
    steps_cycles = steps_input()  
    
    #reads current step number
    step,inc = read_step_inc(current_dir+'STEP.txt')

    #print('############################################')
    print('increment:', inc)

    #if inc == 18:
    #    sys.exit("exit for debuging")
        
    #load dadn and r_curve definitions:
    dadn_mm_table,rcurve_mm_table = load_dadn_r_curve_data(current_dir)
        
    #checks if the current step is a fatigue step:
    bool_fatigue_f = fatigue_load_step(steps_cycles,step)
    ntotal = 0
    n_min_dN = 0
    if bool_fatigue_f:
        ntotal = read_cycles(current_dir+'ncycles.txt')
    #default:
    pnewdt = np.array([1.0E6])

    list_of_names = ['analysis_dict','all_els_dict','all_nodes_dict','node_2vcct_adj_els_dict','vcct_el_dict']
    [analysis_dict,all_els_dict,all_nodes_dict,node_2vcct_adj_els_dict,vcct_el_dict] = load_vcct_data(current_dir+input_file,list_of_names)

    #options:
    PR = analysis_dict['PR']
    max_normalized_inc = analysis_dict['max_normalized_inc']
    unst_gth_tol = analysis_dict['unst_gth_tol']
    cut_back_factor = analysis_dict['cut_back_factor']
    static_scale_f = analysis_dict['static_scale_f']
    stress_crit_tol = analysis_dict['stress_crit_tol']
    
    max_normalized_inc = max_normalized_inc_PR(PR,max_normalized_inc)
    
    #debuf nonzero d values
    for i,info in all_nodes_dict.items():
        if info['d'] < 0.0:
            print('element d < 0:',i)
    
    
    all_nodes_dict,onset_nodes_arr = load_FE_data(all_els_dict,all_nodes_dict)
        
    node_set_critical = []
    vcct_els_front_dict,vcct_els_edge_of_model_dict = VCCT_el_G_samples(all_els_dict,all_nodes_dict,node_2vcct_adj_els_dict,vcct_el_dict,analysis_dict)
    vcct_els_front_dict = self_similar_active(vcct_els_front_dict,vcct_el_dict,all_nodes_dict)
    if len(vcct_els_front_dict) > 0: #crack front active
         vcct_els_front_dict,bool_static_f = f_cr_static_cfront(rcurve_mm_table,vcct_els_front_dict,all_nodes_dict,analysis_dict)
         if bool_static_f:
             vcct_els_front_dict = dadn_static_sample(vcct_els_front_dict,static_scale_f)
         elif bool_fatigue_f:
                 #needs to be coded:
                 vcct_els_front_dict = dadn_fatigue_sample(vcct_els_front_dict,dadn_mm_table,rcurve_mm_table,all_nodes_dict)
                 #initializes fatigue increment cycle count
         elif verbose_vcct_el_inc:
             #compute ERRs etc for output
             vcct_els_front_dict = dadn_static_sample(vcct_els_front_dict,static_scale_f)
         
         if bool_static_f or bool_fatigue_f or verbose_vcct_el_inc:
             if analysis_dict['lead_nodes_bool'] == True: 
                 vcct_els_front_dict = select_sample(vcct_els_front_dict,option = 'max_dadn')
             else:
                 vcct_els_front_dict = select_sample(vcct_els_front_dict,option = 'max_d')
                 
         if bool_static_f or bool_fatigue_f:
             #vcct_els_front_dict = self_similar_active(vcct_els_front_dict,vcct_el_dict,all_nodes_dict)
            
    #         #calculated the number of cycles
             n_min_dN = min_cycles_dN(vcct_els_front_dict,max_normalized_inc,all_nodes_dict,option='extrap')
             
             all_nodes_dict = calc_softening_EXP_new(vcct_els_front_dict,vcct_els_edge_of_model_dict,
                                                 n_min_dN,ntotal,PR,all_nodes_dict,option = 'extrap')
             
             
             all_nodes_dict = update_edge_of_model(all_nodes_dict,vcct_els_edge_of_model_dict)
             if bool_fatigue_f:
                 cycles_to_file(current_dir + 'ncycles.txt',ntotal+n_min_dN)
                 
             if bool_static_f:
                 pnewdt = update_pnewdt(pnewdt,vcct_els_front_dict,cut_back_factor,unst_gth_tol)

    #if len(onset_nodes_arr) != 0:
    #    sys.exit("onset possible - exit for debuging")
    
    
    # if ctod_bool == True:         
    #     no_els = len(vcct_els_front_dict)
    #     el_number = np.zeros(no_els)
    #     el_count = 0
    #     for el,el_info in vcct_els_front_dict.items():
    #         if vcct_el_any_antena_status_bool(el,vcct_el_dict,all_nodes_dict,option='deq1'):
    #             el_number[el_count] = el
    #             el_count += 1
    #     el_number = el_number[:el_count]
    #     nodes_coords_ctod_arr,f_nodes = ctod_calc(el_number,vcct_el_dict,all_nodes_dict)
    #     filename = current_dir + 'ctod'+str(int(step))+'_'+str(int(inc))+'.txt'
    #     with open(filename,'ab') as f_handle:
    #         np.savetxt(f_handle,nodes_coords_ctod_arr,fmt='%f') 
            
        # filename = current_dir + 'f_nodes_debug'+str(int(step))+'_'+str(int(inc))+'.txt'
        # with open(filename,'ab') as f_handle:
        #     np.savetxt(f_handle,f_nodes,fmt='%f') 


    
    if hasattr(inp_vars,'no_fail_zone_char_length'):
        #all_nodes_dict = reset_no_fail_zones(all_nodes_dict)
        
        #all_nodes_dict = setup_no_fail_zones(all_nodes_dict)
    #collects nodes in the user defined no fail zone
        no_fail_nodes_arr = no_fail_zone(vcct_els_front_dict,all_nodes_dict) 
        #print('no_fail_nodes_arr',no_fail_nodes_arr)
#    print('no_fail_nodes_arr:',no_fail_nodes_arr)
    #removes nodes in no fail zone from the user defined no fail zone
#    print('onset_nodes_arr:',onset_nodes_arr)
        # fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(9, 9))
        onset_nodes_arr = remove_nodes_in_no_fail_zone(onset_nodes_arr,no_fail_nodes_arr)
        # for node,nodes_info in all_nodes_dict.items():
        #     coords = all_nodes_dict[int(node)]['coords']
        #     plot_point(coords,[0,1],ax1,'b')
        # for node in no_fail_nodes_arr:
        #     coords = all_nodes_dict[int(node)]['coords']
        #     plot_point(coords,[0,1],ax1)
        # for el,el_info in vcct_els_front_dict.items():
        #     coords = all_nodes_dict[int(el)]['coords']
        #     plot_point(coords,[0,1],ax1,'r')
#    print('onset_nodes_arr - no_fail_nodes_arr:',onset_nodes_arr)
    #assesses onset criterion
    all_nodes_dict = onset_crit(vcct_el_dict,all_els_dict,all_nodes_dict,onset_nodes_arr,inc)
    #onsets cracks 
    all_nodes_dict = onset_cracks(vcct_el_dict,all_nodes_dict,onset_nodes_arr,stress_crit_tol)               

    #island removal -  enable crack merging
    all_nodes_dict = island_remover(vcct_el_dict,all_nodes_dict)
    
    #selects VCCT elements at the delamination front
    vcct_els_cfront = vcct_els_crackfront(vcct_el_dict,all_nodes_dict,option_node = 'deq1',option_a_nodes = 'dgt0')
    
    # print('vcct_els_cfront',vcct_els_cfront)
    # if analysis_dict['ctod_bool'] == True:
    #     vcct_els_cfront = expand_cfront_for_ctod_extraction(vcct_els_cfront,analysis_dict['ref_vect_ctod'],analysis_dict['nsearches_ctod'],vcct_el_dict,all_nodes_dict)
        
    #maps nodal infor to user elements 
    tip_els_I = map_vcctnodes2els(vcct_els_cfront,node_2vcct_adj_els_dict,all_nodes_dict,all_els_dict)
    write_to_file_els(current_dir+'TIP_ELS.txt',tip_els_I)

    write_pnewdt(current_dir+'pnewdt.txt',pnewdt)
    #print('verbose_vcct_el_inc',verbose_vcct_el_inc)
    if verbose_vcct_el_inc == True:
        if asc=='off':
            verbose_output_vcct_txt(step,inc,vcct_els_front_dict,vcct_el_dict,all_nodes_dict,ntotal+n_min_dN,bool_fatigue_f,analysis_dict)
        else:
            verbose_output_vcct_txt_asc(step,inc,vcct_els_front_dict,vcct_el_dict,all_nodes_dict,ntotal+n_min_dN,bool_fatigue_f)

    #updates nodal vcct data:
    #resets forces and displacements
    all_nodes_dict = reset_force_delta(vcct_els_cfront,onset_nodes_arr,all_nodes_dict,node_2vcct_adj_els_dict,all_els_dict)
    list_of_dicts = [all_nodes_dict]
    list_of_names = ['all_nodes_dict']
    save_vcct_data(current_dir+input_file,list_of_dicts,list_of_names)
    
    #return vcct_els_front_dict,vcct_el_dict,all_els_dict,all_nodes_dict,node_2vcct_adj_els_dict
    return


# def along_nearest_mesh_line(vcct_el,ns2crack,a_node_prist_along_normal,vcct_el_dict,all_nodes_dict):    

#     el_type = vcct_el_type(vcct_el,vcct_el_dict)
#     if el_type == 'm':
#         if len(a_node_prist_along_normal) == 0:
#             if all_nodes_dict[vcct_el]['d'] > 0.001:
#                 a_node_next = all_nodes_dict[vcct_el]['a_node_next_p']
#             else:
#                 no_antena_nodes = len(a_nodes_indexes)
#                 a_node_pristine = np.zeros(no_antena_nodes)
#                 vcct_el_nodes = vcct_el_dict[vcct_el]
#                 #order the antena nodes by direction: closest to the normal considered
#                 for j,a_node in enumerate(vcct_el_nodes[a_nodes_indexes]):
#                     if a_node != 0:
#                         if all_nodes_dict[a_node]['d'] < 0.001:
#                             a_node_next = a_node
#                             break
#         else:
#             a_node_next = a_node_prist_along_normal[0]
#         ns2crack = all_nodes_dict[a_node_next]['coords'][:2]  - all_nodes_dict[vcct_el]['coords'][:2]
#         ns2crack = np.array([ns2crack/np.linalg.norm(ns2crack)])
#         sample_xy = np.array([2.0*all_nodes_dict[vcct_el]['coords'][:2]-all_nodes_dict[a_node_next]['coords'][:2] ])  
#         a_node_prist_along_normal = np.array([a_node_next])

#     return sample_xy,ns2crack,a_node_prist_along_normal
    

def verbose_output_vcct_txt_asc(step,inc,vcct_els_front_dict,vcct_el_dict,all_nodes_dict,ntotal,bool_fatigue_f,glob_dir = 'x'):
    
    no_els = len(vcct_els_front_dict)
    G_s = np.zeros((no_els,3))
    xy_node = np.zeros((no_els,2))
    xy_node_d = np.zeros((no_els,2))
    d_val = np.zeros(no_els)
    el_number = np.zeros(no_els)
    f_cr = np.zeros(no_els)
    next_a_node = np.zeros(no_els)
    a_acc = np.zeros(no_els)
    vcct_el2txt = np.zeros((no_els,16))
    cf_v = np.zeros((no_els,4))
    el_count = 0
    for el,el_info in vcct_els_front_dict.items():

        if  vcct_el_any_antena_status_bool(el,vcct_el_dict,all_nodes_dict,option='deq1'):
            el_number[el_count] = el
            G_s[el_count,:] = el_info['G_s']
            xy_node[el_count,:] = all_nodes_dict[el]['coords'][:2]
            d_val[el_count] = all_nodes_dict[el]['d'] 
            for i,v in enumerate( all_nodes_dict[el]['cf_v']):
                cf_v[el_count,2*i:2*(i+1)] = v
    
            if d_val[el_count] > 0.001:
                next_antenna_xy = all_nodes_dict[all_nodes_dict[el]['a_node_next_p']]['coords'][:2]
            else:
                next_antenna_xy = xy_node[el_count,:]
            xy_node_d[el_count,:] = xy_node[el_count,:] + (next_antenna_xy-xy_node[el_count,:])*d_val[el_count]
            f_cr[el_count] = el_info['f_cr_s']
            if 'a_node_next_p' in all_nodes_dict[el]:
                next_a_node[el_count] = all_nodes_dict[el]['a_node_next_p']
            if 'a_acc' in all_nodes_dict[el]:
                a_acc[el_count] =  all_nodes_dict[el]['a_acc']
            el_count +=1
            
            
    vcct_el2txt[:el_count,0] = el_number[:el_count]
    vcct_el2txt[:el_count,1:3] = xy_node[:el_count] 
    vcct_el2txt[:el_count,3:5] = xy_node_d[:el_count] 
    vcct_el2txt[:el_count,5] = d_val[:el_count]

    vcct_el2txt[:el_count,6:9] = G_s[:el_count]
    vcct_el2txt[:el_count,9] = f_cr[:el_count]
    vcct_el2txt[:el_count,10] = next_a_node[:el_count]
    vcct_el2txt[:el_count,11:15] =  cf_v[:el_count,:]
    vcct_el2txt[:el_count,15] = a_acc[:el_count]
    

    vcct_el2txt = vcct_el2txt[:el_count,:] 
#    print('vcct_el2txt',vcct_el2txt)
    filename = current_dir + 'GT'+str(int(step))+'_'+str(int(inc))+'.txt'
    with open(filename,'ab') as f_handle:
        np.savetxt(f_handle,vcct_el2txt,fmt='%f') 

    if bool_fatigue_f:
        filename = current_dir + 'amax_vs_n.txt'
        cycles_vs_crack_length_to_file(filename,ntotal,xy_node_d[:el_count,:],glob_dir)
        

    return


def verbose_output_vcct_txt(step,inc,vcct_els_front_dict,vcct_el_dict,all_nodes_dict,ntotal,bool_fatigue_f,analysis_dict,glob_dir = 'x',):
    
    no_els = len(vcct_els_front_dict)
    G_s = np.zeros((no_els,3))
    xy_node = np.zeros((no_els,3))
    xy_node_d = np.zeros((no_els,3))
    d_val = np.zeros(no_els)
    el_number = np.zeros(no_els)
    f_cr = np.zeros(no_els)
    next_a_node = np.zeros(no_els)
    a_acc = np.zeros(no_els)
    vcct_el2txt = np.zeros((no_els,21))
    cf_v = np.zeros((no_els,6))
    el_count = 0
    for el,el_info in vcct_els_front_dict.items():

        if  vcct_el_any_antena_status_bool(el,vcct_el_dict,all_nodes_dict,option='deq1'):
            el_number[el_count] = el
            G_s[el_count,:] = el_info['G_s']
            xy_node[el_count,:] = all_nodes_dict[el]['coords']
            d_val[el_count] = all_nodes_dict[el]['d'] 
            for i,v in enumerate( all_nodes_dict[el]['cf_v']):
                cf_v[el_count,3*i:3*(i+1)] = v
    
            if d_val[el_count] > 0.001:
                next_antenna_xy = all_nodes_dict[all_nodes_dict[el]['a_node_next_p']]['coords']
            else:
                next_antenna_xy = xy_node[el_count,:]
            xy_node_d[el_count,:] = xy_node[el_count,:] + (next_antenna_xy-xy_node[el_count,:])*d_val[el_count]
            f_cr[el_count] = el_info['f_cr_s']
            if 'a_node_next_p' in all_nodes_dict[el]:
                next_a_node[el_count] = all_nodes_dict[el]['a_node_next_p']
            if 'a_acc' in all_nodes_dict[el]:
                a_acc[el_count] =  all_nodes_dict[el]['a_acc']
            el_count +=1
            
            
    vcct_el2txt[:el_count,0] = el_number[:el_count]
    vcct_el2txt[:el_count,1:4] = xy_node[:el_count] 
    vcct_el2txt[:el_count,4:7] = xy_node_d[:el_count] 
    vcct_el2txt[:el_count,7] = d_val[:el_count]

    vcct_el2txt[:el_count,8:11] = G_s[:el_count]
    vcct_el2txt[:el_count,11] = f_cr[:el_count]
    vcct_el2txt[:el_count,12] = next_a_node[:el_count]
    vcct_el2txt[:el_count,13:19] =  cf_v[:el_count,:]
    vcct_el2txt[:el_count,19] = a_acc[:el_count]
    

    vcct_el2txt = vcct_el2txt[:el_count,:] 
#    print('vcct_el2txt',vcct_el2txt)
    filename = current_dir + 'GT'+str(int(step))+'_'+str(int(inc))+'.txt'
    with open(filename,'ab') as f_handle:
        np.savetxt(f_handle,vcct_el2txt,fmt='%f') 

    if bool_fatigue_f:
        filename = current_dir + 'amax_vs_n.txt'
        cycles_vs_crack_length_to_file(filename,ntotal,xy_node_d[:el_count,:],glob_dir)
        
    if analysis_dict['ctod_bool'] == True:
        nsearches_ctod = analysis_dict['nsearches_ctod']
        ref_vect = analysis_dict['ref_vect_ctod']
        nodes_coords_ctod_arr = ctod_calc(el_number[:el_count],ref_vect,nsearches_ctod,vcct_el_dict,all_nodes_dict)
        filename = current_dir + 'ctod'+str(int(step))+'_'+str(int(inc))+'.txt'
        with open(filename,'ab') as f_handle:
            np.savetxt(f_handle,nodes_coords_ctod_arr,fmt='%f') 
            
    return


# def finds_reference_node_ind(vcct_el,vcct_el_dict,all_nodes_dict):

#     #extracts adjacent vcct els (nodes)
#     a_nodes_indexes = [0,2,4,6]
#     vcct_el_antena_nodes = vcct_el_dict[vcct_el][a_nodes_indexes]    
    
#     if vcct_el == 10608:
#         print('here')
#     if 'a_node_prev' in all_nodes_dict[vcct_el].keys():
#         for i,node in enumerate(vcct_el_antena_nodes):   
#             if node == all_nodes_dict[vcct_el]['a_node_prev']:
#                   ind_a_f = a_nodes_indexes[i] 
    
#     else:
#         ind_a_f_arr = np.array([],dtype=int)
#         delta_a_f_list =[]
#         d_arr = np.zeros(4)
#         for i,node in enumerate(vcct_el_antena_nodes):
#             if node != 0:
#                 d_arr[i] = all_nodes_dict[node]['d']
#                 if all_nodes_dict[node]['d'] == 1.0:
#                     ind_a_f_arr = np.append(ind_a_f_arr,i)
#                     delta_a_f_list +=[np.linalg.norm(all_nodes_dict[node]['delta'])]
                    
#         if len(ind_a_f_arr) == 0.0:
#             ind_a_f = np.argmax(d_arr)
#         if len(ind_a_f_arr) == 1:
#             ind_a_f = ind_a_f_arr[0]
#         if len(ind_a_f_arr) == 2:
#             if delta_a_f_list[0] > delta_a_f_list[1]:
#                 ind_a_f = ind_a_f_arr[0]
#             else:
#                 ind_a_f = ind_a_f_arr[1]
#         if len(ind_a_f_arr) == 3:
#             if np.all(ind_a_f_arr == np.array([0,2,3])):
#                 ind_a_f = 3
#             elif np.all(ind_a_f_arr == np.array([0,1,3])):
#                 ind_a_f = 0
#             elif np.all(ind_a_f_arr == np.array([0,1,2])):
#                 ind_a_f = 1
#             elif np.all(ind_a_f_arr == np.array([1,2,3])):
#                 ind_a_f = 2
#     print(ind_a_f)  
#     return ind_a_f


def next_a_node_ind_by_inc(a_nodes_inds,a_node_ind,inc):
    
    ind = np.nonzero(np.array(a_nodes_inds) == a_node_ind)[0][0]
    
    no_a_nodes = len(a_nodes_inds)
    ind_inc = index_roll(ind+inc,no_a_nodes)
    
    return a_nodes_inds[ind_inc]
    

def finds_reference_node_ind(vcct_el,vcct_el_coords,el_type,vcct_el_dict,
                             all_nodes_dict):    

    vcct_el_antena_nodes = vcct_el_dict[vcct_el][a_nodes_indexes]    
    if 'a_node_ref' in all_nodes_dict[vcct_el].keys():
        for i,node in enumerate(vcct_el_antena_nodes):   
            if node == all_nodes_dict[vcct_el]['a_node_ref']:
                  ind_a_f = i 
    
    elif 'a_node_prev' in all_nodes_dict[vcct_el].keys():
        for i,node in enumerate(vcct_el_antena_nodes):   
            if node == all_nodes_dict[vcct_el]['a_node_prev']:
                  ind_a_f = i 
    else:    
        if el_type == 'e':
            #original:
            # for i,node in enumerate(vcct_el_antena_nodes):
            #     if node == 0:
            #         ind_1 = index_roll(i+1,4)
            #         ind_2 = index_roll(i-1,4)
            #         if (all_nodes_dict[vcct_el_antena_nodes[ind_1]]['d'] > 
            #             all_nodes_dict[vcct_el_antena_nodes[ind_2]]['d']):
            #             ind_a_f = ind_1
            #         else:
            #             ind_a_f = ind_2
            
            
            rank_arr = np.zeros((3,3))
            for i,node in enumerate(vcct_el_antena_nodes[:3]):
            #builds a matrix with 'd' and 'delta_t/dl0'
            #first criterion is 'd' second 'delta_t/dl0'
                if node != 0:
                    node_coords = all_nodes_dict[node]['coords']
                    dl0 = np.linalg.norm(vcct_el_coords - node_coords)
                    rank_arr[i,0] = all_nodes_dict[node]['d']
                    rank_arr[i,1] = np.linalg.norm(all_nodes_dict[node]['delta'])/dl0
                    rank_arr[i,2] = i
            rank_arr = rank_arr[rank_arr[:,0].argsort()]
            rank_arr = rank_arr[rank_arr[:,1].argsort(kind='mergesort')]    
            ind_a_f =  int(rank_arr[-1,2])
#            print('ind_a_f',ind_a_f)
            a_nodes_indexes[ind_a_f]
            
        else:
            rank_arr = np.zeros((4,3))
            for i,node in enumerate(vcct_el_antena_nodes):
            #builds a matrix with 'd' and 'delta_t/dl0'
            #first criterion is 'd' second 'delta_t/dl0'
                if node != 0:
                    node_coords = all_nodes_dict[node]['coords']
                    dl0 = np.linalg.norm(vcct_el_coords - node_coords)
                    rank_arr[i,0] = all_nodes_dict[node]['d']
                    rank_arr[i,1] = np.linalg.norm(all_nodes_dict[node]['delta'])/dl0
                    rank_arr[i,2] = i
            rank_arr = rank_arr[rank_arr[:,0].argsort()]
            rank_arr = rank_arr[rank_arr[:,1].argsort(kind='mergesort')]    
            ind_a_f =  int(rank_arr[-1,2])

    return a_nodes_indexes[ind_a_f]

def load_FE_data(all_els_dict,all_nodes_dict):
    
    if test_main == False:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")   
            el_delta = np.loadtxt(current_dir+'ELDELTA.txt',ndmin=2) 
            el_force = np.loadtxt(current_dir+'ELFORCE.txt',ndmin=2) 
            el_force_onset = np.loadtxt(current_dir+'ELONSET.txt',ndmin=2)  
            all_nodes_dict,onset_node_list = load_el_I_2nodes(el_force,el_force_onset,el_delta,all_nodes_dict,all_els_dict)    
    else:
        all_nodes_dict = tn.test_upd_node_data(all_nodes_dict)
        
    return all_nodes_dict,onset_node_list


def finds_reference_node_ind_asc(vcct_el,vcct_el_dict,all_nodes_dict):    

    #extracts adjacent vcct els (nodes)
    el_type = vcct_el_type(vcct_el,vcct_el_dict)
        
    a_nodes_indexes = [0,2,4,6]
    vcct_el_antena_nodes = vcct_el_dict[vcct_el][a_nodes_indexes]    
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]
    if 'a_node_ref' in all_nodes_dict[vcct_el].keys():
        for i,node in enumerate(vcct_el_antena_nodes):   
            if node == all_nodes_dict[vcct_el]['a_node_ref']:
                  ind_a_f = i 
    
    elif 'a_node_prev' in all_nodes_dict[vcct_el].keys():
        for i,node in enumerate(vcct_el_antena_nodes):   
            if node == all_nodes_dict[vcct_el]['a_node_prev']:
                  ind_a_f = i 
    else:    
        if el_type == 'e':
            for i,node in enumerate(vcct_el_antena_nodes):
                if node == 0:
                    ind_1 = index_roll(i+1,4)
                    ind_2 = index_roll(i-1,4)
                    if (all_nodes_dict[vcct_el_antena_nodes[ind_1]]['d'] > 
                        all_nodes_dict[vcct_el_antena_nodes[ind_2]]['d']):
                        ind_a_f = ind_1
                    else:
                        ind_a_f = ind_2
        else:
            rank_arr = np.zeros((4,3))
            for i,node in enumerate(vcct_el_antena_nodes):
            #builds a matrix with 'd' and 'delta_t'
                if node != 0:
                    node_coords = all_nodes_dict[node]['coords'][:2]
                    dl0 = np.linalg.norm(vcct_el_xy - node_coords)
                    rank_arr[i,0] = all_nodes_dict[node]['d']
                    rank_arr[i,1] = np.linalg.norm(all_nodes_dict[node]['delta'])/dl0
                    rank_arr[i,2] = i
            rank_arr = rank_arr[rank_arr[:,0].argsort()]
            rank_arr = rank_arr[rank_arr[:,1].argsort(kind='mergesort')]    
            ind_a_f =  int(rank_arr[-1,2])

    return ind_a_f


def along_nearest_mesh_line(vcct_el,ns2crack,a_node_prist_along_normal,vcct_el_dict,all_nodes_dict):    

    el_type = vcct_el_type(vcct_el,vcct_el_dict)
    if el_type == 'm':
        if len(a_node_prist_along_normal) == 0:
            if all_nodes_dict[vcct_el]['d'] > 0.001:
                a_node_next = all_nodes_dict[vcct_el]['a_node_next_p']
            else:
                a_nodes_indexes = [0,2,4,6]
                no_antena_nodes = 4
                a_node_pristine = np.zeros(no_antena_nodes)
                vcct_el_nodes = vcct_el_dict[vcct_el]
                #order the antena nodes by direction: closest to the normal considered
                for j,a_node in enumerate(vcct_el_nodes[a_nodes_indexes]):
                    if a_node != 0:
                        if all_nodes_dict[a_node]['d'] < 0.001:
                            a_node_next = a_node
                            break
        else:
            a_node_next = a_node_prist_along_normal[0]
        ns2crack = all_nodes_dict[a_node_next]['coords'][:2]  - all_nodes_dict[vcct_el]['coords'][:2]
        ns2crack = np.array([ns2crack/np.linalg.norm(ns2crack)])
        sample_xy = np.array([2.0*all_nodes_dict[vcct_el]['coords'][:2]-all_nodes_dict[a_node_next]['coords'][:2] ])  
        a_node_prist_along_normal = np.array([a_node_next])

    return sample_xy,ns2crack,a_node_prist_along_normal

def cfront_faces_asc(vcct_el,vcct_el_dict,all_nodes_dict,ind_a_f,option = 'carvalho'):

    a_nodes_indexes = [0,2,4,6]
    ind_side_a_nodes = np.array([a_nodes_indexes[index_roll(ind_a_f + 1,4)],
                              a_nodes_indexes[index_roll(ind_a_f - 1,4)]])
    
    cf_vs = np.empty((0,2))
    
    vcct_adj_els = vcct_el_dict[vcct_el]
    inc_ori = np.array([1,-1])
    vcct_el_xy = all_nodes_dict[vcct_el]['coords'][:2]   
    vcct_el_dstat = all_nodes_dict[vcct_el]['d']
    for i,ind in enumerate(ind_side_a_nodes):
        node_ref = vcct_adj_els[ind]
        if node_ref != 0:       
            d_ref = all_nodes_dict[node_ref]['d'] #damage state of the antenna node
            if d_ref < 0.0:
                d_ref = 0.0
                all_nodes_dict[node_ref]['d'] = 0.0
                print('warning: d registered negative(el,d) - reset to zero:',node_ref,d_ref)
            
            if d_ref == 0.00:    
                if option == 'carvalho':
                    node_d_0 = vcct_adj_els[index_roll(ind - inc_ori[i],8)]   
                elif option == 'carvalho_no_corner':
                    node_d_0 = node_ref
                node_d_1 = node_ref
            elif d_ref > 0.0:
                node_d_0 = node_ref
                node_d_1 = vcct_adj_els[index_roll(ind + inc_ori[i],8)]

            node_d_0_coords = all_nodes_dict[node_d_0]['coords'][:2]
            node_d_1_coords = all_nodes_dict[node_d_1]['coords'][:2]
            cf_v_tip = node_d_0_coords+(node_d_1_coords-node_d_0_coords)*all_nodes_dict[node_d_0]['d']
            
            vcct_el_next_node = vcct_adj_els[a_nodes_indexes[index_roll(ind_a_f + 2,4)]]
            if vcct_el_next_node != 0:
                vcct_el_next_node_xy = all_nodes_dict[vcct_el_next_node]['coords'][:2]
            else:
                vcct_el_next_node_xy = vcct_el_xy
            cf_vs = np.append(cf_vs,np.array([cf_v_tip - (vcct_el_xy+(vcct_el_next_node_xy-vcct_el_xy)*vcct_el_dstat)]),axis=0)
            
    #all_nodes_dict[vcct_el]['cf_v'] = cf_v

    return cf_vs


def VCCT_el_G_samples(all_els_dict,all_nodes_dict,node_2vcct_adj_els_dict,vcct_el_dict,analysis_dict):
    

    vcct_els_cfront = vcct_els_crackfront(vcct_el_dict,all_nodes_dict,option_node = 'dneq1',option_a_nodes='dgt0')

    #option_cf_v: 'boeing' or 'carvalho'
    #'boeing': ignores the crack shape (cf_v)
    #'carvalho': has the crack shape into account, calculates two vectors corresponding to the local crack faces
    #            only considers the nodes along the mesh line parallel to the crack growth
    #option_cf_v = 'carvalho' 
    option_cf_v = 'carvalho'
    
    #option_n2crack: 'boeing' or 'carvalho'
    #'boeing': normals are calculated along the direction of failed antena nodes to the node being assessed
    #'carvalho': calculates the normals having the crack shape into account - requires option_cf_v = 'carvalho'
    #'carvalho2': calculates a single normal corresponding to the bisection
    #considers three cases: 2 normals to each crack face, one along the bisection of the the normals to the crack faces
    option_n2crack = 'carvalho2'

    #option_prist_along_normal: 'boeing' or 'carvalho'
    #'boeing': finds pristine node along normal direction
    #'carvalho': finds pristine node along direction having into account the crack shape - requires  option_cf_v = option_n2crack = 'carvalho'
    option_prist_along_normal = 'boeing'

    #option_disp_sample_pos: 'boeing' or 'carvalho'
    #'boeing': computes sampling points based on the pristine nodes along the normals
    #'carvalho': calculates the sampling points based on cf_v and the normals to the crack front
    #           requires option_cf_v = option_n2crack = option_disp_sample_pos = 'carvalho'
    #option_disp_sample_pos = 'normal'
    option_disp_sample_pos = 'boeing-mod'
    #option_disp_sample_pos = 'normal_along_antena'
    
    #option_force_sample_points: 'next_a_node' or 'sample'
    #'next_a_node': uses the next antena node to compute the force, regardless of orientation, for a case in which ERR is computed explicitly at intermediate locations
    #'sample': samples force based on the displacement sampling point
    option_force_sample_points = 'next_a_node'


    #option_calc_area_sample_points: 'carvalho' or 'boeing'
    #'carvalho': accounts for the crack shape, requires:  option_cf_v = 'carvalho'
    #'boeing': ignores crack shape
    option_calc_area_sample_points = 'carvalho2'

    option_G_sample_points = 'nonlocal'
    vcct_els_cfront_dict = {}
    vcct_els_edge_of_model_dict ={}  
    
    #fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9, 9))
    #color = 'green'

    for i,vcct_el in enumerate(vcct_els_cfront):
        
         
        #coordinates of the vcct element (node)
        vcct_el_coords = all_nodes_dict[vcct_el]['coords']
        #print('vcct_el,vcct_el_coords',vcct_el,vcct_el_coords)
        #adjacent nodes to the vcct element (antena + corner)
        vcct_adj_els = vcct_el_dict[vcct_el]
        
        #regular cohesive elements adjacent to the node 'vcct_el'
        reg_adj_els = node_2vcct_adj_els_dict[vcct_el]
        
        #el_type....
        el_type = vcct_el_type(vcct_el,vcct_el_dict)
        
        #finds wake (failed) node index:
        a_node_wake_index = finds_reference_node_ind(vcct_el,vcct_el_coords,el_type,
                                           vcct_el_dict,all_nodes_dict)

        #for debuging purposes
        #ind_ref = finds_reference_node_ind_asc(vcct_el,vcct_el_dict,all_nodes_dict) 
      
        a_node_wake = vcct_el_dict[vcct_el][a_node_wake_index]

        
        #pristine node along the mesh line:
        a_node_prist_along_normal_ind = next_a_node_ind_by_inc(a_nodes_indexes,a_node_wake_index, 2)
        a_node_prist_along_normal = vcct_el_dict[vcct_el][a_node_prist_along_normal_ind]

        a_node_prist_along_normal_asc = np.array([a_node_prist_along_normal])
        #it seems the calculation is perfoormed even if the node at the wake is not failed
        #considering adding that condition to the if statement for a faster procedure
        

        if a_node_prist_along_normal == 0:
            print('a_node_prist_long_normal',a_node_prist_along_normal)
            #edge of a model, no node ahead of the node being evaluated along the reference node direction
            vcct_els_edge_of_model_dict[vcct_el] = {}
            #reference node number, a_node_ref
            
            vcct_els_edge_of_model_dict[vcct_el]['ref_node'] = a_node_wake
            
        else: #elif all_nodes_dict[a_node_wake]['d'] == 1.0:
            
            #determines the type of the VCCT element: corner:'c', edge:'e' or middle 'm'
            # if vcct_el == 491:
            #     print('vcct_el,491')
            cf_vs = cfront_faces(vcct_el,vcct_el_coords,vcct_adj_els,vcct_el_dict,
                                      all_nodes_dict,a_node_wake_index,a_node_prist_along_normal,option = 'corner_nodes')

        
            area = calc_area_vcct_el_cfront(vcct_el,vcct_el_coords,vcct_adj_els,a_node_wake_index,
                                            a_node_prist_along_normal,cf_vs,
                                         node_2vcct_adj_els_dict,all_nodes_dict,all_els_dict)
             
            #print('area',area)
            #print('vcct_el_coords',vcct_el_coords)
            
        
            all_nodes_dict[vcct_el]['cf_v'] = cf_vs
            ns2crack = calc_normals_2_crack(vcct_el,a_node_prist_along_normal,vcct_el_dict,all_nodes_dict,option_n2crack)


            #determines the sample point
            sp_coords_w,sp_coords_f = sample_points_coords(vcct_el,vcct_el_coords,vcct_adj_els,
                                                a_node_prist_along_normal,a_node_wake,a_node_wake_index,
                                                ns2crack,vcct_el_dict,node_2vcct_adj_els_dict,
                                                all_els_dict,all_nodes_dict,option_disp_sample_pos)

            
            if option_disp_sample_pos == 'boeing-mod':

                delta_sample_points = np.array([all_nodes_dict[a_node_wake]['delta']])
                
            else:
                #
                #calculate displacement jumps at sample points:
                delta_sample_points = eval_sample_point(sp_coords_w,'delta',reg_adj_els,
                                                     all_els_dict,all_nodes_dict)

            
            if option_disp_sample_pos == 'boeing-mod':

                force_sample_points = np.array([all_nodes_dict[a_node_prist_along_normal]['force']])
                
            else:
                  
            #calculate force at crack front projection
                force_sample_points = eval_sample_point(sp_coords_f,'force',reg_adj_els,
                                                     all_els_dict,all_nodes_dict)

           
             
             #calculate err
            #options: not defined yet - needs to be re-coded
            G_sample_points = calc_G_sample_points(vcct_el,delta_sample_points,force_sample_points,area,
                                 all_nodes_dict)


            #arb crack version:
            if asc != 'on':
                vcct_els_cfront_dict[vcct_el] = {'G_s': G_sample_points}
                vcct_els_cfront_dict[vcct_el]['dl0_s'] = dl0_sample_points(sp_coords_w,vcct_el,all_nodes_dict)
                vcct_els_cfront_dict[vcct_el]['delta_s'] = delta_sample_points
                vcct_els_cfront_dict[vcct_el]['a_node_next_p'] = np.array([a_node_prist_along_normal])
                vcct_els_cfront_dict[vcct_el]['a_node_ref'] = np.array([a_node_wake])#*len(a_node_prist_along_normal))

            #if all_nodes_dict[a_node_wake]['d'] == 1.0:  
            #    print('----------vcct el:',vcct_el)
            #    print('force_sample_points',force_sample_points)
            #    print('delta_sample_points',delta_sample_points)
            #print('G_sample_points',G_sample_points)
            #    print('dl0_s',vcct_els_cfront_dict[vcct_el]['dl0_s'])
            #    print('area',area)
            #    print('sp_coords_w',sp_coords_w)
            #    print('------------------------------')
            #    print('------------------------------')
                    

            
            if asc == 'on':
                cf_vs_asc = cfront_faces_asc(vcct_el,vcct_el_dict,all_nodes_dict,ind_ref,option_cf_v)

                all_nodes_dict[vcct_el]['cf_v'] = cf_vs_asc
                ns2crack_asc = calc_normals_2_crack_asc(vcct_el,[a_node_prist_along_normal],vcct_el_dict,all_nodes_dict,option_n2crack)

                a_node_ref = vcct_el_dict[vcct_el][a_nodes_indexes[ind_ref]]
                 
                sp_coords_w_asc = disp_sample_pos(vcct_el,
                                                 a_node_prist_along_normal_asc,a_node_ref,ind_ref,
                                                 ns2crack_asc,vcct_el_dict,node_2vcct_adj_els_dict,
                                                       all_els_dict,all_nodes_dict,'carvalho2')


                ind_inside_el,a_els_with_s_points = adj_els_with_sample_points(vcct_el,
                                                sp_coords_w_asc,node_2vcct_adj_els_dict,all_nodes_dict,all_els_dict)
            
                
                if len(sp_coords_w_asc[ind_inside_el]) == 0:#needs cleanup
                    sp_coords_w_asc,ns2crack,a_node_prist_along_normal = along_nearest_mesh_line(vcct_el,ns2crack,
                                                                             [a_node_prist_along_normal],
                                                                             vcct_el_dict,all_nodes_dict)
                    
                    a_node_prist_along_normal = a_node_prist_along_normal[0]
                    ind_inside_el,a_els_with_s_points = adj_els_with_sample_points(vcct_el,
                                                    sp_coords_w_asc,node_2vcct_adj_els_dict,all_nodes_dict,all_els_dict)
                    
                    
                #calculate displacement jumps at sample points, asc version:
                delta_sample_points_asc = calc_delta_sample_points_asc(sp_coords_w_asc,
                                                                    a_els_with_s_points,all_els_dict,all_nodes_dict)
             
             
                #calculate force at crack front projection, asc version
                force_sample_points_asc = calc_force_sample_points_asc(sp_coords_w_asc,a_node_prist_along_normal_asc,
                                                                    vcct_el,node_2vcct_adj_els_dict,all_els_dict,
                                                                    all_nodes_dict,option_force_sample_points)
             
                #area calculation, asc version
                area_asc = calc_area_sample_points_asc(vcct_el,ind_ref,sp_coords_w_asc,a_node_prist_along_normal_asc,vcct_el_dict,
                                                      all_nodes_dict,all_els_dict,option_calc_area_sample_points)  
                
                
                #options: not defined yet - needs to be re-coded
                G_sample_points_asc = calc_G_sample_points(vcct_el,delta_sample_points_asc,force_sample_points_asc,area_asc,
                                     all_nodes_dict)

         
                vcct_els_cfront_dict[vcct_el] = {'G_s': G_sample_points_asc}
                vcct_els_cfront_dict[vcct_el]['dl0_s'] = dl0_sample_points_asc(sp_coords_w_asc,vcct_el,all_nodes_dict)
                vcct_els_cfront_dict[vcct_el]['delta_s'] = delta_sample_points_asc
                vcct_els_cfront_dict[vcct_el]['a_node_next_p'] = [a_node_prist_along_normal]
                vcct_els_cfront_dict[vcct_el]['a_node_ref'] =  np.array([a_node_ref]*len(a_node_prist_along_normal_asc))
                 
                #if all_nodes_dict[a_node_wake]['d'] == 1.0:
                #     print('----------vcct el:',vcct_el)
                #     print('force_sample_points_asc',force_sample_points_asc)
                #     print('delta_sample_points_asc',delta_sample_points_asc)
                #     print('G_sample_points_asc',G_sample_points_asc)
                #     print('dl0_s_asc',vcct_els_cfront_dict[vcct_el]['dl0_s'])
                #     print('area_asc',area_asc)
                #     print('sp_coords_w_asc',sp_coords_w_asc)
                #     print('------------------------------')
                #     print('------------------------------')

                if abs(np.linalg.norm(G_sample_points) - np.linalg.norm(G_sample_points_asc)) > 0.1:
                    sys.exit("algorithms differing too much - check this case") 
            
            
    return vcct_els_cfront_dict,vcct_els_edge_of_model_dict


#global variables
a_nodes_indexes = [0,2,4,6]
a_nodes_indexes_w_corners = [0,1,2,3,4,5,6,7]
asc = 'off'
if __name__ == "__main__":
    
    current_dir = inp_vars.current_dir
    input_file = inp_vars.input_file
        
    verbose_plot = False
    verbose_vcct_el_inc = True
    test_init = False
    test_main = False
    
    if os.path.isfile(current_dir+input_file+'all_els_dict.pickle') == False:      
        VCCT_init()
    else:
        VCCT_main()
        #vcct_els_cfront_dict,vcct_el_dict,all_els_dict,all_nodes_dict,node_2vcct_adj_els_dict = VCCT_main()
       
        
        
        #for debuging:
       
        # vcct_el = 10608
        # a_nodes_indexes = [0,2,4,6]
        # vcct_el_antena_nodes = vcct_el_dict[vcct_el][a_nodes_indexes]    
        # #find reference antenna node
        # a_sel_arr = np.array([[all_nodes_dict[a_node]['d'],np.linalg.norm(all_nodes_dict[a_node]['delta'])] 
        #                       for a_node in vcct_el_antena_nodes])
                  
        #all_nodes_dict = reset_force_delta(vcct_els_cfront,all_nodes_dict,node_2vcct_adj_els_dict,all_els_dict)
        
        
        
