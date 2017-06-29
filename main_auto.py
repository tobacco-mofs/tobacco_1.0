#################
 #
 # This file is part of
 # ToBaCCo - Topologically-Based Crystal Constructor
 #
 # Copyright 2017 Yamil J. Colon <yamilcolon2015@u.northwestern.edu>
 #                Diego Gomez-Gualdron <dgomezgualdron@mines.edu>
 #                Ben Bucior <ben.bucior@gmail.com>
 #
 # ToBaCCo is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # ToBaCCo is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU Lesser General Public License for more details.
 #
 #################
import os
import numpy as np
import sys
import re 
import fnmatch 
import itertools 
from neighbors import neighbor_edges, neighbor_vertices
from nodes import __node_properties
from transformations import superimposition_matrix
from placing_bb_func import __place_bb
from connect_nodetonode import __node_to_node
from edges import __edge_properties
from placing_edge import __place_edge
from connect_edgeto2nodes import __edge_to_2nodes
from connect_edgetonode import __edge_to_node
from placing_edge_1node import __place_edge_1node
from write2cif import __write_cif
import contextlib

#makes a list of all template files in "templates" directory
template_list = os.listdir("templates")
np.set_printoptions(threshold=sys.maxint)
node_list = os.listdir("nodes_bb")
edge_list = os.listdir("edges_bb")
#iterate through each template file to extract relevant information
for template in template_list:
    template_now = open("templates/"+template, "r") 
    
    template_file = template_now.read()
    tf = template_file.split()
    
    #read topology
    topology = tf[1]
    #read unit cell parameters
    a_vector = np.asarray(map(float, [tf[3], tf[4], tf[5]]))
    b_vector = np.asarray(map(float, [tf[6], tf[7], tf[8]]))
    c_vector = np.asarray(map(float, [tf[9], tf[10],tf[11]]))
    unit_cell = np.row_stack([a_vector, b_vector, c_vector])
    
    #read  keys for vertices and edges 
    types_vertices = int(tf[12])


    if types_vertices == 1:
       n_vertex_one = int(tf[19])
       cn_vertex_one = int(tf[20])
       sym_vertex_one = int(tf[21])
       

      
       types_edges = int(tf[22])
       

       type_edge_one = int(tf[27])
       n_edge_one = int(tf[28])
       length_edge_one = 2 * float(tf[29]) 
          


       
    if types_vertices == 2:
       n_vertex_one = int(tf[19])
       cn_vertex_one = int(tf[20])
       sym_vertex_one = int(tf[21])

       n_vertex_two = int(tf[23])
       cn_vertex_two = int(tf[24])
       sym_vertex_two = int(tf[25])


       

       types_edges = int(tf[26])

       if types_edges == 1:
          type_edge_one = int(tf[31])
          n_edge_one = int(tf[32])
          length_edge_one = 2 * float(tf[33])



       if types_edges == 2:
          type_edge_one = int(tf[31])
          n_edge_one = int(tf[32])
          length_edge_one = 2 * float(tf[33])
          
          type_edge_two = int(tf[34])
          n_edge_two = int(tf[35])
          length_edge_two = 2 *  float(tf[36])

          


 
    #read coordinates for vertices and edges  
    
    #first case is regular nets (i.e. 1 type of vertex  and 1 type of edge)
    if types_vertices == 1  and types_edges == 1 :
       vertex_one_coord=[]
       #read coordinates for vertices
       for i in range(33, 33 + 3 * n_vertex_one) :
          vertex_one_coord.append(tf[i])
       v_one_coord =  np.asarray(map(float,(vertex_one_coord)))
       v1_coord = np.reshape(v_one_coord, (n_vertex_one, 3))

       checkpointindex = 36 + 3 * n_vertex_one
       edge_one_coord=[]
       #read coordinates for edges
       for i in range(checkpointindex, checkpointindex + 3 * n_edge_one): 
          edge_one_coord.append(tf[i])       

       e_one_coord = np.asarray(map(float,(edge_one_coord)))
       e1_coord = np.reshape(e_one_coord, (n_edge_one, 3))

       
    #second case is edge-transitive nets(i.e. 2 types of vertices and 1 type of edges
    if types_vertices == 2 and types_edges == 1 :
       vertex_one_coord=[]
       #read coordinates for first kind of vertices
       for i in range(37, 37 + 3* n_vertex_one) :
          vertex_one_coord.append(tf[i])
       v_one_coord =  np.asarray(map(float,(vertex_one_coord)))
       v1_coord = np.reshape(v_one_coord, (n_vertex_one, 3))

       checkpointindex = 39 + 3 * n_vertex_one
       vertex_two_coord=[]
       #read coordinates for second kind of vertices
       for i in range(checkpointindex, checkpointindex + 3 * n_vertex_two):
          vertex_two_coord.append(tf[i])
       v_two_coord = np.asarray(map(float, (vertex_two_coord)))
       v2_coord = np.reshape(v_two_coord, (n_vertex_two, 3))

       checkpointindex = checkpointindex + 3 * n_vertex_two + 3
       edge_one_coord=[]
       #read coordinates for edges
       for i in range(checkpointindex, checkpointindex + 3 * n_edge_one) :
          edge_one_coord.append(tf[i])
       e_one_coord = np.asarray(map(float,(edge_one_coord)))
       e1_coord = np.reshape(e_one_coord, (n_edge_one, 3))

    #third case is for nets with 2 types of vertices and 2 types of edges
    if types_vertices == 2 and types_edges == 2 :
       vertex_one_coord=[]
       #read coordinates for first kind of vertices
       for i in range(40, 40 + 3* n_vertex_one) :
          vertex_one_coord.append(tf[i])
       v_one_coord =  np.asarray(map(float,(vertex_one_coord)))
       v1_coord = np.reshape(v_one_coord, (n_vertex_one, 3))

       checkpointindex = 42 + 3 * n_vertex_one
       vertex_two_coord=[]
       #read coordinates for second kind of vertices
       for i in range(checkpointindex, checkpointindex + 3 * n_vertex_two):
          vertex_two_coord.append(tf[i])
       v_two_coord = np.asarray(map(float, (vertex_two_coord)))
       v2_coord = np.reshape(v_two_coord, (n_vertex_two, 3))

       checkpointindex = checkpointindex + 3 * n_vertex_two + 3
       edge_one_coord=[]
       #read coordinates for first kind of edges
       for i in range(checkpointindex, checkpointindex + 3 * n_edge_one) :
          edge_one_coord.append(tf[i])
       e_one_coord = np.asarray(map(float,(edge_one_coord)))
       e1_coord = np.reshape(e_one_coord, (n_edge_one, 3))

       checkpointindex = checkpointindex + 3 * n_edge_one + 2
       edge_two_coord=[]
       #read coordinates for second kind of edges
       for i in range(checkpointindex, checkpointindex + 3 * n_edge_two) :
           edge_two_coord.append(tf[i])
       e_two_coord = np.asarray(map(float, (edge_two_coord)))
       e2_coord = np.reshape(e_two_coord, (n_edge_two, 3))

       
    #create list of neighbors for regular net cases
    if types_vertices == 1 and types_edges == 1 :   
       e1_nbors_for_v1 = neighbor_edges(unit_cell, v1_coord, n_vertex_one, e1_coord, n_edge_one, cn_vertex_one,(0.5 * length_edge_one)) 

       v1_nbors_for_e1 = neighbor_vertices(unit_cell, v1_coord, n_vertex_one, e1_coord, n_edge_one, (0.5 * length_edge_one))
       v11_nbors_for_e1 = np.reshape(v1_nbors_for_e1, (n_edge_one, 2))

       reg_node=[]
       for i in range(len(node_list)):
           if fnmatch.fnmatch(node_list[i], 'sym_'+str(sym_vertex_one)+'_'+'*'):
              reg_node.append(node_list[i])
       reg_node=np.asarray(reg_node)#Identify appropriate nodes for this particular topology
       for i in range(len(reg_node)):
           for j in range(len(edge_list)):
               node_properties = __node_properties(reg_node[i])
               distance_node_reg = node_properties[1]
               edge_properties = __edge_properties(edge_list[j])
               distance_edge = edge_properties[1]
               scale = (2*distance_node_reg[0] + 2*distance_edge[0] + 3.0)/length_edge_one
               unit_cell = scale*unit_cell

               
               node_reg = __place_bb(reg_node[i], unit_cell, e1_coord, v1_coord, e1_nbors_for_v1)

               edge = __place_edge_1node(edge_list[j], unit_cell, e1_coord, v1_coord, v11_nbors_for_e1)
               node_1_elements = node_reg[0]
               node_1_frac_coord = node_reg[1]
               node_1_connectivity = node_reg[2] 
               edge_elements = edge[0]
               edge_frac_coord = edge[1]
               edge_connectivity = edge[2]
               connection = __edge_to_node(unit_cell, e1_coord, v11_nbors_for_e1, node_1_elements, edge_elements, node_1_frac_coord, edge_frac_coord, node_1_connectivity, edge_connectivity)
               elements_and_frac_coord = connection[0]
               connectivity = connection[1]
               elements_and_frac_coord =str(elements_and_frac_coord).replace('[','').replace("'", '').replace(']','')
               node_2='_'
               __write_cif(template, reg_node[i], node_2, edge_list[j], unit_cell, elements_and_frac_coord, connectivity)  
               unit_cell = unit_cell/scale

    if types_vertices == 2 and types_edges == 1 :
       e1_nbors_for_v1 = neighbor_edges(unit_cell, v1_coord, n_vertex_one, e1_coord, n_edge_one, cn_vertex_one,(0.5 * length_edge_one))

       e1_nbors_for_v2 = neighbor_edges(unit_cell, v2_coord, n_vertex_two, e1_coord, n_edge_one, cn_vertex_two,(0.5 * length_edge_one))

       v1_nbors_for_e1 = neighbor_vertices(unit_cell, v1_coord, n_vertex_one, e1_coord, n_edge_one, (0.5 * length_edge_one))
       v2_nbors_for_e1 = neighbor_vertices(unit_cell, v2_coord, n_vertex_two, e1_coord, n_edge_one, (0.5 * length_edge_one))
       v12_nbors_for_e1 = np.column_stack((v1_nbors_for_e1, v2_nbors_for_e1))

       node_1_top=[]
       node_2_top=[]   
       for i in range(len(node_list)):
           if fnmatch.fnmatch(node_list[i], 'sym_'+str(sym_vertex_one)+'_'+'*'):
              node_1_top.append(node_list[i])
           if fnmatch.fnmatch(node_list[i], 'sym_'+str(sym_vertex_two)+'_'+'*'):
              node_2_top.append(node_list[i])
       node_1_top=np.asarray(node_1_top)
       node_2_top=np.asarray(node_2_top)   
       for i in range(len(node_1_top)):
           for j in range(len(node_2_top)):
               node_1_properties = __node_properties(node_1_top[i])
               node_2_properties = __node_properties(node_2_top[j])
               distance_node_1 = node_1_properties[1]
               distance_node_2 = node_2_properties[1]
               scale = (distance_node_1[0]+ distance_node_2[0] + 1.5)/length_edge_one
               unit_cell = unit_cell*scale 

               
           ##Adjusting edge coordinates when having two different nodes
               for k in range(len(v12_nbors_for_e1)):
                   node_1_coord = v1_coord[v12_nbors_for_e1[k][0]]
                   node_2_coord = v2_coord[v12_nbors_for_e1[k][1]]
                   edge_coord = e1_coord[k]
                   
                   ##Adjust for PBCs, centered on edge and moving node
                   
                   ##Node 1
                   diffa_1= node_1_coord[0] - edge_coord[0]
                   diffb_1= node_1_coord[1] - edge_coord[1]
                   diffc_1= node_1_coord[2] - edge_coord[2]
                   
                   if diffa_1 > 0.5:
                    node_1_coord[0] = node_1_coord[0] - 1
                   elif diffa_1 < -0.5:
                    node_1_coord[0] = node_1_coord[0] + 1
                   
                   if diffb_1 > 0.5:
                    node_1_coord[1] = node_1_coord[1] - 1
                   elif diffb_1 < -0.5:
                    node_1_coord[1] = node_1_coord[1] + 1
                    
                   if diffc_1 > 0.5:
                    node_1_coord[2] = node_1_coord[2] - 1
                   elif diffc_1 < -0.5:
                    node_1_coord[2] = node_1_coord[2] + 1
                   
                   ##Node 2
                   diffa_2= node_2_coord[0] - edge_coord[0]
                   diffb_2= node_2_coord[1] - edge_coord[1]
                   diffc_2= node_2_coord[2] - edge_coord[2]
                   
                   if diffa_2 > 0.5:
                    node_2_coord[0] = node_2_coord[0] - 1
                   elif diffa_2 < -0.5:
                    node_2_coord[0] = node_2_coord[0] + 1
                   
                   if diffb_2 > 0.5:
                    node_2_coord[1] = node_2_coord[1] - 1
                   elif diffb_2 < -0.5:
                    node_2_coord[1] = node_2_coord[1] + 1
                   
                   if diffc_2 > 0.5:
                    node_2_coord[2] = node_2_coord[2] - 1
                   elif diffc_2 < -0.5:
                    node_2_coord[2] = node_2_coord[2] + 1
                   
                   #Change from fractional to 'real' coordinates
                   node_1_coord_real = np.dot(np.transpose(unit_cell), node_1_coord)
                   node_2_coord_real = np.dot(np.transpose(unit_cell), node_2_coord)
                   edge_coord_real = np.dot(np.transpose(unit_cell), edge_coord)
                   
                   
                   ##Calculate vectors from node to edge , normalize, and multiply by connection distance
                   node_1_to_edge = (edge_coord_real - node_1_coord_real)/np.linalg.norm(edge_coord_real - node_1_coord_real)*distance_node_1[0] + node_1_coord_real
                   node_2_to_edge = (edge_coord_real - node_2_coord_real)/np.linalg.norm(edge_coord_real - node_2_coord_real)*distance_node_2[0] + node_2_coord_real
                   
                   e1_coord[k][0] = (node_1_to_edge[0] + node_2_to_edge[0])/2
                   e1_coord[k][1] = (node_1_to_edge[1] + node_2_to_edge[1])/2
                   e1_coord[k][2] = (node_1_to_edge[2] + node_2_to_edge[2])/2
                   e1_coord[k] = np.dot(np.transpose(np.linalg.inv(unit_cell)), e1_coord[k])   
   
               node_1 = __place_bb(node_1_top[i], unit_cell, e1_coord, v1_coord, e1_nbors_for_v1)
               node_2 = __place_bb(node_2_top[j], unit_cell, e1_coord, v2_coord, e1_nbors_for_v2) 
               node_1_elements = node_1[0]
               node_1_frac_coord = node_1[1]
               node_1_connectivity = node_1[2] 
               node_2_elements = node_2[0]
               node_2_frac_coord = node_2[1]
               node_2_connectivity = node_2[2] 
               connection = __node_to_node(unit_cell, e1_coord, v12_nbors_for_e1, node_1_elements, node_2_elements, node_1_frac_coord, node_2_frac_coord, node_1_connectivity, node_2_connectivity)
               elements_and_frac_coord = connection[0]
               connectivity = connection[1]

               elements_and_frac_coord =str(elements_and_frac_coord).replace('[','').replace("'", '').replace(']','')
               edge = '_.cif'
               __write_cif(template, node_1_top[i], node_2_top[j], edge, unit_cell, elements_and_frac_coord, connectivity)  
               unit_cell = unit_cell/scale
               for m in range(len(edge_list)):
                          node_1_properties = __node_properties(node_1_top[i])
                          node_2_properties = __node_properties(node_2_top[j])
                          edge_properties = __edge_properties(edge_list[m])
                          distance_node_1 = node_1_properties[1]
                          distance_node_2 = node_2_properties[1]

                            
                            
                          distance_edge = edge_properties[1]
                          scale = (distance_node_1[0] + 2*distance_edge[0] + distance_node_2[0] + 3.0)/length_edge_one
                          unit_cell = scale*unit_cell           
                          ##Adjusting edge coordinates when having two different nodes
                          for k in range(len(v12_nbors_for_e1)):
                              node_1_coord = v1_coord[v12_nbors_for_e1[k][0]]
                              node_2_coord = v2_coord[v12_nbors_for_e1[k][1]]
                              edge_coord = e1_coord[k]
                              
                              ##Adjust for PBCs, centered on edge and moving node
                              
                              ##Node 1
                              diffa_1= node_1_coord[0] - edge_coord[0]
                              diffb_1= node_1_coord[1] - edge_coord[1]
                              diffc_1= node_1_coord[2] - edge_coord[2]
                              
                              if diffa_1 > 0.5:
                               node_1_coord[0] = node_1_coord[0] - 1
                              elif diffa_1 < -0.5:
                               node_1_coord[0] = node_1_coord[0] + 1
                              
                              if diffb_1 > 0.5:
                               node_1_coord[1] = node_1_coord[1] - 1
                              elif diffb_1 < -0.5:
                               node_1_coord[1] = node_1_coord[1] + 1
                               
                              if diffc_1 > 0.5:
                               node_1_coord[2] = node_1_coord[2] - 1
                              elif diffc_1 < -0.5:
                               node_1_coord[2] = node_1_coord[2] + 1
                              
                              ##Node 2
                              diffa_2= node_2_coord[0] - edge_coord[0]
                              diffb_2= node_2_coord[1] - edge_coord[1]
                              diffc_2= node_2_coord[2] - edge_coord[2]
                              
                              if diffa_2 > 0.5:
                               node_2_coord[0] = node_2_coord[0] - 1
                              elif diffa_2 < -0.5:
                               node_2_coord[0] = node_2_coord[0] + 1
                              
                              if diffb_2 > 0.5:
                               node_2_coord[1] = node_2_coord[1] - 1
                              elif diffb_2 < -0.5:
                               node_2_coord[1] = node_2_coord[1] + 1
                              
                              if diffc_2 > 0.5:
                               node_2_coord[2] = node_2_coord[2] - 1
                              elif diffc_2 < -0.5:
                               node_2_coord[2] = node_2_coord[2] + 1
                              
                              #Change from fractional to 'real' coordinates
                              node_1_coord_real = np.dot(np.transpose(unit_cell), node_1_coord)
                              node_2_coord_real = np.dot(np.transpose(unit_cell), node_2_coord)
                              edge_coord_real = np.dot(np.transpose(unit_cell), edge_coord)
                              
                              
                              ##Calculate vectors from node to edge , normalize, and multiply by connection distance
                              node_1_to_edge = (edge_coord_real - node_1_coord_real)/np.linalg.norm(edge_coord_real - node_1_coord_real)*distance_node_1[0] + node_1_coord_real
                              node_2_to_edge = (edge_coord_real - node_2_coord_real)/np.linalg.norm(edge_coord_real - node_2_coord_real)*distance_node_2[0] + node_2_coord_real
                              
                              e1_coord[k][0] = (node_1_to_edge[0] + node_2_to_edge[0])/2
                              e1_coord[k][1] = (node_1_to_edge[1] + node_2_to_edge[1])/2
                              e1_coord[k][2] = (node_1_to_edge[2] + node_2_to_edge[2])/2
                              e1_coord[k] = np.dot(np.transpose(np.linalg.inv(unit_cell)), e1_coord[k])
                          node_1 = __place_bb(node_1_top[i], unit_cell, e1_coord, v1_coord, e1_nbors_for_v1)
                          node_2 = __place_bb(node_2_top[j], unit_cell, e1_coord, v2_coord, e1_nbors_for_v2) 
                          edge = __place_edge(edge_list[m], unit_cell, e1_coord, v1_coord, v2_coord, v1_nbors_for_e1, v2_nbors_for_e1)
                          node_1_elements = node_1[0]
                          node_1_frac_coord = node_1[1]
                          node_1_connectivity = node_1[2] 
                          node_2_elements = node_2[0]
                          node_2_frac_coord = node_2[1]
                          node_2_connectivity = node_2[2] 
                          edge_elements = edge[0]
                          edge_frac_coord = edge[1]
                          edge_connectivity = edge[2]
                          connection = __edge_to_2nodes(unit_cell, e1_coord, v12_nbors_for_e1, node_1_elements, node_2_elements, edge_elements, node_1_frac_coord, node_2_frac_coord, edge_frac_coord, node_1_connectivity, node_2_connectivity, edge_connectivity)
                          elements_and_frac_coord = connection[0]
                          connectivity = connection[1]

                          elements_and_frac_coord =str(elements_and_frac_coord).replace('[','').replace("'", '').replace(']','')
                          __write_cif(template, node_1_top[i], node_2_top[j], edge_list[m], unit_cell, elements_and_frac_coord, connectivity)  
                          unit_cell = unit_cell/scale
                      
                      
                      
                      
    template_now.close    


