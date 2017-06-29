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
from operator import itemgetter

##Placing edge building block when topology contains two nodes.

np.set_printoptions(threshold = sys.maxint)

def __place_edge(arg, unit_cell, edge_coord, node_1_coord, node_2_coord, node_1_neighbor, node_2_neighbor):
	edge= __edge_properties(arg)
	
	connection_site_edge = edge[0] 
	distance_connection_site = edge[1]
	angle_connection_site_pair = edge[2]
	connectivity = edge[3]
	elements = edge[4]
	element_coord = edge[5]
	
	
	
	
	
	
	bb_elements=[]
	bb_frac_coordinates=[]
	bb_connectivity=[]
	for j in range(len(edge_coord)):
		node_1_vector =[]
		node_2_vector =[]
		
		##Get coordinates of neighboring nodes to find connection vectors
		for i in range(len(node_1_neighbor[j])):
			node_1_vector.append(node_1_coord[node_1_neighbor[j][i]])
		
		for i in range(len(node_2_neighbor[j])):
			node_2_vector.append(node_2_coord[node_2_neighbor[j][i]])
		
		node_1_vector=np.asarray(node_1_vector) # fract. coord. of neigboring nodes
		node_2_vector=np.asarray(node_2_vector)
		
		
		edge_node_1_vector=[]
		for i in range(len(node_1_vector)):
			diffa = node_1_vector[i][0]-edge_coord[j][0]  
			diffb = node_1_vector[i][1]-edge_coord[j][1]  
			diffc = node_1_vector[i][2]-edge_coord[j][2]  
			
			### PERIODIC BOUNDARY CONDITIONS
			if diffa > 0.5:
				node_1_vector[i][0] = node_1_vector[i][0] - 1
			elif diffa < -0.5:
				node_1_vector[i][0] = node_1_vector[i][0] + 1
			
			if diffb > 0.5:
				node_1_vector[i][1] = node_1_vector[i][1] - 1
			elif diffb < -0.5:
				node_1_vector[i][1] = node_1_vector[i][1] + 1
			
			if diffc > 0.5:
				node_1_vector[i][2] = node_1_vector[i][2] - 1
			elif diffc < -0.5:
				node_1_vector[i][2] = node_1_vector[i][2] + 1
			
			
			edge_node_1_vector.append(node_1_vector[i] - edge_coord[j])
		
		edge_node_1_vector =np.asarray(edge_node_1_vector) ## fract vectors edge to node adjusted for PBCs
		
		node_1_vector_real=[]
		for i in range(len(edge_node_1_vector)):
			vector_1_real = np.dot(np.transpose(unit_cell), edge_node_1_vector[i])
			node_1_vector_real.append(vector_1_real)
		
			
		
		node_1_vector_real = np.asarray(node_1_vector_real) # real (not fractional) coord vector of network
		
		edge_node_2_vector=[]
		for i in range(len(node_2_vector)):
			diffa = node_2_vector[i][0]-edge_coord[j][0]  
			diffb = node_2_vector[i][1]-edge_coord[j][1]  
			diffc = node_2_vector[i][2]-edge_coord[j][2]  
			
			### PERIODIC BOUNDARY CONDITIONS
			if diffa > 0.5:
				node_2_vector[i][0] = node_2_vector[i][0] - 1
			elif diffa < -0.5:
				node_2_vector[i][0] = node_2_vector[i][0] + 1
			
			if diffb > 0.5:
				node_2_vector[i][1] = node_2_vector[i][1] - 1
			elif diffb < -0.5:
				node_2_vector[i][1] = node_2_vector[i][1] + 1
			
			if diffc > 0.5:
				node_2_vector[i][2] = node_2_vector[i][2] - 1
			elif diffc < -0.5:
				node_2_vector[i][2] = node_2_vector[i][2] + 1
			
			
			edge_node_2_vector.append(node_2_vector[i] - edge_coord[j])
		
		edge_node_2_vector =np.asarray(edge_node_2_vector) ## fract vectors from edge to node 2 adjusted for PBCs
		
		node_2_vector_real=[]
		for i in range(len(edge_node_2_vector)):
			vector_2_real = np.dot(np.transpose(unit_cell), edge_node_2_vector[i])
			node_2_vector_real.append(vector_2_real)
		
		
		
		node_2_vector_real = np.asarray(node_2_vector_real) # real (not fractional) coord vector of network
		
		edge_coord_real = np.dot(np.transpose(unit_cell), edge_coord[j]) # real coord of network edge (centroid of edge)
		
		
		
		norm_node_1_vector_real=[]
		for i in range(len(node_1_vector_real)):
			norm = node_1_vector_real[i]/np.linalg.norm(node_1_vector_real[i])
			norm_node_1_vector_real.append(norm)
		
		norm_node_1_vector_real = np.asarray(norm_node_1_vector_real) # normalized network vectors
		
		norm_node_2_vector_real=[]
		for i in range(len(node_2_vector_real)):
			norm = node_2_vector_real[i]/np.linalg.norm(node_2_vector_real[i])
			norm_node_2_vector_real.append(norm)
		
		norm_node_2_vector_real = np.asarray(norm_node_2_vector_real) # normalized network vectors
	
		connection_edge=[]
		connection_1 = norm_node_1_vector_real[0]*distance_connection_site[0]
		connection_2 = norm_node_2_vector_real[0]*distance_connection_site[1] ##coordinates to where edge connection sites should be placed
		
		connection_edge.append(connection_1)
		connection_edge.append(connection_2)
		connection_edge.append(np.cross(connection_1, connection_2))
		connection_edge = np.asarray(connection_edge)

		connection_site=[]
		connection_site.append(connection_site_edge[0])
		connection_site.append(connection_site_edge[1])
		connection_site.append(np.cross(connection_site_edge[0], connection_site_edge[1]))
		connection_site = np.asarray(connection_site)
		
		
		
	
	
		
		
		##Calculate transformation matrix, using quaternions, to map building block vectors onto network vectors
		tfmatrix = superimposition_matrix(np.transpose(connection_site), np.transpose(connection_edge), usesvd=False)
		
		
		connection_site_plusone = np.append(connection_site, np.ones([len(connection_site),1]),1) # add a column of ones for dimension agreement
		
		tf_connection_site =[]
		for i in range(len(connection_site)):
			new_sites = np.dot(tfmatrix, connection_site_plusone[i]) #apply transformation matrix to each building block vector
			tf_connection_site.append(new_sites)
		
		tf_connection_site = np.asarray(tf_connection_site) #coordinates of building block connection sites, mapped onto network node sites
		tf_connection_site = tf_connection_site[:, :-1] #remove the column of ones, to obtain final set of coordinates
		
		###Apply transformation matrix to all atoms in building block
		element_coord_plusone = np.append(element_coord, np.ones([len(element_coord),1]),1)
		
		tf_element_coord=[]
		for i in range(len(element_coord)):
			new_element= np.dot(tfmatrix, element_coord_plusone[i])
			tf_element_coord.append(new_element)
		tf_element_coord = np.asarray(tf_element_coord)
		tf_element_coord = tf_element_coord[:, :-1]
		
		
		tf_frac_element_coord=[]
		for i in range(len(tf_element_coord)):
			frac_element_coord = np.dot(np.transpose(np.linalg.inv(unit_cell)), tf_element_coord[i])
			frac_element_coord = frac_element_coord + edge_coord[j]
			tf_frac_element_coord.append(frac_element_coord)
		
		tf_frac_element_coord = np.asarray(tf_frac_element_coord) #edge after transformation, in frac coords
		
		bb_elements.append(elements)
		bb_frac_coordinates.append(tf_frac_element_coord)
		bb_connectivity.append(connectivity)
	bb_elements = np.asarray(bb_elements)
	bb_frac_coordinates = np.asarray(bb_frac_coordinates)
	bb_connectivity=np.asarray(bb_connectivity)
	
	return bb_elements, bb_frac_coordinates, bb_connectivity


















