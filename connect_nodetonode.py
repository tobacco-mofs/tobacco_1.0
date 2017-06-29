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
import re 
import fnmatch 
import itertools 
from neighbors import neighbor_edges, neighbor_vertices
from nodes import __node_properties
from transformations import superimposition_matrix
from placing_bb_func import __place_bb
from operator import itemgetter



#Connect two nodes to each other, no edge.  


def __node_to_node(unit_cell, e1_coord, v12_nbors_for_e1, node_1_elements, node_2_elements, node_1_frac_coord, node_2_frac_coord, node_1_connectivity, node_2_connectivity):
	elements=[]
	frac_coords=[]
	#Update indexes of atoms in structure.  
	for j in range(len(node_1_elements)):
		if len(node_1_elements[j])==1:
			index = re.split('(\d+)', node_1_elements[j][0][0])
			index_number = int(index[1])+len(node_1_elements[0])*j
			index_number = str(index_number)
			updated_index = index[0]+index_number
			elements.append(index[0]+index_number + ' ' + node_1_elements[j][0][1])
			frac_coords.append(node_1_frac_coord[j])
		elif not len(node_1_elements[j])==1:
			for i in range(len(node_1_elements[j])):
				index = re.split('(\d+)', node_1_elements[j][i][0])
				index_number = int(index[1])+len(node_1_elements[0])*j
				index_number = str(index_number)
				updated_index = index[0]+index_number
				elements.append(index[0]+index_number + ' ' + node_1_elements[j][i][1])
				frac_coords.append(node_1_frac_coord[j][i])
			
	
	for j in range(len(node_2_elements)):
		if len(node_2_elements[j])==1:
			index = re.split('(\d+)', node_2_elements[j][0][0])
			index_number = int(index[1])+len(node_2_elements[0])*j + len(node_1_elements)*len(node_1_elements[0])
			index_number = str(index_number)
			updated_index = index[0]+index_number
			elements.append(index[0]+index_number + ' ' + node_2_elements[j][0][1])
			frac_coords.append(node_2_frac_coord[j])
		elif not len(node_2_elements[j])==1:
			for i in range(len(node_2_elements[j])):
				index = re.split('(\d+)', node_2_elements[j][i][0])
				index_number = int(index[1])+len(node_2_elements[0])*j + len(node_1_elements)*len(node_1_elements[0])
				index_number = str(index_number)
				updated_index = index[0]+index_number
				elements.append(index[0]+index_number + ' ' + node_2_elements[j][i][1])
				frac_coords.append(node_2_frac_coord[j][i])
			
			
	#Update indexes of atoms in connectivity
	connectivity=[]
	for j in range(len(node_1_connectivity)):
		for i in range(len(node_1_connectivity[j])):
			index_1 = re.split('(\d+)', node_1_connectivity[j][i][0])
			index_number_one = int(index_1[1])+len(node_1_elements[0])*j
			index_number_one = str(index_number_one)
			updated_index_one = index_1[0]+index_number_one
			index_2 = re.split('(\d+)', node_1_connectivity[j][i][1])
			index_number_two = int(index_2[1])+len(node_1_elements[0])*j
			index_number_two = str(index_number_two)
			updated_index_two = index_2[0]+index_number_two
			connectivity.append(updated_index_one + ' ' + updated_index_two + ' ' + node_1_connectivity[j][i][2] + ' ' + node_1_connectivity[j][i][3] + ' ' + node_1_connectivity[j][i][4])
	
	for j in range(len(node_2_connectivity)):
		for i in range(len(node_2_connectivity[j])):
			index_1 = re.split('(\d+)', node_2_connectivity[j][i][0])
			index_number_one = int(index_1[1])+len(node_2_elements[0])*j + len(node_1_elements)*len(node_1_elements[0])
			index_number_one = str(index_number_one)
			updated_index_one = index_1[0]+index_number_one
			index_2 = re.split('(\d+)', node_2_connectivity[j][i][1])
			index_number_two = int(index_2[1])+len(node_2_elements[0])*j + len(node_1_elements)*len(node_1_elements[0])
			index_number_two = str(index_number_two)
			updated_index_two = index_2[0]+index_number_two
			connectivity.append(updated_index_one + ' ' + updated_index_two + ' ' + node_2_connectivity[j][i][2] + ' ' + node_2_connectivity[j][i][3] + ' ' + node_2_connectivity[j][i][4])
	
	
	
	
	
	####Establishing a connection when there is no edge
	for k in range(len(v12_nbors_for_e1)):
		node_1_connection = node_1_elements[v12_nbors_for_e1[k][0]] 
		node_2_connection = node_2_elements[v12_nbors_for_e1[k][1]]
		
		
		
		
		
		location_1=[]
		connection_frac_coord_1=[]
		location_2=[]
		connection_frac_coord_2=[]
		for i in [i for i, x in enumerate(node_1_connection) if x[0][0][0] == 'X']:
			location_1.append(i)
			if len(node_1_connection)==1:
				connection_frac_coord_1.append(node_1_frac_coord[v12_nbors_for_e1[k][0]])
			elif not len(node_1_connection)==1:
				connection_frac_coord_1.append(node_1_frac_coord[v12_nbors_for_e1[k][0]][i]) # grab connection coord of first node
		for i in [i for i, x in enumerate(node_2_connection) if x[0][0][0] == 'X']:
			location_2.append(i)
			if len(node_2_connection)==1:
				connection_frac_coord_2.append(node_2_frac_coord[v12_nbors_for_e1[k][1]])
			elif not len(node_2_connection)==1:
				connection_frac_coord_2.append(node_2_frac_coord[v12_nbors_for_e1[k][1]][i]) # grab connection coord of second node


		
		
		###Calculate vectors from node to edge, and from node to connection points
		centroid_x_1=[]
		centroid_y_1=[]
		centroid_z_1=[]
		for j in range(len(connection_frac_coord_1)):
			centroid_x_1.append(connection_frac_coord_1[j][0])
			centroid_y_1.append(connection_frac_coord_1[j][1])
			centroid_z_1.append(connection_frac_coord_1[j][2])
		centroid_x_1 = sum(centroid_x_1)/len(connection_frac_coord_1)
		centroid_y_1 = sum(centroid_y_1)/len(connection_frac_coord_1)
		centroid_z_1 = sum(centroid_z_1)/len(connection_frac_coord_1)

		centroid_x_2=[]
		centroid_y_2=[]
		centroid_z_2=[]
		for j in range(len(connection_frac_coord_2)):
			centroid_x_2.append(connection_frac_coord_2[j][0])
			centroid_y_2.append(connection_frac_coord_2[j][1])
			centroid_z_2.append(connection_frac_coord_2[j][2])
		centroid_x_2 = sum(centroid_x_2)/len(connection_frac_coord_2)
		centroid_y_2 = sum(centroid_y_2)/len(connection_frac_coord_2)
		centroid_z_2 = sum(centroid_z_2)/len(connection_frac_coord_2)
		
		centroid_1 = [centroid_x_1, centroid_y_1, centroid_z_1]
		centroid_2 = [centroid_x_2, centroid_y_2, centroid_z_2]
		##Adjust edge_coord for PBCs, node_1
		
		edge_pos_1 = e1_coord[k]
		
		diffa = edge_pos_1[0]-centroid_1[0]
		diffb = edge_pos_1[1]-centroid_1[1]
		diffc = edge_pos_1[2]-centroid_1[2]
		
		### PERIODIC BOUNDARY CONDITIONS
		if diffa > 0.5:
			edge_pos_1[0] = edge_pos_1[0] - 1
		elif diffa < -0.5:
			edge_pos_1[0] = edge_pos_1[0] + 1
		
		if diffb > 0.5:
			edge_pos_1[1] = edge_pos_1[1] - 1
		elif diffb < -0.5:
			edge_pos_1[1] = edge_pos_1[1] + 1
		
		if diffc > 0.5:
			edge_pos_1[2] = edge_pos_1[2] - 1
		elif diffc < -0.5:
			edge_pos_1[2] = edge_pos_1[2] + 1
		
		
		##Vector node to edge, node_1
		node_1_to_edge = edge_pos_1-centroid_1
		
		##Adjust edge_coord for PBCs, node_1
		
		edge_pos_2 = e1_coord[k]
		
		diffa = edge_pos_2[0]-centroid_2[0]
		diffb = edge_pos_2[1]-centroid_2[1]
		diffc = edge_pos_2[2]-centroid_2[2]
		
		### PERIODIC BOUNDARY CONDITIONS
		if diffa > 0.5:
			edge_pos_2[0] = edge_pos_2[0] - 1
		elif diffa < -0.5:
			edge_pos_2[0] = edge_pos_2[0] + 1
		
		if diffb > 0.5:
			edge_pos_2[1] = edge_pos_2[1] - 1
		elif diffb < -0.5:
			edge_pos_2[1] = edge_pos_2[1] + 1
		
		if diffc > 0.5:
			edge_pos_2[2] = edge_pos_2[2] - 1
		elif diffc < -0.5:
			edge_pos_2[2] = edge_pos_2[2] + 1
		
		
		##Vector node to edge, node_2
		
		node_2_to_edge = edge_pos_2-centroid_2
		
		##Vectors node to connection_site, node_1
		connection_1_vector=[]
		for i in range(len(connection_frac_coord_1)):
			connection_1_vector.append(connection_frac_coord_1[i]-centroid_1)
		
		
		connection_2_vector=[]
		for i in range(len(connection_frac_coord_2)):
			connection_2_vector.append(connection_frac_coord_2[i]-centroid_2)
		
		
		
		##Calculate angle between node-connection site and node-edge vectors, the smallest angles correspond to the sites to be connected
		angle_1=[]
		if len(node_1_connection)==1:
			angle_1.append(0)
		elif not len(node_1_connection)==1:
			for i in range(len(connection_1_vector)):
				angle=np.arccos(np.dot(node_1_to_edge, connection_1_vector[i])/(np.linalg.norm(node_1_to_edge)*np.linalg.norm(connection_1_vector[i])))*180/np.pi
				if np.isnan(angle)==True:
					angle= np.arccos(round(np.dot(node_1_to_edge, connection_1_vector[i])/(np.linalg.norm(node_1_to_edge)*np.linalg.norm(connection_1_vector[i]))))*180/np.pi
				angle_1.append(angle)
		angle_2=[]
		if len(node_2_connection)==1:
			angle_2.append(0)
		elif not len(node_2_connection)==1:
			for i in range(len(connection_2_vector)):
				angle=np.arccos(np.dot(node_2_to_edge, connection_2_vector[i])/(np.linalg.norm(node_2_to_edge)*np.linalg.norm(connection_2_vector[i])))*180/np.pi
				if np.isnan(angle)==True:
					angle= np.arccos(round(np.dot(node_2_to_edge, connection_2_vector[i])/(np.linalg.norm(node_2_to_edge)*np.linalg.norm(connection_2_vector[i]))))*180/np.pi
				angle_2.append(angle)
		
		index_of_connection_1 = min(enumerate(angle_1), key=itemgetter(1))[0] #connection site of node_1
		index_of_connection_2 = min(enumerate(angle_2), key=itemgetter(1))[0] #connection site of node_2
		
		
		
		##Coordinates of connection sites to be connected
		
		if len(node_1_connection)==1:
			connection_node_1_coord = node_1_frac_coord[v12_nbors_for_e1[k][0]]
		elif not len(node_1_connection)==1:
			connection_node_1_coord = node_1_frac_coord[v12_nbors_for_e1[k][0]][location_1[index_of_connection_1]]
		if len(node_2_connection)==1:
			connection_node_2_coord = node_2_frac_coord[v12_nbors_for_e1[k][1]]
		elif not len(node_2_connection)==1:
			connection_node_2_coord = node_2_frac_coord[v12_nbors_for_e1[k][1]][location_2[index_of_connection_2]]
		diffa_conn = connection_node_2_coord[0]-connection_node_1_coord[0]
		diffb_conn = connection_node_2_coord[1]-connection_node_1_coord[1]
		diffc_conn = connection_node_2_coord[2]-connection_node_1_coord[2]
		symmetry = [5,5,5]
		symm_alt = [5,5,5]
		### PERIODIC BOUNDARY CONDITIONS
		if diffa_conn > 0.5:
			symmetry[0]=4
			symm_alt[0] =6
		elif diffa_conn < -0.5:
			symmetry[0]=6
			symm_alt[0] =4
		if diffb_conn > 0.5:
			symmetry[1]=4
			symm_alt[1] = 6
		elif diffb_conn < -0.5:
			symmetry[1]=6
			symm_alt[1] = 4
		if diffc_conn > 0.5:
			symmetry[2]=4
			symm_alt[2] = 6
		elif diffc_conn < -0.5:
			symmetry[2]=6
			symm_alt[2]=4
		
		connection_bond = np.linalg.norm(np.dot(np.transpose(unit_cell), connection_node_1_coord) - np.dot(np.transpose(unit_cell), connection_node_2_coord))
		
		
		
		index_conn_1 = location_1[index_of_connection_1] + 1
		
		atom_index_conn_1 = index_conn_1 + len(node_1_elements[0])*v12_nbors_for_e1[k][0]
		
		index_conn_2 = location_2[index_of_connection_2] + 1
		atom_index_conn_2 = index_conn_2 + len(node_2_elements[0])*v12_nbors_for_e1[k][1] + len(node_1_elements)*len(node_1_elements[0])
		
		
		
		
		##Add connection to connectivity information, taking into account periodic boundary conditions
		if all(x==5 for x in symmetry):
			connectivity.append(node_1_connection[location_1[index_of_connection_1]][0][0] + str(atom_index_conn_1) + ' ' + node_2_connection[location_2[index_of_connection_2]][0][0] + str(atom_index_conn_2) + ' ' + str(connection_bond) + ' ' + '.' + ' ' + 'S')
		elif not all(x==5 for x in symmetry):#For when connection crosses boundaries
			connectivity.append(node_1_connection[location_1[index_of_connection_1]][0][0] + str(atom_index_conn_1) + ' ' + node_2_connection[location_2[index_of_connection_2]][0][0] + str(atom_index_conn_2) + ' ' + str(connection_bond) + ' ' + '1_' + str(symmetry[0])+str(symmetry[1]) + str(symmetry[2]) + ' ' + 'S')
			connectivity.append(node_2_connection[location_2[index_of_connection_2]][0][0] + str(atom_index_conn_2) + ' ' + node_1_connection[location_1[index_of_connection_1]][0][0] + str(atom_index_conn_1) + ' ' + str(connection_bond) + ' ' + '1_' + str(symm_alt[0])+str(symm_alt[1]) + str(symm_alt[2]) + ' ' + 'S')

	

	#Change to .cif format
	new_frac_coords=[]
	for j in range(len(frac_coords)):
		for i in range(len(frac_coords[j])):
			new_frac_coords.append('{:f}'.format(float(frac_coords[j][i])))
	new_frac_coords = np.asarray(new_frac_coords)
	new_frac_coords = np.asarray(np.split(new_frac_coords, len(elements)))
	element_and_frac_coord = np.asarray(np.column_stack((elements, new_frac_coords)))
	def fmtcols(mylist, cols):
		lines = ("\t".join(mylist[i:i+cols]) for i in xrange(0,len(mylist),cols))
		return '\n'.join(lines)
	new_connectivity = fmtcols(connectivity, 1)
	return element_and_frac_coord, new_connectivity

