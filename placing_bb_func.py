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
from operator import itemgetter
import sys
import contextlib
import collections


#Placing node.


def __place_bb(arg, unit_cell, edge_coord, vertex_coord, edge_neighbor_vertex):
	node = __node_properties(arg)
	np.set_printoptions(threshold=sys.maxint)
	connect_site = node[0]
	distance_connection_site = node[1]
	angle_connection_site_pair = node[2]
	connectivity = node[3]
	elements = node[4]
	element_coord = node[5]
	
	
	
	
	
	
	bb_elements=[]
	bb_frac_coordinates=[]
	bb_connectivity=[]
	for j in range(len(vertex_coord)):
		
		if len(elements)==1: ### Special case when node building block is a single atom
			bb_elements.append(elements)
			bb_frac_coordinates.append(vertex_coord[j])
			continue
		edge_vector =[]
		
		##Get coordinates of neighboring edges to find connection vectors
		for i in range(len(edge_neighbor_vertex[j])):
			edge_vector.append(edge_coord[edge_neighbor_vertex[j][i]])
		
		edge_vector=np.asarray(edge_vector) # fract. coord. of neigboring edges
		
		
		
		node_vector=[]
		for i in range(len(edge_vector)):
			diffa = edge_vector[i][0]-vertex_coord[j][0]  
			diffb = edge_vector[i][1]-vertex_coord[j][1]  
			diffc = edge_vector[i][2]-vertex_coord[j][2]  
			
			### PERIODIC BOUNDARY CONDITIONS
			if diffa > 0.5:
				edge_vector[i][0] = edge_vector[i][0] - 1
			elif diffa < -0.5:
				edge_vector[i][0] = edge_vector[i][0] + 1
			
			if diffb > 0.5:
				edge_vector[i][1] = edge_vector[i][1] - 1
			elif diffb < -0.5:
				edge_vector[i][1] = edge_vector[i][1] + 1
			
			if diffc > 0.5:
				edge_vector[i][2] = edge_vector[i][2] - 1
			elif diffc < -0.5:
				edge_vector[i][2] = edge_vector[i][2] + 1
			
			
			node_vector.append(edge_vector[i] - vertex_coord[j])
		
		node_vector =np.asarray(node_vector) ## fract vectors from node to edge adjusted for PBCs
		
		node_vector_real=[]
		for i in range(len(node_vector)):
			vector_real = np.dot(np.transpose(unit_cell), node_vector[i])
			node_vector_real.append(vector_real)
		
			
		
		node_vector_real = np.asarray(node_vector_real) # real (not fractional) coord vector of network
		
		node_coord_real = np.dot(np.transpose(unit_cell), vertex_coord[j]) # real coord of network node (centroid of node)
		norm_node_vector_real=[]
		for i in range(len(node_vector_real)):
			norm = node_vector_real[i]/np.linalg.norm(node_vector_real[i])
			norm_node_vector_real.append(norm)
		
		norm_node_vector_real = np.asarray(norm_node_vector_real) # normalized network vectors
		
		
		connect_node=[]
		connection_node=[]
		for i in range(len(norm_node_vector_real)):
			connect = norm_node_vector_real[i]*distance_connection_site[i]
			connect_node.append(connect)
			connection_node.append(connect)
		
		
		
		connection_node=np.asarray(connection_node) ## coordinates to where node connection sites should be placed
		
		connection_site = []
		for i in range(len(connect_site)):
			connection_site.append(connect_site[i])
		connection_site = np.asarray(connection_site)
		
		
		### To deal with nodes with ONLY two connections.
		if len(connection_site)==2:
			bi_connection_site=[]
			bi_connection_node=[]
			#test_vector=[0, 0, 0]
			for i in range(len(connection_site)):
				bi_connection_site.append(connection_site[i])
				bi_connection_node.append(connection_node[i])
			#bi_connection_site.append(-connection_site[0])
			#bi_connection_site.append(-connection_site[1])
			#bi_connection_node.append(-connection_node[0])
			#bi_connection_node.append(-connection_node[1])
			bi_connection_site.append(np.cross(connection_site[0], connection_site[1]))
			bi_connection_site.append(np.cross(connection_site[1], connection_site[0]))
			bi_connection_node.append(np.cross(connection_node[1], connection_node[0]))
			bi_connection_node.append(np.cross(connection_node[0], connection_node[1]))
			#bi_connection_site.append(test_vector)
			#bi_connection_node.append(test_vector)
			connection_site=np.asarray(bi_connection_site)
			connection_node=np.asarray(bi_connection_node)
			#print "again", connection_site, len(connection_site)
			#print connection_node
		
		### To deal with *bct* topologies
		if len(connection_site)==10:
			angle_site_sum=[]
			angle_node_sum=[]
			distance_site_sum=[]
			distance_node_sum=[]
			for i in range(len(connection_site)):
				angle_site=[]
				angle_node=[]
				distance_site=[]
				distance_node=[]
				for k in range(len(connection_site)):
					angle_s=np.arccos(np.dot(connection_site[i], connection_site[k])/(np.linalg.norm(connection_site[i])*np.linalg.norm(connection_site[k])))*180/np.pi
					angle_n=np.arccos(np.dot(connection_node[i], connection_node[k])/(np.linalg.norm(connection_node[i])*np.linalg.norm(connection_node[k])))*180/np.pi
					dist_s = np.linalg.norm(connection_site[i] - connection_site[k])
					dist_n = np.linalg.norm(connection_node[i] - connection_node[k])
					
					if np.isnan(angle_s)==True:
						angle_s=np.arccos(round(np.dot(connection_site[i], connection_site[k])/(np.linalg.norm(connection_site[i])*np.linalg.norm(connection_site[k]))))*180/np.pi
					if np.isnan(angle_n)==True:
						angle_n=np.arccos(round(np.dot(connection_node[i], connection_node[k])/(np.linalg.norm(connection_node[i])*np.linalg.norm(connection_node[k]))))*180/np.pi
					angle_site.append(angle_s)
					angle_node.append(angle_n)
					distance_site.append(dist_s)
					distance_node.append(dist_n)
					counter_site = collections.Counter(np.around(distance_site,1))
					counter_node = collections.Counter(np.around(distance_node,1))
				angle_site_sum.append(sum(angle_site))
				angle_node_sum.append(sum(angle_node))
				distance_site_sum.append(sum(distance_site))
				distance_node_sum.append(sum(distance_node))
			location_dist_site=[]
			location_dist_node=[]
			index_dist_site = min(enumerate(distance_site_sum), key=itemgetter(1))[0]
			index_dist_node = min(enumerate(distance_node_sum), key=itemgetter(1))[0]
			location_dist_site.append(index_dist_site)
			location_dist_node.append(index_dist_node)
			distance_site_sum[index_dist_site]=1000
			distance_node_sum[index_dist_node]=1000
			index_dist_site = min(enumerate(distance_site_sum), key=itemgetter(1))[0]
			index_dist_node = min(enumerate(distance_node_sum), key=itemgetter(1))[0]
			location_dist_site.append(index_dist_site)
			location_dist_node.append(index_dist_node)
			location_dist_site = np.sort(location_dist_site)
			location_dist_node = np.sort(location_dist_node)
			
			
			index_site = max(enumerate(angle_site_sum), key=itemgetter(1))[0]
			angle_site_sum[index_site]=0
			index_site_1 = max(enumerate(angle_site_sum), key=itemgetter(1))[0]
			index_node = max(enumerate(angle_node_sum), key=itemgetter(1))[0]
			angle_node_sum[index_node]=0
			index_node_1 = max(enumerate(angle_node_sum), key=itemgetter(1))[0]
			dist_site =[]
			dist_node=[]
			
			location_add_site=[]
			location_add_node=[]
			for i in range(len(connection_site)):
				dist_site.append(np.linalg.norm(connection_site[location_dist_site[0]] - connection_site[i]))
				dist_node.append(np.linalg.norm(connection_node[location_dist_node[0]] - connection_node[i]))
			counter_site = collections.Counter(np.around(dist_site,0))
			counter_node = collections.Counter(np.around(dist_node,0))
			
			
			site_criterion = counter_site.most_common(1)[0][0]
			
			if site_criterion == counter_node.most_common(1)[0][0]:
				node_criterion = counter_node.most_common(1)[0][0]
			else:
				node_criterion = counter_node.most_common(2)[1][0]
			
			for i in range(len(dist_site)):
				if np.around(dist_site[i],0) == site_criterion:
					location_add_site.append(i)
				if np.around(dist_node[i],0) == node_criterion:
					location_add_node.append(i)

			connection_site = [connection_site[location_dist_site[1]], connection_site[location_dist_site[0]], connection_site[location_add_site[0]], connection_site[location_add_site[1]], connection_site[location_add_site[2]], connection_site[location_add_site[3]]]
			connection_node = [connection_node[location_dist_node[1]], connection_node[location_dist_node[0]], connection_node[location_add_node[0]], connection_node[location_add_node[1]], connection_node[location_add_node[2]], connection_node[location_add_node[3]]]


		
		
		
		## This part of the code orders vectors and then takes the ratio to find opposite vectors and, if they exist, perpendiculars.
		## This is to deal with topologies with have nodes with more than 8 connection sites.
		list_a = []
		list_a_1 =[]
		for i in range(len(connection_site)):
			list_a.append(np.dot(connection_site[0], connection_site[i]))
			list_a_1.append(np.dot(connection_site[1], connection_site[i]))
		
		list_a = np.asarray(list_a)
		list_a_1 = np.asarray(list_a_1)
		list_b = []
		list_b_1 =[]
		for i in range(len(connection_node)):
			list_b.append(np.dot(connection_node[0], connection_node[i]))
			list_b_1.append(np.dot(connection_node[1], connection_node[i]))
			
		list_b = np.asarray(list_b)
		list_b_1 = np.asarray(list_b_1)
		
		
		sigma = np.sort(list_a)
		sigma_1 = np.sort(list_a_1)
		tau = np.sort(list_b)
		tau_1 = np.sort(list_b_1)

		sorted_a = np.argsort(list_a)
		sorted_a_1 = np.argsort(list_a_1)
		sorted_b = np.argsort(list_b)
		sorted_b_1 = np.argsort(list_b_1)
		
		
		inner_distance_site = []
		for i in range(len(connection_site)):
			inner_distance_site.append(np.linalg.norm(connection_site[i]-connection_site[0]))
			
			
		
		
		sorted_sites = []
		sorted_sites_1 =[]
		sorted_nodes_1 =[]
		for i in range(len(connection_site)):
			sorted_sites.append(connection_site[i])
			sorted_sites_1.append(connection_site[i])
			sorted_nodes_1.append(connection_node[i])
			
		sorted_sites = np.asarray(sorted_sites)
		sorted_sites_1 = np.asarray(sorted_sites_1)
		sorted_nodes_1 = np.asarray(sorted_nodes_1)
		
		for i in range(len(sorted_sites)):
			sorted_sites[sorted_b[i]]=connection_site[sorted_a[i]]
			sorted_sites_1[sorted_b_1[i]] = connection_site[sorted_a_1[i]]
			sorted_nodes_1[sorted_a_1[i]] = connection_node[sorted_b_1[i]]
		
		connection_site = sorted_sites
		
		inner_dot_site_sorted=[]
		inner_dot_site_sorted_1 = []
		inner_dot_node_sorted=[]
		inner_dot_node_sorted_1=[]
		
		for i in range(len(sorted_sites)):
			inner_dot_site_sorted.append(np.dot(sorted_sites[0], sorted_sites[i]))
			inner_dot_node_sorted.append(np.dot(connection_node[0], connection_node[i]))
			inner_dot_site_sorted_1.append(np.dot(sorted_sites_1[1], sorted_sites_1[i]))
			inner_dot_node_sorted_1.append(np.dot(connection_node[1], connection_node[i]))

		ratio_sorted = np.divide(inner_dot_site_sorted, inner_dot_node_sorted)
		ratio_sorted_1 = np.divide(inner_dot_site_sorted_1, inner_dot_node_sorted_1)

		location_sorted=[]
		location_sorted_1=[]
		for i in range(len(ratio_sorted)):
			if round(ratio_sorted[i],2)==1:
				location_sorted.append(i)
			if round(ratio_sorted_1[i],2)==1:
				location_sorted_1.append(i)
		if len(connection_node)>8 and len(connection_node)<24:
			location_sortednan=[]
			location_sortednan_1=[]
			for i in range(len(ratio_sorted)):
				if np.isnan(ratio_sorted[i])==True or round(ratio_sorted[i],2)==0:
					location_sortednan.append(i)
				if np.isnan(ratio_sorted_1[i])==True:
					location_sortednan_1.append(i)
		
		if len(location_sorted)==1:
			location_sorted=[]
			location_sorted.append(0)
			difference = []
			for i in range(1, len(ratio_sorted)):
				if ratio_sorted[i] < 1:
					difference.append(10000)
				elif ratio_sorted[i] >1:
					difference.append(abs(1 - ratio_sorted[i]))
			index_ratio = min(enumerate(difference), key=itemgetter(1))[0]
			location_sorted.append(index_ratio +1)

		tfflag=0

		if len(connection_node)>10 and len(connection_node)<24: ## to deal with fcu and ftw
			if len(location_sortednan)<2:
				connection_node_spec = [connection_node[location_sorted[0]], connection_node[location_sorted[1]], np.cross(connection_node[location_sorted[0]], connection_node[location_sorted[1]])]
				connection_site_spec = [sorted_sites[location_sorted[0]], sorted_sites[location_sorted[1]], np.cross(sorted_sites[location_sorted[0]], sorted_sites[location_sorted[1]])]
			elif len(location_sortednan)>=2:
				connection_node_spec = [connection_node[location_sorted[0]], connection_node[location_sorted[1]], connection_node[location_sortednan[0]], connection_node[location_sortednan[1]]]
				connection_site_spec = [sorted_sites[location_sorted[0]], sorted_sites[location_sorted[1]], sorted_sites[location_sortednan[0]], sorted_sites[location_sortednan[1]]]
				
			
			connection_node_spec = np.asarray(connection_node_spec, dtype=np.float64)
			connection_site_spec = np.asarray(connection_site_spec, dtype=np.float64)
					
						
			tfflag=1
		
		if len(connection_node)>12: ## to deal with rht
			connection_node = [connection_node[location_sorted[0]], connection_node[location_sorted[1]], connection_node[location_sorted_1[0]], connection_node[location_sorted_1[1]]]
			connection_site = [sorted_sites[location_sorted[0]], sorted_sites[location_sorted[1]], sorted_sites_1[location_sorted_1[0]], sorted_sites_1[location_sorted_1[1]]]
			
		if tfflag==0:## if number of connection points in topology is 8 or less, except *bct*.  
			perm = np.asarray(list(itertools.permutations(connection_site)))#permutations of connection sites
			
			
			node_site_distance=[]
			for i in range(len(perm)):
				trans_matrix = superimposition_matrix(np.transpose(perm[i]), np.transpose(connection_node), usesvd=False)
				perm_plus_one = np.append(perm[i], np.ones([len(perm[i]),1]),1)
				trial_sites=[]
				for k in range(len(perm_plus_one)):
					test_sites=np.dot(trans_matrix, perm_plus_one[k])
					trial_sites.append(test_sites)
				perm_sites = np.asarray(trial_sites)
				perm_sites = perm_sites[:, :-1]
				site_distance=[]
				for k in range(len(perm_sites)):
					site_distance.append(np.linalg.norm(perm_sites[k]-connection_node[k]))
				node_site_distance.append(sum(site_distance))
				#if node_site_distance[i] < 1:
					#break
			
			index_perm = min(enumerate(node_site_distance), key=itemgetter(1))[0]#pick permutation that fits best
		elif tfflag==1:# if connection points in topology is more than 8, except *bct*.  Number of connection points has been decreased.
			perm = np.asarray(list(itertools.permutations(connection_site_spec)))
			node_site_distance=[]
			for i in range(len(perm)):
				trans_matrix = superimposition_matrix(np.transpose(perm[i]), np.transpose(connection_node_spec), usesvd=False)
				perm_plus_one = np.append(perm[i], np.ones([len(perm[i]),1]),1)
				trial_sites=[]
				for k in range(len(perm_plus_one)):
					test_sites=np.dot(trans_matrix, perm_plus_one[k])
					trial_sites.append(test_sites)
				perm_sites=np.asarray(trial_sites)
				perm_sites = perm_sites[:, :-1]
				site_distance=[]
				for k in range(len(perm_sites)):
					site_distance.append(np.linalg.norm(perm_sites[k] - connection_node_spec[k]))
				node_site_distance.append(sum(site_distance))
			index_perm = min(enumerate(node_site_distance), key=itemgetter(1))[0]#pick permutation that fits best

		connection_site = perm[index_perm]
		
		#print index_perm
		##Calculate transformation matrix, using quaternions, to map building block vectors onto network vectors
		if tfflag==0:
			tfmatrix = superimposition_matrix(np.transpose(connection_site), np.transpose(connection_node), usesvd=False)
		elif tfflag==1:
			tfmatrix = superimposition_matrix(np.transpose(connection_site), np.transpose(connection_node_spec), usesvd=False)

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
			frac_element_coord = frac_element_coord + vertex_coord[j]
			tf_frac_element_coord.append(frac_element_coord)
		
		tf_frac_element_coord = np.asarray(tf_frac_element_coord) #building block after transformation, in frac coords
		
		bb_elements.append(elements)
		bb_frac_coordinates.append(tf_frac_element_coord)
		bb_connectivity.append(connectivity)
	bb_elements = np.asarray(bb_elements)
	bb_frac_coordinates = np.asarray(bb_frac_coordinates)
	bb_connectivity=np.asarray(bb_connectivity)
	return bb_elements, bb_frac_coordinates, bb_connectivity











