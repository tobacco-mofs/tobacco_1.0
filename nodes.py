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
import math



## This code will read a building block cif file for a node. 
## It will identify unit cell length, connection sites, distances from 
## center to connection site, angles between connection sites.
## It will also identify atoms and connectivity to be used later.






def __node_properties(arg):
	grab_lines = False
	grab_connectivity = False
	atom_info=[]
	connectivity_info=[]
	element_index=[]
	element=[]
	x_frac_coord=[]
	y_frac_coord=[]
	z_frac_coord=[]
	bond_site1=[]
	bond_site2=[]
	bond_distance=[]
	bond_symm=[]
	bond_type=[]
	connection_site=[]
	distance_connection_site=[]
	angle_connection_site_pair=[]
	with open("nodes_bb/"+arg, 'r') as node_cif:
		node = node_cif.read()
	for line in node.split("\n"):
		if "_cell_length_a" in line:
			length_a = line.split()[1]#unit cell vector
			length_a =float(length_a)#string to float
		if "_cell_length_b" in line:
			length_b = line.split()[1]
			length_b = float(length_b)
		if "_cell_length_c" in line:
			length_c= line.split()[1]
			length_c= float(length_c)
		if "_cell_angle_alpha" in line:
			alpha = line.split()[1]
			alpha = float(alpha)
		if "_cell_angle_beta" in line:
			beta= line.split()[1]
			beta= float(beta)
		if "_cell_angle_gamma" in line:
			gamma = line.split()[1]
			gamma = float(gamma)
		if line.startswith('_atom_site_occupancy'): #grabbing atoms and coordinates
			grab_lines = True
			continue
		if line.startswith('loop_'):
			grab_lines = False
		if grab_lines:
			atom_info.append(line)
		if line.startswith('_ccdc_geom_bond_type'): # grabbing connectivity
			grab_connectivity = True
			continue
		if line.startswith(' '):
			grab_connectivity = False
		if grab_connectivity:
			connectivity_info.append(line) 
	for line in atom_info:
		split = line.split()
		atom_element_index = np.asarray([split[0]]) # element and index
		atom_element= np.asarray([split[1]]) # element 
		x_frac = np.asarray([split[2]]) # x frac coord
		y_frac = np.asarray([split[3]]) # y frac coord
		z_frac = np.asarray([split[4]]) # z frac coord
		element_index.append(atom_element_index) #build the list with index info
		element.append(atom_element) #build list with element info
		x_frac_coord.append(x_frac) # build list of x frac coords
		y_frac_coord.append(y_frac) #build list of y frac coords
		z_frac_coord.append(z_frac) #build list of z frac coords
	
	for line in connectivity_info:
		split = line.split()
		site1 = np.asarray([split[0]])
		site2 = np.asarray([split[1]])
		distance = np.asarray([split[2]])
		symm = np.asarray([split[3]])
		type = np.asarray([split[4]])
		bond_site1.append(site1)
		bond_site2.append(site2)
		bond_distance.append(distance)
		bond_symm.append(symm)
		bond_type.append(type)
		
		
		
	# create unit cell vectors from a, b, c. alpha, beta, gamma information
	ax = length_a
	ay = 0.0
	az = 0.0
	bx = length_b * np.cos(gamma * np.pi / 180.0)
	by = length_b * np.sin(gamma * np.pi / 180.0)
	bz = 0.0
	cx = length_c * np.cos(beta * np.pi / 180.0)
	cy = length_c * length_b * math.cos(alpha * np.pi /180.0) - bx * cx / by
	cz = (length_c ** 2 - cx ** 2 - cy ** 2) ** 0.5
	
	unit_cell =  np.asarray([[ax, ay, az],[bx, by, bz], [cx, cy, cz]])
##### TURN LISTS INTO ARRAYS AND BUILD MATRIX
	element_index=np.asarray(element_index) 
	element=np.asarray(element)
	x_frac_coord=np.asarray(x_frac_coord)
	y_frac_coord=np.asarray(y_frac_coord)
	z_frac_coord=np.asarray(z_frac_coord)
	elements = np.column_stack([element_index, element]) 
	element_frac_coord= np.column_stack([x_frac_coord, y_frac_coord, z_frac_coord])#combine arrays to make matrix
	element_frac_coord=element_frac_coord.astype(float)
	
	element_coord=[]
	for i in range(len(element_frac_coord)):
		coord = np.dot(np.transpose(unit_cell), element_frac_coord[i])
		element_coord.append(coord)
	element_coord=np.asarray(element_coord)

	connectivity=np.column_stack([bond_site1, bond_site2, bond_distance, bond_symm, bond_type]) #connectivity data
	
	
	if len(element) == 1:
		connection_site=[]
		distance_connection_site=[0.0]
		angle_connection_site_pair=[]
		connectivity=[]
		return connection_site, distance_connection_site, angle_connection_site_pair, connectivity, elements, element_coord




	###Identify connection sites.
	for i in [i for i, x in enumerate(element_index) if x[0][0]=='X']:
		connection=element_coord[i]
		connection_site.append(connection)
	connection_site=np.asarray(connection_site)
	
	
	centroid_x=[]
	centroid_y=[]
	centroid_z=[]
	for j in range(len(connection_site)):
		centroid_x.append(connection_site[j][0])
		centroid_y.append(connection_site[j][1])
		centroid_z.append(connection_site[j][2])
	centroid_x = sum(centroid_x)/len(connection_site)
	centroid_y = sum(centroid_y)/len(connection_site)
	centroid_z = sum(centroid_z)/len(connection_site)
		
	centroid = np.column_stack([centroid_x, centroid_y, centroid_z])
	
	if len(connection_site)==2:
		centroid_x=[]
		centroid_y=[]
		centroid_z=[]
		for j in range(len(element_coord)):
			centroid_x.append(element_coord[j][0])
			centroid_y.append(element_coord[j][1])
			centroid_z.append(element_coord[j][2])
		centroid_x=sum(centroid_x)/len(element_coord)
		centroid_y=sum(centroid_y)/len(element_coord)
		centroid_z=sum(centroid_z)/len(element_coord)
		
		centroid=np.column_stack([centroid_x, centroid_y, centroid_z])
	
	
	for i in range(len(element_coord)):
		element_coord[i] = element_coord[i]-centroid
	for i in range(len(connection_site)):
		connection_site[i] = connection_site[i] - centroid

	
	for i in range(len(connection_site)):
		distance = np.linalg.norm(connection_site[i] - 0)
		distance_connection_site.append(distance)
	distance_connection_site=np.asarray(distance_connection_site)

	
	
	## Angles between connection points

	
	for pair in itertools.combinations(connection_site, 2):

		angle = round(np.arccos(np.dot(*pair)/(np.linalg.norm(pair[0])*np.linalg.norm(pair[1])))*180/np.pi)
		angle_connection_site_pair.append(angle)
	angle_connection_site_pair=np.asarray(angle_connection_site_pair)

	return connection_site, distance_connection_site, angle_connection_site_pair, connectivity, elements, element_coord


