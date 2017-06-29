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
import contextlib



#Write .cif file of structure 



def __write_cif(template, node_1, node_2, edge, unit_cell, elements_and_frac_coord, connectivity):
 
 
 a = unit_cell[0]
 b = unit_cell[1]
 c = unit_cell[2]
 
 length_a = np.linalg.norm(a)
 length_b = np.linalg.norm(b)
 length_c = np.linalg.norm(c)
 
 alpha= round(np.arccos(np.dot(b,c)/(np.linalg.norm(b)*np.linalg.norm(c)))*180/np.pi)
 beta= round(np.arccos(np.dot(a, c)/(np.linalg.norm(a)*np.linalg.norm(c)))*180/np.pi)
 gamma= round(np.arccos(np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)))*180/np.pi)


 if not os.path.exists("output_structures"):
  os.makedirs("output_structures")
 with open("output_structures/"+template[-50:-9]+'_'+node_1[-50:-4]+'_'+node_2[-50:-4]+'_'+edge[-50:-4]+".cif", "w") as out_file:
  out_file.write("""data_top_down
  _audit_creation_date              2014-05-22
  _audit_creation_method            'Top_Down_generator'
  _symmetry_space_group_name_H-M    'P1'
  _symmetry_Int_Tables_number       1
  _symmetry_cell_setting            triclinic
  loop_
  _symmetry_equiv_pos_as_xyz""" + '\n'
  ' ' + "x,y,z" + '\n'
  "_cell_length_a" +' ' +      str("{0:.4f}".format(length_a)) + '\n'
  "_cell_length_b" + ' ' + str("{0:.4f}".format(length_b)) + '\n'
  "_cell_length_c" + ' ' + str("{0:.4f}".format(length_c)) + '\n'
  "_cell_angle_alpha" + ' ' + str(alpha) + '\n'
  "_cell_angle_beta" + ' ' + str(beta) + '\n'
  "_cell_angle_gamma" + ' ' + str(gamma) + '\n'
  """loop_
  _atom_site_label
  _atom_site_type_symbol
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z""" + '\n'
  + str(elements_and_frac_coord) + '\n'
  """loop_
  _geom_bond_atom_site_label_1
  _geom_bond_atom_site_label_2
  _geom_bond_distance
  _geom_bond_site_symmetry_2
  _ccdc_geom_bond_type""" + '\n'
  + str(connectivity) + '\n')	
