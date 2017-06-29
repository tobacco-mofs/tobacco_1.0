Manual for ToBaCCo
==================

ToBaCCo (**To**pologically-**Ba**sed **C**rystal **Co**nstructor) is a code to make crystalline structures based on a given topology, nodes, and edges.  This current version of the code only deals with edge-transitive nets, which were obtained from RCSR.  Briefly, the code reads a net, chooses appropriate nodes and edges, scales the unit cell to fit the building blocks and subsequently, places and connects them, generating a crystal in .cif format. It uses python 2.7. .


Topology
========

A topology, or net, template file contains the name of the topology, the unit cell vectors, the number of types of vertices, or nodes, as well as the number of types of edges. It also contains the coordination and symmetry key information of each type of node as well as the total number and the fractional coordinates of each building block. Below is an example of a net template:

```
topology SHE 
unit_cell_vectors
28.284 0.0000  0.0000
0.0000  28.284 0.0000
0.00000  0.0000 28.284               
  
2 types_of_vertices  
type number coordination  symmetry_key
1       12        4             5    
2       8         6             0       

1 types_of_edges   
type   number half_length
12      48         5.00

fract_coordinates_vertices
first type 
0.2500   0.0000   0.5000
0.7500   0.0000   0.5000
0.2500   0.5000   0.0000
0.5000   0.2500   0.0000
0.5000   0.7500   0.0000
0.0000   0.5000   0.2500
0.0000   0.5000   0.7500
0.5000   0.0000   0.2500
0.0000   0.2500   0.5000
0.0000   0.7500   0.5000
0.7500   0.5000   0.0000
0.5000   0.0000   0.7500

second type
0.2500   0.2500   0.2500
0.7500   0.7500   0.2500
0.7500   0.2500   0.7500
0.2500   0.7500   0.7500
0.2500   0.2500   0.7500
0.7500   0.7500   0.7500
0.2500   0.7500   0.2500
0.7500   0.2500   0.2500

fract_coordinates_edges
first type
0.2500   0.1250   0.3750
0.7500   0.8750   0.3750
0.7500   0.1250   0.6250
0.2500   0.8750   0.6250
0.3750   0.2500   0.1250
0.3750   0.7500   0.8750
0.6250   0.7500   0.1250
0.6250   0.2500   0.8750
0.1250   0.3750   0.2500
0.8750   0.3750   0.7500
0.1250   0.6250   0.7500
0.8750   0.6250   0.2500
0.1250   0.2500   0.6250
0.8750   0.7500   0.6250
0.1250   0.7500   0.3750
0.8750   0.2500   0.3750
0.2500   0.3750   0.8750
0.7500   0.3750   0.1250
0.7500   0.6250   0.8750
0.2500   0.6250   0.1250
0.3750   0.1250   0.7500
0.3750   0.8750   0.2500
0.6250   0.1250   0.2500
0.6250   0.8750   0.7500
0.7500   0.8750   0.6250
0.2500   0.3750   0.1250
0.2500   0.1250   0.6250
0.2500   0.8750   0.3750
0.7500   0.1250   0.3750
0.6250   0.7500   0.8750
0.1250   0.2500   0.3750
0.6250   0.2500   0.1250
0.3750   0.2500   0.8750
0.3750   0.7500   0.1250
0.8750   0.6250   0.7500
0.3750   0.1250   0.2500
0.1250   0.6250   0.2500
0.8750   0.3750   0.2500
0.1250   0.3750   0.7500
0.8750   0.7500   0.3750
0.8750   0.2500   0.6250
0.1250   0.7500   0.6250
0.7500   0.6250   0.1250
0.2500   0.6250   0.8750
0.7500   0.3750   0.8750
0.6250   0.8750   0.2500
0.6250   0.1250   0.7500
0.3750   0.8750   0.7500
```


Building Blocks
===============

A building block can be an edge or a node.  They should be in a .cif file and the connection points should be marked with X’s in the first column.  It is important to note that *these files cannot have an empty line at the end of the file*.  Below is an example of a .cif file for a node.

```cif
data_sym_4_mc_1
_audit_creation_date              2014-09-12
_audit_creation_method            'Materials Studio'
_symmetry_space_group_name_H-M    'P1'
_symmetry_Int_Tables_number       1
_symmetry_cell_setting            triclinic
loop_
_symmetry_equiv_pos_as_xyz
  x,y,z
_cell_length_a                    20.0000
_cell_length_b                    20.0000
_cell_length_c                    20.0000
_cell_angle_alpha                 90.0000
_cell_angle_beta                  90.0000
_cell_angle_gamma                 90.0000
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
_atom_site_adp_type
_atom_site_occupancy
C1     C    -0.10182   0.05995  -0.07511   0.00000  Uani   1.00
H2     H    -0.09897   0.09775  -0.03653   0.16000  Uiso   1.00
C3     C    -0.14698   0.06723  -0.12789   0.00000  Uani   1.00
H4     H    -0.17950   0.11037  -0.13046   0.17800  Uiso   1.00
C5     C    -0.06270  -0.04251  -0.12036   0.00000  Uani   1.00
H6     H    -0.02914  -0.08471  -0.11793   0.16000  Uiso   1.00
C7     C    -0.10721  -0.03695  -0.17393   0.00000  Uani   1.00
H8     H    -0.10868  -0.07489  -0.21245   0.17800  Uiso   1.00
N9     N    -0.06133   0.00531  -0.07188   0.00000  Uani   1.00
C10    C    -0.00999  -0.01534   0.13924   0.00000  Uani   1.00
H11    H     0.04400  -0.01632   0.13603   0.16000  Uiso   1.00
C12    C    -0.04143  -0.02215   0.20123   0.00000  Uani   1.00
H13    H    -0.01220  -0.02817   0.24634   0.17800  Uiso   1.00
C14    C    -0.11509  -0.00723   0.08480   0.00000  Uani   1.00
H15    H    -0.14374  -0.00246   0.03931   0.16000  Uiso   1.00
C16    C    -0.14832  -0.01393   0.14589   0.00000  Uani   1.00
H17    H    -0.20237  -0.01366   0.14792   0.17800  Uiso   1.00
Zn18   Zn   -0.00049   0.00005   0.00051   0.00000  Uani   1.00
N19    N    -0.04701  -0.00732   0.08267   0.00000  Uani   1.00
C20    C     0.05389   0.12110   0.05226   0.00000  Uiso   1.00
N21    N     0.05195   0.07889  -0.00115   0.00000  Uiso   1.00
C22    C     0.08556   0.09318  -0.05860   0.00000  Uiso   1.00
C23    C     0.12421   0.15116  -0.06329   0.00000  Uiso   1.00
C24    C     0.09200   0.17959   0.04944   0.00000  Uiso   1.00
C25    C     0.02427  -0.13758  -0.00349   0.00000  Uiso   1.00
C26    C     0.06223  -0.19604  -0.00802   0.00000  Uiso   1.00
C27    C     0.16156  -0.12851  -0.01938   0.00000  Uiso   1.00
C28    C     0.12193  -0.07116  -0.01468   0.00000  Uiso   1.00
N29    N     0.05442  -0.07668  -0.00761   0.00000  Uiso   1.00
H30    H     0.02565   0.10988   0.09688   0.00000  Uiso   1.00
H31    H     0.08232   0.05935  -0.10069   0.00000  Uiso   1.00
H32    H     0.15129   0.16253  -0.10871   0.00000  Uiso   1.00
H33    H     0.09390   0.21319   0.09178   0.00000  Uiso   1.00
H34    H    -0.02932  -0.14058   0.00321   0.00000  Uiso   1.00
H35    H     0.03844  -0.24453  -0.00514   0.00000  Uiso   1.00
H36    H     0.21517  -0.12444  -0.02522   0.00000  Uiso   1.00
H37    H     0.14526  -0.02251  -0.01622   0.00000  Uiso   1.00
X38    C    -0.11109  -0.02135   0.20440   0.00000  Uiso   1.00
X39    C    -0.14961   0.01831  -0.17752   0.00000  Uiso   1.00
X40    C     0.13133  -0.19129  -0.01608   0.00000  Uiso   1.00
X41    C     0.12740   0.19453  -0.00877   0.00000  Uiso   1.00
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
_geom_bond_site_symmetry_2
_ccdc_geom_bond_type
C1     N9      1.362   .     S
C1     C3      1.397   .     D
C1     H2      1.082   .     S
C3     H4      1.082   .     S
C3     X39     1.395   .     S
C5     N9      1.362   .     D
C5     C7      1.397   .     S
C5     H6      1.079   .     S
C7     H8      1.082   .     S
C7     X39     1.395   .     D
N9     Zn18    1.894   .     S
C10    N19     1.362   .     S
C10    C12     1.397   .     D
C10    H11     1.082   .     S
C12    H13     1.082   .     S
C12    X38     1.395   .     S
C14    N19     1.362   .     D
C14    C16     1.397   .     S
C14    H15     1.079   .     S
C16    H17     1.082   .     S
C16    X38     1.395   .     D
Zn18   N19     1.894   .     S
Zn18   N21     1.894   .     S
Zn18   N29     1.894   .     S
C20    N21     1.362   .     A
C20    C24     1.397   .     A
C20    H30     1.080   .     S
N21    C22     1.362   .     A
C22    C23     1.397   .     A
C22    H31     1.082   .     S
C23    H32     1.082   .     S
C23    X41     1.395   .     D
C24    H33     1.082   .     S
C24    X41     1.395   .     S
C25    C26     1.397   .     A
C25    N29     1.362   .     A
C25    H34     1.082   .     S
C26    H35     1.082   .     S
C26    X40     1.395   .     D
C27    C28     1.397   .     A
C27    H36     1.082   .     S
C27    X40     1.395   .     S
C28    N29     1.362   .     A
C28    H37     1.080   .     S
```


Running the code
================

The code is written in python.  In the directory of choice, simply run `main_auto`.  This will iterate through all the topologies found in the folder `templates`, fitting in appropriate nodes and edges found in the folders `nodes_bb` and `edges_bb`, respectively. All these folders should be in the same directory as `main_auto` along with the modules:  `connect_edgeto2nodes`, `connect_nodetonode`, `edges`, `neighbors`, `nodes`, `placing_bb_func`, `placing_edge`, `placing_edge_1node`, `transformations`, and `write2cif`.  If successful, the code will generate a folder called `output_structures` where the .cif files with the structures will be placed.  The name of the structures will contain the topology as well as node, or nodes, and edge used.  
