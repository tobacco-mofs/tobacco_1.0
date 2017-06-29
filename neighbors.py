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
import numpy as np
def neighbor_edges(uc, v1, nv1, e1, ne1, cn, hl):
    tolerance = 0.1
    e1_nabors = []
    for i in range(0,nv1) :
       v = np.asarray([v1[i,0], v1[i,1], v1[i,2]])
       for j in range(0, ne1):
          e = np.asarray([e1[j,0], e1[j,1], e1[j,2]])
          diffa = v[0] - e[0]
          diffb = v[1] - e[1]
          diffc = v[2] - e[2]

          if diffa > 0.5:
             e[0] = e[0] + 1.0 
          elif diffa < -0.5:
             e[0] = e[0] - 1.0 

          if diffb > 0.5:
             e[1] = e[1] + 1.0
          elif diffb < -0.5:
             e[1] = e[1] - 1.0
 
          if diffc > 0.5:
             e[2] = e[2] + 1.0
          elif diffc < -0.5:
             e[2] = e[2] - 1.0
          v_xyz = np.dot(np.transpose(uc), v)
          e_xyz = np.dot(np.transpose(uc), e)

          distance = np.linalg.norm(v_xyz - e_xyz)
          if distance < (hl + tolerance) and distance > (hl - tolerance):
             e1_nabors.append([j])
    e1_nabors = np.asarray(e1_nabors)
    e1_nabors = np.reshape(e1_nabors, (nv1, cn))
    
    return e1_nabors

def neighbor_vertices(uc, v1, nv1, e1, ne1, hl):
    tolerance = 0.1
    v1_nabors = []
    for i in range(0,ne1) :
       e = np.asarray([e1[i,0], e1[i,1], e1[i,2]])
       for j in range(0, nv1):
          v = np.asarray([v1[j,0], v1[j,1], v1[j,2]])
          diffa = e[0] - v[0]
          diffb = e[1] - v[1]
          diffc = e[2] - v[2]

          if diffa > 0.5:
             v[0] = v[0] + 1.0
          elif diffa < -0.5:
             v[0] = v[0] - 1.0

          if diffb > 0.5:
             v[1] = v[1] + 1.0
          elif diffb < -0.5:
             v[1] = v[1] - 1.0
        
          if diffc > 0.5:
             v[2] = v[2] + 1.0 
          elif diffc < -0.5:
             v[2] = v[2] - 1.0
          v_xyz = np.dot(np.transpose(uc), v)
          e_xyz = np.dot(np.transpose(uc), e)
    
          distance = np.linalg.norm(v_xyz - e_xyz)
          if distance < (hl + tolerance) and distance > (hl - tolerance):
             v1_nabors.append([j])
         
    v1_nabors = np.asarray(v1_nabors)

    return v1_nabors

