# -*- coding: utf-8 -*-
"""
Created on Fri May 13 09:59:31 2022

@author: ashley

gets the connectivity of a network or a sparse adjacency matrix
"""

from diffpy.structure import Atom, Lattice, Structure
import numpy as np
from orix.crystal_map import CrystalMap
from orix.quaternion.rotation import Rotation

import tempfile

from diffpy.structure import Atom, Lattice, Structure
import matplotlib.pyplot as plt
import numpy as np

import maxflow
from orix import data, io, plot
from orix.crystal_map import CrystalMap, Phase, PhaseList
from orix.quaternion import Orientation, Rotation, symmetry
from orix.vector import Vector3d
plt.close('all')

#%%
################## load materials and image 
# this whole section will be replaced by graph cut function eventually
path = r'/Users/paytone/Documents/GitHub/maxflow_for_matsci/Data/steel_ebsd.ang'

# Read each column from the file
euler1, euler2, euler3, x, y, iq, dp, phase_id, sem, fit  = np.loadtxt(path, unpack=True)

# Create a Rotation object from Euler angles
euler_angles = np.column_stack((euler1, euler2, euler3))
rotations = Rotation.from_euler(euler_angles)

# Create a property dictionary
properties = dict(iq=iq, dp=dp)

# Create unit cells of the phases
structures = [
    Structure(
        title="ferrite",
        atoms=[Atom("fe", [0] * 3)],
        lattice=Lattice(0.287, 0.287, 0.287, 90, 90, 90)
    ),
]
phase_list = PhaseList(
    names=["ferrite"],
    point_groups=["432"],
    structures=structures,
)

# Create a CrystalMap instance
xmap2 = CrystalMap(
    rotations=rotations,
    phase_id=phase_id,
    x=x,
    y=y,
    phase_list=phase_list,
    prop=properties,
)
xmap2.scan_unit = "um"


ckey_m3m = plot.IPFColorKeyTSL(xmap2.phases["ferrite"].point_group, direction=Vector3d.zvector())
rgb_fe = ckey_m3m.orientation2color(xmap2["ferrite"].orientations)

fer_x = np.round(2*(xmap2['ferrite'].x))
fer_y = np.round(2*(xmap2['ferrite'].y))

#%% start graph cut
ipw = 0.01  #inplane weight
g = maxflow.GraphFloat()
nodeids = g.add_grid_nodes((305, 305)) #int(np.sqrt(len(dp)))

# arrange x and y corrdinates for ferrite phase and rgb values 
for_network = fer_x, fer_y, rgb_fe[:,0], rgb_fe[:,1], rgb_fe[:,2]
for_network = np.asarray(for_network).T

#replace extra phase with zeros and place ferrite phase in correct node spots in image
img = np.zeros((305,305))
for xx in range(len(for_network)):
    coordx, coordy = int(for_network[xx,0]), int(for_network[xx,1])
    img[coordx, coordy] = for_network[xx,2] #only using r channel for now--will be replaced with somehting else later
img = img.T #it gets turned around

# start creating the network for graph cut
g.add_grid_edges(nodeids, ipw, symmetric=True)
g.add_grid_tedges(nodeids, img, 1-img)
g.maxflow() #get graph cut


#%% get connectivity information from network
C = g.get_nx_graph()

import networkx as nx
import scipy.sparse as sp
AA = nx.to_scipy_sparse_array(C)
sparseUpperArr = sp.triu(AA)


u,v,wt = sp.find(sparseUpperArr)
connectivity = np.asanyarray([u,v])

#%%
sink = np.amax(connectivity.ravel())
source = sink-1
source_edges = np.any(connectivity==source,axis=0)
sink_edges = np.any(connectivity==sink, axis=0)
out_plane_loc = np.any(np.vstack([source_edges, sink_edges]), axis=0)
connectivity2 = connectivity[:, ~out_plane_loc]

#%%
o1 = xmap2.rotations[connectivity2[0,:]]
o2 = xmap2.rotations[connectivity2[1,:]]
m = (~o1).outer(o2)

#%%
m2 = Orientation(m, symmetry.Oh).map_into_symmetry_reduced_zone()
D2 = m2.angle

'''
NEXT STEP: ASSIGN D2 AS EDGE WEIGHTS, THEN PERFORM GRAPH CUT ON IT
'''


#%% plot graph cut
sgm = g.get_grid_segments(nodeids)
img2 = np.int_(np.logical_not(sgm))
from matplotlib import pyplot as ppl
ppl.imshow(img2)
ppl.show()


