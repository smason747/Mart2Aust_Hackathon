import numpy as np
import matplotlib.pyplot as plt
import maxflow
from math import sqrt

from diffpy.structure import Atom, Lattice, Structure

from orix import plot
from orix.crystal_map import CrystalMap, Phase, PhaseList
from orix.io import load, save
from orix.quaternion import Orientation, Rotation, symmetry
from orix.vector import Vector3d

from loadAng import loadAng, getLibraryPhaseList

dataPath = 'data/'
# fileName = 'ferrite.ang' ### hexgrid
fileName = 'steel_ebsd.ang' ### squaregrid
target = dataPath + fileName
MaterialName = 'Ferrite'

# xmap = load(dataPath + fileName) ### currently does not handle 

GRID, NCOLS_ODD, NCOLS_EVEN, NROWS, numEBSDColumns = loadAng(target)
NCOLS_ODD = int(NCOLS_ODD)
NCOLS_EVEN = int(NCOLS_EVEN)
NROWS = int(NROWS)
numEBSDColumns = int(numEBSDColumns)

print(GRID, NCOLS_ODD, NCOLS_EVEN, NROWS, numEBSDColumns)

if numEBSDColumns == 10:
    euler1, euler2, euler3, x, y, iq, ci, phase_id, sem, fit = np.loadtxt(target, unpack=True)
else:
    print('TODO handle more/fewer columns')

euler_angles = np.column_stack((euler1, euler2, euler3))
rotations = Rotation.from_euler(euler_angles)
properties = dict(iq=iq, dp=ci)

print(euler_angles.shape)

# #############################################
# ### Hardcoded ORIX crystalmap setup
# # Create unit cells of the phases
# structures = [
#     Structure(
#         title="ferrite",
#         atoms=[Atom("fe", [0] * 3)],
#         lattice=Lattice(0.287, 0.287, 0.287, 90, 90, 90)
#     ),
# ]
# phase_list = PhaseList(
#     names=["ferrite"],
#     point_groups=["432"],
#     structures=structures,
# )

phase_list = getLibraryPhaseList('Ferrite')

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

print(xmap2)

ckey_m3m = plot.IPFColorKeyTSL(xmap2.phases["ferrite"].point_group, direction=Vector3d.zvector())
rgb_fe = ckey_m3m.orientation2color(xmap2["ferrite"].orientations)

################################
### start building structure to feed into network
### build square array of hex grid points for graph cut
### populate pixels with R channel of orientation image
### leave non-indexed or other phases as zeros

if GRID == 'HexGrid':

    nodeDataArr = np.zeros((NROWS, NCOLS_ODD))

    cRow = 0
    cID = 0

    for pair in np.arange(int(NROWS/2)):

        nodeDataArr[cRow, :NCOLS_EVEN] = rgb_fe[cID:cID+NCOLS_EVEN][:,0] ### even
        cRow += 1
        cID = cID+NCOLS_EVEN

        nodeDataArr[cRow, :NCOLS_ODD] = rgb_fe[cID:cID+NCOLS_ODD][:,0] ### odd
        cRow += 1
        cID = cID+NCOLS_ODD

        # print(cRow, pair)

    plt.imshow(nodeDataArr)
    axes=plt.gca()
    axes.set_aspect(0.5)
    plt.show()

elif GRID == 'SqrGrid':
    ### squaregrid data for graphcut

    nodeData2 = np.zeros((xmap2.id.shape[0], 1))
    ferriteIDs =  xmap2['ferrite'].id
    ferriteGray = rgb_fe[:,0]
    failCounter = 0

    for i in xmap2.id:
        # print(i)
        try:
            # print(i)
            idx = np.where(ferriteIDs == i)[0][0]
            nodeData2[i] = ferriteGray[idx]
        except:
            # print('no value', failCounter)
            failCounter += 1

    n = np.round(sqrt(xmap2.id.shape[0])).astype(int)
    nodeDataArr = np.reshape(nodeData2, (n,n))

    plt.imshow(nodeDataArr)
    plt.show()


### start graphcut

ip_weight = 0.01 ### max > ip_weight > min of data as a place to start

g = maxflow.GraphFloat()
if GRID == 'HexGrid':
    nodeids = g.add_grid_nodes((NROWS, NCOLS_ODD))
    structure = np.array([[0, 1, 1],
                      [1, 0, 1],
                      [0, 1, 1]]) ### struct for hex grid
elif GRID == 'SqrGrid':
    nodeids = g.add_grid_nodes((n,n))
    structure = maxflow.vonNeumann_structure(ndim=2, directed=True) ### square grid structure
    
g.add_grid_edges(nodeids, ip_weight)
g.add_grid_tedges(nodeids, nodeDataArr, 1-nodeDataArr)

g.maxflow()

sgm = g.get_grid_segments(nodeids)
print(np.sum(sgm))
img2 = np.int_(np.logical_not(sgm))

plt.imshow(img2)
if GRID == 'HexGrid':
    axes=plt.gca()
    axes.set_aspect(0.5)
plt.show()