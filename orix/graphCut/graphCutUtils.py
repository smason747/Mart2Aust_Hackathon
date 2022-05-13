import numpy as np
import matplotlib.pyplot as plt
import maxflow

###TODO make graph cut object

class graphCut(object):
    def __init__(self):
        self.g = maxflow.GraphFloat()


def graphCut(nodeDataArr, ip_weight, NROWS, NCOLS_ODD, GRID):

#ip_weight = 0.01 ### max > ip_weight > min of data as a place to start

### start graphcut

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

    return img2