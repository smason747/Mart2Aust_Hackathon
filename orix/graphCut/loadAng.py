import numpy as np
import matplotlib.pyplot as plt

from diffpy.structure import Atom, Lattice, Structure

from orix import plot
from orix.crystal_map import CrystalMap, Phase, PhaseList
from orix.io import load, save
from orix.quaternion import Orientation, Rotation, symmetry
from orix.vector import Vector3d

from glob import glob


def loadAng(filePath):

    strings = ['GRID', 'NCOLS_ODD', 'NCOLS_EVEN', 'NROWS']
    # strings = ['GRID', 'NCOLS_ODD', 'NCOLS_EVEN', 'NROWS', 'Phase', 'MaterialName', 'Formula' 'Symmetry', 'LatticeConstants']
    # secondaryStrings = ['Phase', 'Formula' 'Symmetry', 'LatticeConstants']
    stringsOut = []
    # secondStringsOut = []

    file = open(filePath, 'r')
    lineNum = 0
    line = file.readline()
    while line[0] == '#':
        print(lineNum, line)
        for subs in strings:
            if subs in line:
                start = line.index(':')
                substring = line[start+2:-1] ### TODO this is bad way to get string i want
                stringsOut.append(substring)
        line = file.readline()
        lineNum += 1

    lastLine = file.readline()
    splitLast = lastLine.split(' ')
    splitLast.remove('\n')
    while("" in splitLast) :
        splitLast.remove("")
    print(splitLast)
    print(len(splitLast))

    stringsOut.append(len(splitLast))

    return stringsOut

def getLibraryPhaseList(MaterialName):

    if MaterialName == 'Ferrite':

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
    
    else:
        print('error: no corresponding material in library')

    return phase_list