"""
Author: Cameron Hargreaves

A simple class to perform a tiling of a two dimensional lattice, label each of
these points, and store them in the variable PointCloudGenerator.simplex_array

Example usage:
import numpy as np
from scipy.spatial import Delaunay

single_point = ((0, 0),)
lattice = {}
lattice['a'] = 1
lattice['b'] = 2
lattice['alpha'] = 90

x = PointCloudGenerator(motif, lattice)
points = np.array(x.pointset)

tri = Delaunay(points)
simplices = tri.simplices
simplex_points = points[tri.simplices]

x.fit_simplices(simplex_points)
"""

import os

import numpy as np
from math import sqrt, cos
from PeriodicSimplexTrie import *

class PointCloudGenerator():
    """
    A tiling class for simple two-dimensional lattices and motifs

    Params
    motif: A list of points within the unit cell as a tuple of tuples
           eg: motif = ((0, 0),)
    lattice: A dictionary of basic 2D lattice parameters, "a", "b", and "alpha"
    """
    def __init__(self, motif, lattice):
        self.a = lattice['a']
        self.b = lattice['b']
        self.alpha = lattice['alpha']

        self.motif = motif
        self.num_motif = len(motif)

        self.generate_supercell()

    def generate_supercell(self):
        """
        Make a supercell of the lattice points that are inputted in the object
        instantiation.
        Additionally create a list of labels for each of the points within the
        supercell, giving them the numeric count of 0-len(motif) and labelling
        these red if they are in the central unit cell of the supercell, and
        black otherwise
        """
        w = self.a
        if self.alpha == 90:
            h = self.b
        else:
            h = cos(self.alpha) * self.b

        pointset = list()
        labelled_pointset = list()

        for i in [-1, 0, 1, 2]:
            for j in [-1, 0, 1, 2]:
                next_supercell = [(x+(h*i), y+(w*j)) for (x, y) in self.motif]

                if i == 0 and j == 0:
                    next_labels = [str(l) + "_red" for l in range(len(self.motif))]
                else:
                    next_labels = [str(l) + "_blk" for l in range(len(self.motif))]

                pointset.append(next_supercell)
                labelled_pointset.append(list(zip(next_supercell, next_labels)))

        # Flatten the nested lists
        pointset = [item for sublist in pointset for item in sublist]
        labelled_pointset = [item for sublist in labelled_pointset for item in sublist]

        self.pointset = pointset
        self.labelled_pointset = labelled_pointset

        return pointset

    def fit_simplices(self, simplices):
        """
        Take each of the points in the labelled pointset and create a simplex
        object with these labels, as long as at least one point is inside the
        central unit cell of the complete supercell.

        There may be multiple repetitions of the same periodic simplex in this
        list
        """
        simplex_array = []
        for simplex in simplices:
            labelled_simplex = []
            for point in simplex:
                point = tuple(point)
                for label in self.labelled_pointset:
                    if point == label[0]:
                        labelled_simplex.append((label))

            simplex_rep = Simplex([])

            for i, simplex in enumerate(labelled_simplex):
                curr_point = simplex[0]
                try:
                    next_point = labelled_simplex[i+1][0]
                except IndexError:
                    next_point = labelled_simplex[0][0]

                vector = [i - j for i, j in zip(next_point, curr_point)]
                vector_to_parent = Vector(vector)
                simplex_rep.append(Simplex((simplex[1], vector_to_parent)))

            # If we only have black points in the simplex pass
            if all(item[0][-4:] == "_blk" for item in simplex_rep):
                pass

            # Else remove the point labels and add this to the list
            else:
                for i, item in enumerate(simplex_rep):
                    simplex_rep[i] = (item[0][:-4], item[1])
                simplex_array.append(simplex_rep)
        self.simplex_array = simplex_array


if __name__ == "__main__":
    lattice = {}
    lattice['a'] = 1
    lattice['b'] = 2
    lattice['alpha'] = 90

    motif = ((4.75, 4.75), (5.25, 5.25))

    x = PointCloudGenerator(motif, lattice)
