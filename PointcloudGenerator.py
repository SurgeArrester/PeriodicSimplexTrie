import os
os.chdir("/home/cameron/Dropbox/University/PhD/SimplexTree")

import numpy as np
from math import sqrt, cos
from PeriodicSimplexTrie import *

class PointCloudGenerator():
    def __init__(self, motif, lattice):
        self.a = lattice['a']
        self.b = lattice['b']
        self.alpha = lattice['alpha']

        self.motif = motif
        self.num_motif = len(motif)

        self.generate_supercell()

    def generate_supercell(self):
        w = self.a

        if self.alpha == 90:
            h = self.b

        else:
            h = cos(self.alpha) * self.b

        pointset = list()
        labelled_pointset = list()
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                next_supercell = [(x+(h*i), y+(w*j)) for (x, y) in self.motif]
                next_labels = range(len(self.motif))
                pointset.append(next_supercell)
                labelled_pointset.append(list(zip(next_supercell, next_labels)))
                print()
        # Flatten the nested lists
        pointset = [item for sublist in pointset for item in sublist]
        labelled_pointset = [item for sublist in labelled_pointset for item in sublist]
        self.pointset = pointset
        self.labelled_pointset = labelled_pointset

        return pointset

    def fit_simplices(self, simplices):
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

            simplex_array.append(simplex_rep)
            print()

        self.simplex_array = simplex_array


if __name__ == "__main__":
    lattice = {}
    lattice['a'] = 1
    lattice['b'] = 2
    lattice['alpha'] = 90

    motif = ((4.75, 4.75),(5.25, 5.25))

    x = PointCloudGenerator(motif, lattice)
