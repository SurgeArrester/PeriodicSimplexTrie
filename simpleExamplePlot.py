"""
Author: Cameron Hargreaves

A small test script to print the simplex tree and associated filtration for a
pointset defined by the user
"""

import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import numpy as np

from PeriodicSimplexTrie import *
from PointcloudGenerator import PointCloudGenerator

lattice = {}
lattice['a'] = 1
lattice['b'] = 2
lattice['alpha'] = 90

single_point = ((0, 0),)
dual_point = ((-0.25, -0.25), (0.25, 0.25))
benzene = ((0, 0.5), (0.25, 0), (0.75, 0), (1, 0.5), (0.25, 1), (0.75, 1))

motif = dual_point

x = PointCloudGenerator(motif, lattice)

points = np.array(x.pointset)
tri = Delaunay(points)

simplices = tri.simplices
simplex_points = points[tri.simplices]

x.fit_simplices(simplex_points)

trie = PeriodicSimplexTrie()

for simplex in x.simplex_array:
    trie.insert(simplex)

print("\nPeriodic Simplex Trie Representation")
print("---------------")
print(trie.string_rep())

full_filt = trie.get_filtration()

print("\nFull Filtration")
print("---------------")
label_list = {}
i = 0
for k, v in full_filt.items():
    for simp in v:
        if str(simp) in label_list.keys():
            label_list[str(simp)] = i
            i += 1

        print(f"{simp} -> {k}")

plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
plt.plot(points[:,0], points[:,1], 'o')
plt.show()

