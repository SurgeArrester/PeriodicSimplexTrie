import os
os.chdir('/home/cameron/Dropbox/University/PhD/SimplexTree')

import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import numpy as np

from PeriodicSimplexTrie import *
from PointcloudGenerator import PointCloudGenerator

lattice = {}
lattice['a'] = 1
lattice['b'] = 2
lattice['alpha'] = 90

motif = ((4.75, 4.75), (5.25, 5.25))

x = PointCloudGenerator(motif, lattice)

points = x.pointset

points = np.array(points)
tri = Delaunay(points)

simplices = tri.simplices
simplex_points = points[tri.simplices]

x.fit_simplices(simplex_points)

trie = PeriodicSimplexTrie()

print()

for simplex in x.simplex_array:
    trie.insert(simplex)

print("---------------")
print(trie.string_rep())
print("---------------")

full_filt = trie.get_filtration()

for k, v in full_filt.items():
    for simp in v:
        print(f"{simp} -> {k}")
print()

plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
plt.plot(points[:,0], points[:,1], 'o')
plt.show()
