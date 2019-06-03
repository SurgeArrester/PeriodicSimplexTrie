#!/usr/bin/env python
"""
Author: Cameron Hargreaves

An implementation of the PeriodicSimplexTrie and supporting classes, Vector,
Simplex, and SimplexTrieNode

Example Usage
-------------
x = PeriodicSimplexTrie()
a = Vector([1, 0])
b = Vector([0, 2])
y = Simplex([("p", a), ("p", b), ("p", -a-b)])
x.insert(y)
print(x.string_rep())
full_filtration = x.get_filtration()

"""
import collections
import reprlib
from math import sqrt
from copy import deepcopy

__author__ = "Cameron Hargreaves"
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Cameron Hargreaves"
__email__ = "cameron.h@rgreaves.me.uk"
__status__ = "Development"

simplex = []

def main():
    x = PeriodicSimplexTrie()
    a = Vector([1, 0])
    b = Vector([0, 2])

    y = Simplex([("p", a), ("p", b), ("p", -a-b)])
    x.insert(y)
    x.string_rep()
    print("------------")

    y = Simplex([("p", a+b), ("p", -a), ("p", -b)])
    x.insert(y)
    x.string_rep()
    print("------------")

    y = Simplex([("p", b), ("p", -a), ("p", -b+a)])
    x.insert(y)
    x.string_rep()
    print("------------")

    y = Simplex([("p", b-a), ("p", -b), ("p", a)])
    x.insert(y)
    x.string_rep()
    print("------------")

    y = Simplex([("p", -a), ("p", -b), ("p", b+a)])
    x.insert(y)
    x.string_rep()
    print("------------")

    y = Simplex([("p", -b-a), ("p", a), ("p", b)])
    x.insert(y)
    x.string_rep()
    print("------------")

    y = Simplex([("p", -b), ("p", a), ("p", b-a)])
    x.insert(y)
    x.string_rep()
    print("------------")

    y = Simplex([("p", -b+a), ("p", b), ("p", -a)])
    x.insert(y)

    print(x.string_rep())

    full_filt = x.get_filtration()

    for k, v in full_filt.items():
        for simp in v:
            print(f"{simp} -> {k}")



class NodeDict(collections.defaultdict):
    """
    Custom dict so that print's a bit prettier
    """
    def __repr__(self):
        children = len(self.keys())
        return_str = f"This node has {children} children: \n"
        for k, v in self.items():
            return_str += f"{k} : {v}\n"
        return return_str

class SimplexTrieNode:
    """
    A node of the Simplex Trie
    """
    def __init__(self):
        self.parent = None
        self.children = NodeDict(SimplexTrieNode)
        self.simplex = None
        self.label = None
        self.vector_to_parent = None
        self.is_word = False
        self.level = None
        self.simplex_list = []

    def print_node(self, level=0):
        """Recursively create a string from this node of its children"""
        ret = "\t"*level + str(self.simplex) + "\n"

        for _, node in self.children.items():
            ret += node.print_node(level=level+1)

        return ret

    def apply_filtration(self, filt):
        """
        Recursively go through each child and return each list that passes the
        filtration
        TODO: Remove nasty globals temporary workaround
        """
        global simplex

        for _, node in self.children.items():
            current_simp = Simplex(node.simplex)
            if hasattr(current_simp.labelled_vertices[0][1], "magnitude"):
                magnitude = current_simp.labelled_vertices[0][1].magnitude
            else:
                magnitude = 0

            if len(node.children) == 0 and magnitude <=filt:
                # If we've hit a leaf and
                x = self.get_full_simplex(node)
                simplex.append(x)

            elif magnitude > filt:
                x = self.get_full_simplex(node.parent)
                simplex.append(x)

            else:
                node.apply_filtration(filt)

        return simplex

    def get_full_simplex(self, node, simplex=None):
        """
        Return a full simplex from a point, working it's way
        recursively up through parents
        """
        if node.parent == None:
            if self.simplex == None:
                return
            elif simplex == None:
                simplex = Simplex([(self.simplex)])
            return simplex

        elif simplex == None:
            simplex = Simplex([(node.simplex)])

        elif isinstance(simplex, Simplex):
            simplex.append(Simplex([(node.simplex)]))

        if hasattr(node, "parent"):
            simplex = self.get_full_simplex(node.parent, simplex=simplex)

        return simplex

    def add_filt_points(self, filt_points):
        """
        Recursively call to each nodes children and add the
        magnitudes of each of the vectors for us in future filtrations
        """
        for _, node in self.children.items():
            current_simp = Simplex([(node.simplex)])
            if hasattr(current_simp.labelled_vertices[0][1], "magnitude"):
                filt_points.add(current_simp.labelled_vertices[0][1].magnitude)
            else:
                filt_points.add(0)
            node.add_filt_points(filt_points)
        return

    def __repr__(self):
        return f"SimplexTrieNode: {self.simplex}, children={len(self.children)}"

    def __str__(self):
        return f"SimplexTrieNode(vertex={self.label}, vector to parent={self.vector_to_parent}, children={len(self.children)})"

class Vector:
    """
    Generic vector class, input: list of numerical scalars
    """
    def __init__(self, dims):
        assert isinstance(dims, (list,))
        # TODO double check we actually want magnitude of vectors here
        self.dims = dims # [abs(dim) for dim in  dims]
        self.magnitude = self.mag()

    def mag(self):
        sum = 0
        for x in self.dims:
            sum += x ** 2
        return sqrt(sum)

    def __str__(self):
        return f"Vector({self.dims}, mag={self.magnitude:.2f})"

    def __repr__(self):
        return f"Vector({reprlib.repr(self.dims)})"

    def __hash__(self):
        return hash(tuple(self.dims))

    def __eq__(self, other):
        """Returns True if both vectors have same orientation"""
        if (isinstance(other, Vector)):
            other = other.dims

        if other is None and self is not None:
            return False

        if self.dims == other:
            return True

        elif (-self).dims == other:
            return True

        else:
            return False

    def __add__(self, other):
        if (isinstance(other, Vector)):
            other = other.dims
        return_dim = []
        for i, dim in enumerate(self.dims):
            return_dim.append(dim + other[i])

        return Vector(return_dim)

    def __sub__(self, other):
        if (isinstance(other, Vector)):
            other = other.dims
        return_dim = []
        for i, dim in enumerate(self.dims):
            return_dim.append(dim - other[i])

        return Vector(return_dim)

    def __neg__(self):
        dim = [-x for x in self.dims]
        return Vector(dim)

    def __lt__(self, other):
        """Sort based on magnitude, then on string representation"""
        if isinstance(other, Vector):
            if self.magnitude != other.magnitude:
                return self.magnitude < other.magnitude
            else:
                return str(self.dims) < str(other.dims)

class Simplex:
    """
    Collection of labelled vectors making up a simplex
    Simplex([(label, Vector(to_next_point))])
    Sorts the points in order of magnitude from smallest to largest.

    Keyword arguments:
    labelled_vertices: A tuple or list of tuples in the
                       form [(label, Vector(to_next_point))]
    """
    def __init__(self, labelled_vertices):
        if isinstance(labelled_vertices, list):
            self.labelled_vertices = sorted(labelled_vertices, key=lambda x: x[1])
        else:
            self.labelled_vertices = [labelled_vertices]


    def __getitem__(self, i):
        return self.labelled_vertices[i]

    def __setitem__(self, key, item):
        self.labelled_vertices[key] = item

    def __len__(self):
        return len(self.labelled_vertices)

    def __iter__(self):
        for vertex in self.labelled_vertices:
            yield vertex
        return

    def __str__(self):
        x = []
        if self.labelled_vertices is None:
            return "Simplex()"
        for vertex in self.labelled_vertices:
            x.append(vertex)
        return f"Simplex({x})"

    def __repr__(self):
        x = []
        if self.labelled_vertices is None:
            return "Simplex()"
        for vertex in self.labelled_vertices:
            x.append(vertex)
        return f"Simplex({reprlib.repr(x)})"

    def __eq__(self, other):
        if (other == None) or (other == []) or not (hasattr(other, "labelled_vertices")) or (other[0] is None):
            return False

        for i, vertex in enumerate(self.labelled_vertices):
            if vertex not in other.labelled_vertices:
                return False
        return True

    def __hash__(self):
        Simplex([('p', Vector([1, 0]))]) == Simplex([('p', Vector([1, 0]))])
        return hash(tuple(self.labelled_vertices))

    def append(self, other):
        """Adds labelled vector, (label, Vector([])) to the simplex"""
        this = self.labelled_vertices

        if isinstance(other, Simplex):
            other = other.labelled_vertices
            if isinstance(other, tuple):
                this = this + [other]
            else:
                this = this + other

            self.labelled_vertices = sorted(this, key=lambda x: x[1])
        return

    def clear(self):
        self.labelled_vertices = []

class PeriodicSimplexTrie:
    """
    Periodic simplex trie, has structure of a regular Trie with additional
    methods for calculating the filtration based on the magnitude of the
    vectors to each subsequent point in the Trie
    """
    def __init__(self):
        self.root = SimplexTrieNode()

    def insert(self, labelled_simplices):
        """Insert an ordered, labelled simplex into the Trie"""
        current_node = self.root
        level = 0
        for labelled_simplex in labelled_simplices:
            for existing_label, child in current_node.children.items():
                if labelled_simplex == existing_label:
                    labelled_simplex = existing_label

            # Label the current_node parent for reconstructing the tree
            current_node.children[labelled_simplex].parent = current_node

            # Change to the child and add labels
            current_node = current_node.children[labelled_simplex]
            current_node.simplex = labelled_simplex
            current_node.label = labelled_simplex[0]
            current_node.vector_to_parent = labelled_simplex[1]
            current_node.level = level
            level += 1

        # Mark as a "word" if it is a complete insertion
        current_node.is_word = True

    def search(self, other_simplex):
        """Return True if simplex already exists in the Trie"""
        current = self.root
        for labelled_simplices in other_simplex:
            current = current.children.get(labelled_simplices)
            if current is None:
                return False
        return current.is_word

    def string_rep(self):
        """Returns string of the nested tree structure"""
        x = ""
        for label, node in self.root.children.items():
            x += (node.print_node())
        return x

    def filtration(self, filt):
        """All simplexes that are <= to the filt value passed"""
        current = self.root
        filtration = current.apply_filtration(filt)
        filtration = [x for x in filtration if x is not None]
        if len(filtration) == 0:
            point_list = [k for k, v in current.children]
            filtration = list(set(point_list))

        Simplex([('0', Vector([0, 1])), ('0', Vector([-2, 0]))]) == Simplex([('0', Vector([0, 1])), ('0', Vector([-2, 0]))])
        try:
            filtration[1] == Simplex([('0', Vector([0, 1])), ('0', Vector([-2, 0]))])
        except:
            pass
        return set(filtration)

    def get_filtration_points(self):
        """
        Return the critical points of the
        filtration where new complexes form
        """
        filt_points = set()
        filt_points.add(0)
        current = self.root
        current.add_filt_points(filt_points)
        self.filt_points = sorted(filt_points)
        return sorted(filt_points)

    def get_filtration(self):
        """
        Apply a full filtration to the simplex Trie
        """
        if not hasattr(self, "filt_points"):
            _ = self.get_filtration_points()

        filt = {}

        # Calculate filtrations for each of the points
        for p in self.filt_points:
            filt[p] = self.filtration(p)

        points = list(reversed(self.filt_points))

        # Take out the filtrations from the previous filtration
        for i, p in enumerate(points):
            if i != len(points) - 1:
                filt[p] = [x for x in filt[p] if x not in filt[points[i+1]]]

        self.full_filtration = filt

        return filt

    def __str__(self):
        return self.string_rep()

if __name__ == "__main__":
    main()
