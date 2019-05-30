"""
Implement a trie with insert, search, and startsWith methods.
Note:
You may assume that all inputs are consist of lowercase letters a-z.
"""
import collections
import reprlib
from math import sqrt
from copy import deepcopy

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
    Custom dict that prints a bit prettier
    """
    def __repr__(self):
        children = len(self.keys())
        return_str = f"This node has {children} children: \n"
        for k, v in self.items():
            return_str += f"{k} : {v}\n"
        return return_str

class SimplexTrieNode:
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
        ret = "\t"*level + str(self.simplex) + "\n"

        for _, node in self.children.items():
            ret += node.print_node(level=level+1)

        return ret

    def apply_filtration(self, filt):
        """
        Recursively go through each child and return each list that passes the
        filtration
        """
        global simplex

        for _, node in self.children.items():
            current_simp = Simplex(node.simplex)
            if hasattr(current_simp.vertices[0][1], "magnitude"):
                magnitude = current_simp.vertices[0][1].magnitude
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
        for _, node in self.children.items():
            current_simp = Simplex([(node.simplex)])
            if hasattr(current_simp.vertices[0][1], "magnitude"):
                filt_points.add(current_simp.vertices[0][1].magnitude)
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
    generic vector class, input a list of dimensional labels
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
        if (isinstance(other, Vector)):
            other = other.dims

        if other is None and self is not None:
            return False

        for i, dim in enumerate(self.dims):
            if (str(dim) != str(other[i])):
                return False
        return True

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
        """
        Sort based on magnitude, then on string representation if mags are the
        same
        """
        if isinstance(other, Vector):
            if self.magnitude != other.magnitude:
                return self.magnitude < other.magnitude
            else:
                return str(self.dims) < str(other.dims)

class Simplex:
    def __init__(self, vertices):
        # Sort the vertices based on their labels? DEPRECATED
        if isinstance(vertices, list):
            self.vertices = sorted(vertices, key=lambda x: x[1])
        else:
            self.vertices = [vertices]


    def __getitem__(self, index):
        return self.vertices[index]

    def __len__(self):
        return len(self.vertices)

    def __iter__(self):
        for vertex in self.vertices:
            yield vertex
        return

    def __str__(self):
        x = []
        if self.vertices is None:
            return "Simplex()"
        for vertex in self.vertices:
            x.append(vertex)
        return f"Simplex({x})"

    def __repr__(self):
        x = []
        if self.vertices is None:
            return "Simplex()"
        for vertex in self.vertices:
            x.append(vertex)
        return f"Simplex({reprlib.repr(x)})"

    def __eq__(self, other):
        if other == None or other == [] or hasattr(other, "vertices") == False or other[0] is None:
            return False

        for i, vertex in enumerate(self.vertices):
            if vertex != other.vertices[0][i]:
                return False
        return True

    def __hash__(self):
        return hash(self.vertices)

    def append(self, other):
        this = self.vertices

        if isinstance(other, Simplex):
            other = other.vertices
            if isinstance(other, tuple):
                this = this + [other]
            else:
                this = this + other

            self.vertices = sorted(this, key=lambda x: x[1])
        return

    def clear(self):
        self.vertices = []

class PeriodicSimplexTrie:
    def __init__(self):
        self.root = SimplexTrieNode()

    def insert(self, labelled_simplices):
        current = self.root
        level = 0
        for labelled_simplex in labelled_simplices:
            current.children[labelled_simplex].parent = current
            current = current.children[labelled_simplex]
            current.simplex = labelled_simplex
            current.label = labelled_simplex[0]
            current.vector_to_parent = labelled_simplex[1]
            current.level = level
            level += 1
        current.is_word = True

    def search(self, other_simplex):
        current = self.root
        for labelled_simplices in other_simplex:
            current = current.children.get(labelled_simplices)
            if current is None:
                return False
        return current.is_word

    def starts_with(self, prefix):
        current = self.root
        for labelled_simplex in prefix:
            current = current.children.get(labelled_simplex)
            if current is None:
                return False
        return True

    def string_rep(self):
        """
        Prints the nested tree and the unique labels
        """
        x = ""
        for label, node in self.root.children.items():
            x += (node.print_node())
        return x

    def filtration(self, filt):
        current = self.root
        filtration = current.apply_filtration(filt)
        filtration = [x for x in filtration if x is not None]
        if len(filtration) == 0:
            point_list = [k for k, v in current.children]
            filtration = list(set(point_list))
        return filtration

    def get_filtration_points(self):
        filt_points = set()
        filt_points.add(0)
        current = self.root
        current.add_filt_points(filt_points)
        self.filt_points = sorted(filt_points)
        return sorted(filt_points)

    def get_filtration(self):
        if not hasattr(self, "filt_points"):
            _ = self.get_filtration_points()

        filt = {}

        for p in self.filt_points:
            filt[p] = self.filtration(p)

        points = list(reversed(self.filt_points))

        for i, p in enumerate(points):
            if i != len(points) - 1:
                filt[p] = [x for x in filt[p] if x not in filt[points[i+1]]]

        self.full_filtration = filt

        return filt

    def __str__(self):
        return self.string_rep()

if __name__ == "__main__":
    main()
