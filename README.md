Example Usage
-------------
```python
from PeriodicSimplexTrie import *

x = PeriodicSimplexTrie()

a = Vector([1, 0])
b = Vector([0, 2])
y = Simplex([("p", a), ("p", b), ("p", -a-b)])

x.insert(y)

print(x.string_rep())

full_filtration = x.get_filtration()

for filt_level, simplex_list in full_filtration.items():
    for simplex in simplex_list:
        print(f"{simplex} -> {filt_level}")
```
