# transform
Python2 module for crystallographic transformations.

transform defines the UCell class, which represents a Bravais lattice by its basis vectors. This can be done by secifying the lattice parametters (a, b, c, alpha, beta, gamma) in the constructor, or by passing specific basis vectors to the UCell.from_basis_vectors method. Upon instance creation, UCell willl calculate the reciprocal basis vectors.

For example, to create the convencional rock salt (NaCl) unit cell:

```python
a_NaCl = 5.6402 # Ångstroms
rock_salt = UCell(a_NaCl, a_NaCl, a_NaCl)
print(rock_salt)

Direct basis:
a= 5.6402       b= 5.6402       c= 5.6402
alpha= 90.00    beta= 90.00     gamma= 90.00
vol= 179.4252

Reciprocal basis:
a= 0.1773       b= 0.1773       c= 0.1773
alpha= 90.00    beta= 90.00     gamma= 90.00
vol= 0.0056
```

UCell objects provide a series of (self explanatory) methods to calculate distances, angles, etc. Most accept one or two numpy.array's parameters, which are interpreted as the coordinates of direct (directions) or reciprocal (planes) vectors in the associated basis. You can also transform coordinates to and from the canonical (cartessian) basis. For example, to constrauct the primitive unit cell of the rock salt structure:

```python
a_prim = rock_salt.dir_to_cart(np.array([0.5, 0, 0.5]))
b_prim = rock_salt.dir_to_cart(np.array([0.5, 0.5, 0]))
c_prim = rock_salt.dir_to_cart(np.array([0, 0.5, 0.5]))
primitive = UCell().from_basis_vectors(a_prim, b_prim, c_prim)
print(primitive)

Direct basis:
a= 3.9882       b= 3.9882       c= 3.9882
alpha= 60.00    beta= 60.00     gamma= 60.00
vol= 44.8563

Reciprocal basis:
a= 0.3071       b= 0.3071       c= 0.3071
alpha= 109.47   beta= 109.47    gamma= 109.47
vol= 0.0223
alpha= 90.00    beta= 90.00     gamma= 90.00
```

To tranform coordinates between bases, just transform to and from cartessian. For example:

```python
rock_salt.plane_from_cart(primitive.plane_to_cart(np.array([1,1,0])))

array([ 2.00000000e+00, -3.40823177e-17, -3.40823177e-17])

```

Finally, the derived class Hexag is specialized to hexagnal unit cells. It includes additional methods to work with directions and planes in 4-coordinate notation.

Feel free to look through the code and experiment. Enjoy!
