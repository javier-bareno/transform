# transform
Python2 module for crystallographic transformations.

transform defines the UCell class, which represents a Bravais lattice by its basis vectors. This can be done by secifying the lattice parametters (a, b, c, alpha, beta, gamma) in the constructor, or by passing specific basis vectors to the UCell.from_basis_vectors. 

''' Python
rocksalt = UCell(a,b,c)
tetrag = UCell().from_basis_vectors()
'''
