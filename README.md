This software project creates a set of small organic molecules on a diamond-based grid, i.e. beginning with a pure carbon 
layout. The molecule is then subjected to a series of random operations to make it more "medicinal" and ensure that it can 
bind well to a particular pharmacophore. The output molecules are written to an SQLite database in MDL MOL format, after 
which OpenBabel can be used to minimize their geometry via a Python script. The software supposes that the development version 
of the SQLite library is installed and its Python scripts use Python 3 along with the tkinter, sqlite3 and openbabel modules. 
