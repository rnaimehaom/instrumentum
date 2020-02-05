This C++ software project creates a set of small organic molecules on a diamond-based grid, i.e. beginning with a pure carbon 
layout, and consists of roughly 5,500 lines of code. The molecules are subjected to a series of random operations to make them 
more "medicinal" and ensure they can bind well to a particular pharmacophore. The output molecules are initially written to an 
[SQLite database](sqlite.org) in MDL MOL format, after which the [Open Babel](openbabel.org) software can be used to minimize 
their geometry via a Python 3 script. The software supposes that the development version of the SQLite library is installed and 
its Python scripts depend on the <code>tkinter</code>, <code>sqlite3</code> and <code>openbabel</code> modules. The C++ code 
has been parallelized using the threading constructs introduced in the 2011 C++ standard along with other aspects of modern C++ 
and it therefore requires a compiler compatible with the 2014 C++ standard. The software reads its parameters from a plain text 
file, an example of which is the <code>parameters.txt</code> file. The Instrumentum software is released under the GNU Public 
License version 3.0; see the <code>LICENSE.txt</code> file or the [Free Software Foundation](www.fsf.org/licensing) for more 
details.    
