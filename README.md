This C++ software project creates a set of small organic molecules on a diamond-based grid, i.e. beginning with a pure carbon 
layout, and consists of roughly 5,500 lines of code. The molecules are subjected to a series of random operations to make them 
more "medicinal" and ensure they can bind well to a particular pharmacophore. The output molecules are initially written to an 
[SQLite database](https://www.sqlite.org) in MDL MOL format, after which the [Open Babel](http://www.openbabel.org) software can 
be used to minimize their geometry via a Python 3 script, <code>post-processing.py</code>. The software supposes that the 
development version of the SQLite 3 library is installed and the post-processing script depends on the <code>sqlite3</code> and 
<code>openbabel</code> modules. The C++ source code has been parallelized using the threading constructs introduced in the 2011 
C++ standard along with other aspects of modern C++ and it therefore requires a compiler compatible with the 2014 C++ standard. 
The Instrumentum software reads its parameters from an XML file, an example of which is the included <code>parameters.xml</code> 
file, and there is also a Python 3 script, <code>editor.py</code>, which provides a portable graphical interface based on 
[Tkinter](https://en.wikipedia.org/wiki/Tkinter) for editing these parameter files which are used by the software. The Instrumentum 
software is released under the GNU Public License version 3.0; see the <code>LICENSE.txt</code> file or the 
[Free Software Foundation](https://www.fsf.org/licensing) for more details.    
