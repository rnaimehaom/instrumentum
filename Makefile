OBJECTS = instrumentum.o global.o random.o node.o grid.o molecule.o molecule_disk.o molecule_ring.o molecule_substitution.o molecular_assembler.o

CXX = g++ -std=c++11

OPT = -O3 -march=native -fstrict-aliasing -ffast-math -ftree-vectorize -funroll-loops

DEBUG = -g -Wall -DDEBUG -DVERBOSE

CXX_FLAGS = $(OPT)
#CXX_FLAGS = $(DEBUG)

LD_FLAGS = $(CXX_FLAGS) -pthread

LIBS = -lsqlite3 -lpugixml -lm

instrumentum: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o instrumentum $(OBJECTS) $(LIBS)

instrumentum.o: instrumentum.cxx instrumentum.h molecular_assembler.h
	$(CXX) $(CXX_FLAGS) -c instrumentum.cxx

global.o: global.cxx instrumentum.h
	$(CXX) $(CXX_FLAGS) -c global.cxx

grid.o: grid.cxx grid.h node.h molecule.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c grid.cxx

node.o: node.cxx node.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c node.cxx

random.o: random.cxx random.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c random.cxx

molecule.o: molecule.cxx molecule.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c molecule.cxx

molecule_disk.o: molecule_disk.cxx molecule.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c molecule_disk.cxx

molecule_ring.o: molecule_ring.cxx molecule.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c molecule_ring.cxx

molecule_substitution.o: molecule_substitution.cxx molecule.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c molecule_substitution.cxx

molecular_assembler.o: molecular_assembler.cxx molecular_assembler.h grid.h molecule.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c molecular_assembler.cxx

docs:
	doxygen docs.config

clean:
	rm -f $(OBJECTS)
	rm -f instrumentum
	rm -rf documentation/*













