OBJECTS = instrumentum.o global.o random.o node.o grid.o molecule.o molecular_assembler.o

CXX = g++ -std=c++17 -fopenmp

OPT = -O3 -fstrict-aliasing -ffast-math -ftree-vectorize -funroll-loops

DEBUG = -g -Wall -DVERBOSE -DDEBUG

#CXX_FLAGS = $(OPT)  
CXX_FLAGS = $(DEBUG)

LIBS = -lsqlite3 -lm

instrumentum: $(OBJECTS)
	$(CXX) $(CXX_FLAGS) -o instrumentum $(OBJECTS) $(LIBS)  

instrumentum.o: instrumentum.cpp instrumentum.h molecular_assembler.h
	$(CXX) $(CXX_FLAGS) -c instrumentum.cpp

global.o: global.cpp instrumentum.h
	$(CXX) $(CXX_FLAGS) -c global.cpp

grid.o: grid.cpp grid.h node.h molecule.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c grid.cpp

node.o: node.cpp node.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c node.cpp

random.o: random.cpp random.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c random.cpp

molecule.o: molecule.cpp molecule.h instrumentum.h 
	$(CXX) $(CXX_FLAGS) -c molecule.cpp

molecular_assembler.o: molecular_assembler.cpp molecular_assembler.h grid.h molecule.h instrumentum.h
	$(CXX) $(CXX_FLAGS) -c molecular_assembler.cpp

clean:
	rm -f $(OBJECTS)
	rm -f instrumentum













