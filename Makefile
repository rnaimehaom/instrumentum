OBJECTS = instrumentum.o random.o node.o grid.o molecule.o molecular_assembler.o

CXX = g++

DEBUG = -g -Wall

OPENMP = -fopenmp

#CXX_FLAGS += $(CXX_OPT) $(OPENMP) 
CXX_FLAGS += $(DEBUG) $(OPENMP)

#LD_FLAGS += $(CXX_OPT) $(OPENMP)
LD_FLAGS += $(DEBUG) $(OPENMP)

LIBS = -lboost_filesystem -lboost_system -lsqlite3 -lm

instrumentum: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o instrumentum $(OBJECTS) $(LIBS)  

instrumentum.o: instrumentum.cpp instrumentum.h molecular_assembler.h
	$(CXX) $(CXX_FLAGS) -c instrumentum.cpp

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













