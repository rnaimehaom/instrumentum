OBJECTS = main.o node.o grid.o molecule.o molecular_assembler.o global.o

#CXX_FLAGS += $(CXX_OPT) $(OPENMP) 
CXX_FLAGS += $(DEBUG) $(OPENMP)

#LD_FLAGS += $(CXX_OPT) $(OPENMP)
LD_FLAGS += $(DEBUG) $(OPENMP)

LIBS = -lsqlite3 -lm

instrumentum: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o instrumentum $(OBJECTS) $(LIBS)  

main.o: main.cpp molecular_assembler.h
	$(CXX) $(CXX_FLAGS) -c main.cpp

grid.o: grid.cpp grid.h atom.h node.h global.h
	$(CXX) $(CXX_FLAGS) -c grid.cpp

node.o: node.cpp atom.h node.h global.h
	$(CXX) $(CXX_FLAGS) -c node.cpp

atom.o: atom.cpp atom.h global.h
	$(CXX) $(CXX_FLAGS) -c atom.cpp

global.o: global.cpp global.h 
	$(CXX) $(CXX_FLAGS) -c global.cpp 

molecule.o: molecule.cpp molecule.h node.h atom.h global.h 
	$(CXX) $(CXX_FLAGS) -c molecule.cpp

molecular_assembler.o: molecular_assembler.cpp molecular_assembler.h grid.h molecule.h global.h
	$(CXX) $(CXX_FLAGS) -c molecular_assembler.cpp

clean:
	rm -f $(OBJECTS)
	rm -f instrumentum













