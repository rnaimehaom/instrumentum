OBJECTS = main.o global.o random.o node.o grid.o molecule.o molecular_assembler.o

#CXX_FLAGS += $(CXX_OPT) $(OPENMP) 
CXX_FLAGS += $(DEBUG) $(OPENMP)

#LD_FLAGS += $(CXX_OPT) $(OPENMP)
LD_FLAGS += $(DEBUG) $(OPENMP)

LIBS = $(BOOST_FILESYSTEM) $(BOOST_SYSTEM) -lsqlite3 -lm

instrumentum: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o instrumentum $(OBJECTS) $(LIBS)  

main.o: main.cpp molecular_assembler.h
	$(CXX) $(CXX_FLAGS) -c main.cpp

grid.o: grid.cpp grid.h node.h global.h
	$(CXX) $(CXX_FLAGS) -c grid.cpp

node.o: node.cpp node.h global.h
	$(CXX) $(CXX_FLAGS) -c node.cpp

global.o: global.cpp global.h 
	$(CXX) $(CXX_FLAGS) -c global.cpp 

random.o: random.cpp random.h global.h
	$(CXX) $(CXX_FLAGS) -c random.cpp

molecule.o: molecule.cpp molecule.h node.h global.h 
	$(CXX) $(CXX_FLAGS) -c molecule.cpp

molecular_assembler.o: molecular_assembler.cpp molecular_assembler.h grid.h molecule.h global.h
	$(CXX) $(CXX_FLAGS) -c molecular_assembler.cpp

clean:
	rm -f $(OBJECTS)
	rm -f instrumentum













