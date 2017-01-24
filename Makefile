OBJECTS = main.o node.o grid.o molecule.o masm.o global.o 

#CXX_FLAGS += $(CXX_OPT) $(OPENMP) -DBASE
CXX_FLAGS += $(DEBUG) $(OPENMP) -DDBASE

#LD_FLAGS += $(CXX_OPT) $(OPENMP)
LD_FLAGS += $(DEBUG) $(OPENMP)

LIBS = $(BOOST_SYSTEM) -lpqxx -lm

instrumentum: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o instrumentum $(OBJECTS) $(LIBS)  

main.o: main.cpp grid.cpp node.cpp atom.cpp atom.h node.h grid.h
	$(CXX) $(CXX_FLAGS) -c main.cpp

grid.o: grid.cpp atom.cpp node.cpp global.cpp grid.h atom.h node.h global.h
	$(CXX) $(CXX_FLAGS) -c grid.cpp

node.o: node.cpp atom.cpp atom.h node.h
	$(CXX) $(CXX_FLAGS) -c node.cpp

atom.o: atom.cpp atom.h 
	$(CXX) $(CXX_FLAGS) -c atom.cpp

global.o: global.h global.cpp 
	$(CXX) $(CXX_FLAGS) -c global.cpp 

molecule.o: molecule.cpp node.cpp atom.cpp global.cpp molecule.h node.h atom.h global.h 
	$(CXX) $(CXX_FLAGS) -c molecule.cpp

masm.o: grid.cpp masm.cpp molecule.cpp masm.h grid.h molecule.h node.h 
	$(CXX) $(CXX_FLAGS) -c masm.cpp

token.o: token.cpp token.h
	$(CXX) $(CXX_FLAGS) -c token.cpp

clean:
	rm -f $(OBJECTS)
	rm -f instrumentum













