CXX=g++
CXXFLAGS = -std=c++11 -O3 -g -Wall

# Note that armadillo is dependent on LAPACK/BLAS
LIBS= -larmadillo

CXX_SOURCES= JAGWAS.cpp
CXX_OBJECTS= JAGWAS.o
OBJS= $(CXX_OBJECTS)

#targets
all: JAGWAS

JAGWAS: $(CXX_OBJECTS) 
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm -f $(OBJS)

