CXX=g++
CXXFLAGS = -std=c++11 -O3 -g -Wall

# Note that armadillo is dependent on LAPACK/BLAS
LIBS= -larmadillo

C_SOURCES= $(/HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/machine/JAGWAS/JAGWAS/JAGWAS/thirdparty/R_pgamma/*.c) $(/HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/machine/JAGWAS/JAGWAS/JAGWAS/thirdparty/R_pgamma/thirdparty/*.h) $(/HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/machine/JAGWAS/JAGWAS/JAGWAS/thirdparty/R_pgamma/R_ext/*.h)
CXX_SOURCES= JAGWAS.cpp
CXX_OBJECTS= JAGWAS.o
OBJS= $(CXX_OBJECTS)

#targets
all: JAGWAS

JAGWAS: $(CXX_OBJECTS) $(C_OBJECTS)
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJS)

