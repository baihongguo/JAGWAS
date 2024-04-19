CXX=g++
CXXFLAGS = -std=c++11 -O3 -g -Wall

# Note that armadillo is dependent on LAPACK/BLAS
LIBS= -larmadillo

C_SOURCES= $(/HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/machine/JAGWAS/thirdparty/R_pgamma/*.c)
CXX_SOURCES= JAGWAS.cpp
C_OBJECTS = $(C_SOURCES:.c=.o)
CXX_OBJECTS= JAGWAS.o
OBJS= $(CXX_OBJECTS) $(C_OBJECTS)

#targets
all: JAGWAS

JAGWAS: $(CXX_OBJECTS) $(C_OBJECTS)
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJS)

