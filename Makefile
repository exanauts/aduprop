CFLAGS=-g -O0 -Wall -g -std=c++11 -I${CODI_DIR}/include

CXX=clang++

all: powerad

powerad: power.o
	$(CXX) -o powerad power.o -llapack 

power.o: power.cpp linsolve.hpp
	$(CXX) $(CFLAGS) -c power.cpp

clean:
	$(RM) *.o


