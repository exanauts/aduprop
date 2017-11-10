CFLAGS=-O3 -Wall -g -std=c++11 -I${CODI_DIR}/include

all: powerad

powerad: power.o
	$(CXX) -o powerad power.o

power.o: power.cpp
	$(CXX) $(CFLAGS) -c power.cpp

clean:
	$(RM) *.o


