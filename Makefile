CFLAGS=-g -O0 -Wall -g -std=c++11 -I${CODI_DIR}/include

CXX=clang++

all: powerad

powerad: power.o
	$(CXX) -o powerad power.o -llapack -lblas

power.o: power.cpp ad.hpp
	$(CXX) $(CFLAGS) -c power.cpp

doc: power.cpp ad.hpp doc/Doxyfile.in doc/mainpage.md
	cd doc ; doxygen Doxyfile.in ; cd ..

clean:
	$(RM) *.o


