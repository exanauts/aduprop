CXX = clang++
CFLAGS = -O0 -Wall -g -std=c++11 -I${CODI_DIR}/include
LDLIBS = -llapack -lblas


HEADERS = linsolve.hpp alg.hpp


all: powerad

tests: alg_test

# MAIN PROGRAM

powerad: power.o
	$(CXX) -o powerad power.o $(LDLIBS)

power.o: power.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c power.cpp

alg.o: alg.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c alg.cpp

# TESTS

alg_test: alg_test.o alg.o
	$(CXX) -o alg_test alg_test.o alg.o $(LDLIBS)

alg_test.o: alg_test.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c alg_test.cpp

# DOCUMENTATION

doc: power.cpp linsolve.hpp doc/Doxyfile.in doc/mainpage.md
	cd doc ; doxygen Doxyfile.in ; cd ..

clean:
	$(RM) *.o


