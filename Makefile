CXX = clang++
CFLAGS = -O0 -Wall -g -std=c++11 -I${CODI_DIR}/include
LDLIBS = -llapack -lblas


HEADERS = ad.hpp user.hpp alg.hpp linsolve.hpp


all: powerad

tests: alg_test

# MAIN PROGRAM

powerad: power.o alg.o
	$(CXX) -o powerad power.o alg.o $(LDLIBS)

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

doc: power.cpp ad.hpp user.hpp doc/Doxyfile.in doc/mainpage.md
	cd doc ; doxygen Doxyfile.in ; cd ..

clean:
	$(RM) *.o


