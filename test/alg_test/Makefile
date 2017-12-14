include ../../Makefile.inc

test: alg_test

alg_test: alg_test.o
	$(CXX) -o alg_test alg_test.o $(LDLIBS)

alg_test.o: alg_test.cpp 
	$(CXX) $(CFLAGS) -c alg_test.cpp

clean:
	$(RM) *.o alg_test
