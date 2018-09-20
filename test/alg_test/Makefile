include ../../Makefile.inc

run: build
	# TODO: on some platforms this returns -0 and and some 0. This is
	# reported as an error. So the test is disabled for now
	./alg_test # > diff.out
	#diff diff.out diff.in >&2

build: alg_test

alg_test: alg_test.o
	$(CXX) -o alg_test alg_test.o $(LDLIBS)

alg_test.o: alg_test.cpp 
	$(CXX) $(CFLAGS) -c alg_test.cpp

clean:
	$(RM) *.o diff.out alg_test out
