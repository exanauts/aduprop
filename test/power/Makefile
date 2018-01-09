include ../../Makefile.inc

CFLAGS+= -I$(PWD)

run:
	./power --test_jac
	./power --tensor1
	./power --tensor2
	./power --tensor3

build: user.hpp power.o
	$(CXX) -o power power.o $(LDLIBS)

power.o: $(HEADERS) power.cpp
	$(CXX) $(CFLAGS) -c power.cpp

clean:
	$(RM) *.o power 


