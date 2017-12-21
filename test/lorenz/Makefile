include ../../Makefile.inc

CFLAGS+= -I$(PWD)

run:
	./lorenz --test_jac
	./lorenz --tensor1
	./lorenz --tensor2
	./lorenz --tensor3

build: user.hpp lorenz.o
	$(CXX) -o lorenz lorenz.o $(LDLIBS)

lorenz.o: $(HEADERS) lorenz.cpp
	$(CXX) $(CFLAGS) -c lorenz.cpp

clean:
	$(RM) *.o lorenz 


