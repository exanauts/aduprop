include ../../Makefile.inc

CFLAGS+= -I$(PWD)

lorenz: user.hpp main.o
	$(CXX) -o lorenz main.o $(LDLIBS)

main.o: main.cpp
	$(CXX) $(CFLAGS) -c main.cpp

clean:
	$(RM) *.o lorenz


