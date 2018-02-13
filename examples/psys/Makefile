include ../../Makefile.inc

CFLAGS+= -I$(PWD)

.PHONY: clean

power: user.hpp main.o $(HEADERS) 
	$(CXX) -o power main.o $(LDLIBS)

main.o: main.cpp
	$(CXX) $(CFLAGS) -c main.cpp

clean:
	$(RM) *.o power


