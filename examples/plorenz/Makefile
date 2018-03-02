include ../../Makefile.inc

CFLAGS+= -I$(PWD)

all: user.hpp main.o
	$(CXX) -o plorenz main.o $(LDLIBS)

main.o: $(HEADERS) main.cpp
	$(CXX) $(CFLAGS) -c main.cpp

clean:
	$(RM) *.o plorenz output.txt 


