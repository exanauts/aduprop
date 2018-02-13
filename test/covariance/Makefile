include ../../Makefile.inc

CFLAGS+= -I$(PWD)

run:
	./power

build: user.hpp main.o
	$(CXX) -o power main.o $(LDLIBS)

main.o: $(HEADERS) main.cpp
	$(CXX) $(CFLAGS) -c main.cpp

clean:
	$(RM) *.o power output.txt 


