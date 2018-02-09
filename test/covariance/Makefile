include ../../Makefile.inc

CFLAGS+= -I$(PWD)

run:
	mkdir -p output
	./power

build: user.hpp main.o
	$(CXX) -o power main.o $(LDLIBS)

main.o: $(HEADERS) main.cpp
	$(CXX) $(CFLAGS) -c main.cpp

clean:
	$(RM) *.o power 


