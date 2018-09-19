include ../../Makefile.inc

CFLAGS+= -I$(PWD)

all: user.hpp main.o
	$(CXX) -o plorenz main.o $(LDLIBS)

main.o: $(HEADERS) main.cpp
	$(CXX) $(CFLAGS) -c main.cpp

clean:
	$(RM) *.o plorenz datas data_cov1 data_cov2 data_cov3 


