include ../../Makefile.inc

power: user.hpp aduprop.o

CFLAGS+= -I$(PWD)

power: user.hpp aduprop.o
	$(CXX) -o lorenz aduprop.o $(LDLIBS)

aduprop.o: $(SRC_DIR)/aduprop.cpp
	$(CXX) $(CFLAGS) -c $(SRC_DIR)/aduprop.cpp


