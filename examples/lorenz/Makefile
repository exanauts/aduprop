include ../../Makefile.inc

lorenz: user.hpp aduprop.o

CFLAGS+= -I$(PWD)

lorenz: user.hpp aduprop.o
	$(CXX) -o lorenz aduprop.o $(LDLIBS)

aduprop.o: $(SRC_DIR)/aduprop.cpp
	$(CXX) $(CFLAGS) -c $(SRC_DIR)/aduprop.cpp


