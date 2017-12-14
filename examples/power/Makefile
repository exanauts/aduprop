include ../../Makefile.inc

CFLAGS+= -I$(PWD)

power: user.hpp aduprop.o

power: user.hpp aduprop.o
	$(CXX) -o power aduprop.o $(LDLIBS)

aduprop.o: $(SRC_DIR)/aduprop.cpp
	$(CXX) $(CFLAGS) -c $(SRC_DIR)/aduprop.cpp

clean:
	$(RM) *.o power


