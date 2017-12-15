include ../../Makefile.inc

CFLAGS+= -I$(PWD)

power: user.hpp aduprop.o
	$(CXX) -o power aduprop.o $(LDLIBS)

aduprop.o: $(HEADERS) $(SRC_DIR)/aduprop.cpp
	$(CXX) $(CFLAGS) -c $(SRC_DIR)/aduprop.cpp

clean:
	$(RM) *.o power


