include ../../Makefile.inc

CFLAGS+= -I$(PWD)

run:
	./tensor_test --test_jac

build: user.hpp tensor_test.o
	$(CXX) -o tensor_test tensor_test.o $(LDLIBS)

tensor_test.o: $(HEADERS) $(SRC_DIR)/tensor_test.cpp
	$(CXX) $(CFLAGS) -c $(SRC_DIR)/tensor_test.cpp

clean:
	$(RM) *.o tensor_test


