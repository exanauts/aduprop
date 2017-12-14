include Makefile.inc

all: lib

# debug options
debug: CFLAGS += -DDBUG
debug: all

# tests
tests: 
	cd test ; for i in `ls` ; do cd $$i ; make test ; done

# LIBRARY

lib: alg.o 
	$(AR) $(AROPT) libaduprop.a alg.o

alg.o: $(SRC_DIR)/alg.cpp
	$(CXX) $(CFLAGS) -c $(SRC_DIR)/alg.cpp

# DOCUMENTATION

doc: doc/Doxyfile.in doc/mainpage.md
	cd doc ; doxygen Doxyfile.in ; cd ..

clean:
	$(RM) *.o


