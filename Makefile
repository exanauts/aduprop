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

# CLEANING
 
clean:
	$(RM) *.o libaduprop.a

cleanall:
	$(RM) *.o libaduprop.a
	cd test ; for i in `ls` ; do cd $$i ; make clean ; done
	cd examples ; for i in `ls` ; do cd $$i ; make clean ; done

# BUILD EXAMPLES

examples:
	cd examples ; for i in `ls` ; do cd $$i ; make ; done

# BUILDALL

all: lib examples tests

.PHONY: clean cleanall all examples tests lib doc





