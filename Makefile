include Makefile.inc

all: lib

# debug options
debug: CFLAGS += -DDBUG
debug: all

# optimized
opt: CFLAGS := $(filter-out -O0,$(CFLAGS))
opt: CFLAGS += -O3
opt: all

# tests
tests: 
	cd test ; for i in `ls` ; do cd $$i ; make test ; done

# LIBRARY

lib: alg.o $(HEADERS) 
	$(AR) $(AROPT) libaduprop.a alg.o

alg.o: $(SRC_DIR)/alg.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c $(SRC_DIR)/alg.cpp

# DOCUMENTATION

doc: $(HEADERS) doc/Doxyfile.in doc/mainpage.md
	cd doc ; doxygen Doxyfile.in ; cd ..

# CLEANING
 
clean:
	$(RM) *.o libaduprop.a

cleanall:
	$(RM) *.o libaduprop.a
	cd test ; for i in `ls` ; do cd $$i ; make clean ; cd .. ; done
	cd examples ; for i in `ls` ; do cd $$i ; make clean ; cd .. ; done

# BUILD EXAMPLES

examples:
	cd examples ; for i in `ls` ; do cd $$i ; make ; cd .. ; done

# BUILDALL

all: lib examples tests

.PHONY: clean cleanall all examples tests lib doc





