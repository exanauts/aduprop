include Makefile.inc

ifndef CODI_DIR
 $(error Environment variable CODI_DIR is undefined)
endif

CFLAGS += -std=c++11 -I$(CODI_DIR)/include
LDLIBS = $(MATH_LIBS)

HEADERS = ad.hpp user.hpp alg.hpp linsolve.hpp tensor.hpp

.PHONY: all debug tests doc clean

# HDF5 support
ifneq ($(HDF_INSTALL),)
	CFLAGS += -I$(HDF_INSTALL)/include
	CFLAGS += -DHDF5
	LDLIBS += -L$(HDF_INSTALL)/lib
	LDLIBS += $(HDF_INSTALL)/lib/libhdf5.a
	LDLIBS += $(HDF_INSTALL)/lib/libhdf5_hl.a
	LDLIBS += -lsz -lz -lm
endif

all: powerad

# debug options
debug: CFLAGS += -DDBUG
debug: all

# optimized
opt: CFLAGS := $(filter-out -O0,$(CFLAGS))
opt: CFLAGS += -O3
opt: all

# tests
tests: alg_test

# MAIN PROGRAM

powerad: power.o alg.o
	$(CXX) -o powerad power.o alg.o $(LDLIBS)

power.o: power.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c power.cpp

alg.o: alg.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c alg.cpp

# TESTS

alg_test: alg_test.o alg.o
	$(CXX) -o alg_test alg_test.o alg.o $(LDLIBS)

alg_test.o: alg_test.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c alg_test.cpp

# DOCUMENTATION

doc: power.cpp ad.hpp user.hpp doc/Doxyfile.in doc/mainpage.md
	cd doc ; doxygen Doxyfile.in ; cd ..

clean:
	$(RM) *.o
	$(RM) powerad


