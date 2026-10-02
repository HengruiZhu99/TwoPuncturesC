# Makefile for TwoPunctures

BASE=$(shell /bin/pwd)
SRCD=$(BASE)/src
INCD=$(BASE)/include
OBJD=$(BASE)/obj
LIBD=$(BASE)/lib


NAME=TwoPuncturesRun
LIBNAME=TwoPunctures
EXE=$(NAME).x
SRC=$(wildcard $(SRCD)/*.c)
OBJ=$(patsubst $(SRCD)/%.c,$(OBJD)/%.o,$(SRC))
LIBOBJ=$(filter-out $(OBJD)/TwoPuncturesRun.o,$(OBJ))
INC=$(wildcard $(INCD)/*.h)
INC_PARAMS=$(foreach d, $(INCD), -I$d)

LIB=$(LIBD)/lib$(LIBNAME).so
STATIC_LIB=$(LIBD)/lib$(LIBNAME).a
TEST_EXE=$(BASE)/tests/test_physical_api.x


# mandatory flags

CC = gcc
LD = ld
AR = ar

CFLAGS = -std=c99 -fPIC -pedantic $(shell gsl-config --cflags)
LFLAGS = $(shell gsl-config --libs)

CFLAGS += -O3

# old flag information [Leave alone]
##CFLAGS + =`gsl-config --cflags`# GSL
##CFLAGS += -fopenmp # activate OMP (opt)
##LFLAGS=`gsl-config --libs`


all: $(EXE) $(LIB) $(STATIC_LIB)
	@echo "All done"

$(EXE): $(OBJ)
	@echo "Building $@ ..."
	@echo

	$(CC) $(CFLAGS) $(OBJ) $(LFLAGS) -o $@
        
#	# this doesn't need to be a shared object
#	@rm $(OBJD)/TwoPuncturesRun.o > /dev/null 2>&1

$(OBJD)/%.o: $(SRCD)/%.c $(INC)
	@echo "Building objects ..."
	@echo

	@mkdir -p $(dir $@)
#	@echo "Compiling $< ... -> $@"
	$(CC) $(CFLAGS) $(INC_PARAMS) -c $< -o $@

$(LIB): $(LIBOBJ) Makefile
	@echo "Making libraries... "
	@echo

	@mkdir -p $(LIBD)
	$(CC) $(CFLAGS) -shared $(LIBOBJ) $(LFLAGS) -o $@

$(STATIC_LIB): $(LIBOBJ) Makefile
	@mkdir -p $(LIBD)
	@rm -f $@
	$(AR) rcs $@ $(LIBOBJ)
#	$(AR) rcs $(LIBD)/libTwoPunctures_static.a $@ $^

$(TEST_EXE): $(BASE)/tests/test_physical_api.c $(STATIC_LIB) $(INC)
	$(CC) $(CFLAGS) $(INC_PARAMS) $< $(STATIC_LIB) $(LFLAGS) -o $@

test: $(TEST_EXE)
	$(TEST_EXE)

clean:
	@echo "Cleaning ..."
	@rm -rf $(OBJD)
	@rm -rf $(LIBD)
	@rm -rf $(EXE)
	@rm -f $(TEST_EXE)
	@echo "... done"

.PHONY: all clean test

# Independent HiSpID backend. No edits to legacy equations/global parameters.
CXX ?= c++
HISPID_DIR = $(BASE)/build-hispid
HISPID_CPP = $(SRCD)/HiSpID_geometry.cpp $(SRCD)/HiSpID_solver.cpp
HISPID_OBJ = $(patsubst $(SRCD)/%.cpp,$(HISPID_DIR)/%.o,$(HISPID_CPP))
HISPID_LIB = $(HISPID_DIR)/libHiSpID.so
HISPID_FLAGS = -std=c++17 -O3 -fPIC -Wall -Wextra $(shell gsl-config --cflags)

$(HISPID_DIR)/%.o: $(SRCD)/%.cpp $(INCD)/HiSpID.h $(SRCD)/HiSpID_jets.hpp $(SRCD)/HiSpID_internal.hpp $(SRCD)/HiSpID_spectral.hpp
	@mkdir -p $(HISPID_DIR)
	$(CXX) $(HISPID_FLAGS) $(INC_PARAMS) -c $< -o $@

$(HISPID_LIB): $(HISPID_OBJ) $(STATIC_LIB)
	$(CXX) -shared $(HISPID_OBJ) $(STATIC_LIB) $(LFLAGS) -o $@

hispid: $(HISPID_LIB)

$(HISPID_DIR)/test_geometry.x: tests/test_hispid_geometry.cpp $(HISPID_LIB)
	$(CXX) $(HISPID_FLAGS) $(INC_PARAMS) $< $(HISPID_OBJ) $(STATIC_LIB) $(LFLAGS) -o $@

test-hispid-native: $(HISPID_DIR)/test_geometry.x
	$(HISPID_DIR)/test_geometry.x

PYTHON ?= python3
test-hispid: test-hispid-native
	OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=$(BASE)/python:$(BASE)/validation:$(BASE)/examples HISPID_LIBRARY=$(HISPID_LIB) $(PYTHON) -m unittest discover -s tests -p test_hispid.py -v

clean-hispid:
	rm -rf $(HISPID_DIR)

.PHONY: hispid test-hispid-native test-hispid clean-hispid
