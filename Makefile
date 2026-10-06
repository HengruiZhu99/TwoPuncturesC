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
INC=$(wildcard $(INCD)/*.h) $(wildcard $(SRCD)/*.inc)
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
HISPID_FLAGS += $(HISPID_MAP_FLAGS)
HISPID_FLAGS += $(HISPID_EXPERIMENT_FLAGS)

HISPID_KERNEL_HEADERS = $(SRCD)/HiSpID_geometry_types.hpp $(SRCD)/HiSpID_trumpet.hpp $(SRCD)/HiSpID_geometry_kernels.hpp $(SRCD)/HiSpID_cache_kernels.hpp
$(HISPID_DIR)/%.o: $(SRCD)/%.cpp $(INCD)/HiSpID.h $(INCD)/PunctureKrylov.h $(INCD)/PunctureExecution.h $(SRCD)/HiSpID_jets.hpp $(SRCD)/HiSpID_internal.hpp $(SRCD)/HiSpID_spectral.hpp $(SRCD)/HiSpID_axis.hpp $(HISPID_KERNEL_HEADERS)
	@mkdir -p $(HISPID_DIR)
	$(CXX) $(HISPID_FLAGS) $(INC_PARAMS) -c $< -o $@

$(HISPID_LIB): $(HISPID_OBJ) $(STATIC_LIB)
	$(CXX) -shared $(HISPID_OBJ) $(STATIC_LIB) $(LFLAGS) -o $@

hispid: $(HISPID_LIB)

$(HISPID_DIR)/test_geometry.x: tests/test_hispid_geometry.cpp $(HISPID_LIB)
	$(CXX) $(HISPID_FLAGS) $(INC_PARAMS) $< $(HISPID_OBJ) $(STATIC_LIB) $(LFLAGS) -o $@

$(HISPID_DIR)/test_axis.x: tests/test_hispid_axis.cpp $(SRCD)/HiSpID_axis.hpp $(SRCD)/HiSpID_spectral.hpp
	@mkdir -p $(HISPID_DIR)
	$(CXX) $(HISPID_FLAGS) $(INC_PARAMS) $< -o $@

$(HISPID_DIR)/test_solver.x: tests/test_hispid_solver.cpp $(SRCD)/HiSpID_solver.cpp $(SRCD)/HiSpID_axis.hpp $(SRCD)/HiSpID_spectral.hpp $(HISPID_DIR)/HiSpID_geometry.o $(STATIC_LIB)
	$(CXX) $(HISPID_FLAGS) $(INC_PARAMS) $< $(HISPID_DIR)/HiSpID_geometry.o $(STATIC_LIB) $(LFLAGS) -o $@

$(HISPID_DIR)/test_charge_rings.x: tests/test_hispid_charge_rings.cpp $(SRCD)/HiSpID_solver.cpp $(SRCD)/HiSpID_axis.hpp $(SRCD)/HiSpID_spectral.hpp $(HISPID_DIR)/HiSpID_geometry.o $(STATIC_LIB)
	$(CXX) $(HISPID_FLAGS) $(INC_PARAMS) $< $(HISPID_DIR)/HiSpID_geometry.o $(STATIC_LIB) $(LFLAGS) -o $@

$(HISPID_DIR)/test_preconditioner_reuse.x: tests/test_hispid_preconditioner_reuse.cpp $(SRCD)/HiSpID_solver.cpp $(SRCD)/HiSpID_axis.hpp $(SRCD)/HiSpID_spectral.hpp $(HISPID_DIR)/HiSpID_geometry.o $(STATIC_LIB)
	$(CXX) $(HISPID_FLAGS) $(INC_PARAMS) $< $(HISPID_DIR)/HiSpID_geometry.o $(STATIC_LIB) $(LFLAGS) -o $@

test-hispid-native: $(HISPID_DIR)/test_geometry.x $(HISPID_DIR)/test_axis.x $(HISPID_DIR)/test_solver.x $(HISPID_DIR)/test_charge_rings.x $(HISPID_DIR)/test_preconditioner_reuse.x
	$(HISPID_DIR)/test_geometry.x
	$(HISPID_DIR)/test_axis.x
	$(HISPID_DIR)/test_solver.x
	$(HISPID_DIR)/test_charge_rings.x
	$(HISPID_DIR)/test_preconditioner_reuse.x

PYTHON ?= python3
test-hispid: test-hispid-native
	OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=$(BASE)/python:$(BASE)/validation:$(BASE)/examples HISPID_LIBRARY=$(HISPID_LIB) $(PYTHON) -m unittest discover -s tests -p 'test_*.py' -v

clean-hispid:
	rm -rf $(HISPID_DIR)

.PHONY: hispid test-hispid-native test-hispid clean-hispid

$(HISPID_DIR)/test_krylov_memory.x: tests/test_hispid_krylov_memory.cpp $(SRCD)/HiSpID_solver.cpp $(HISPID_DIR)/HiSpID_geometry.o $(STATIC_LIB)
	$(CXX) $(HISPID_FLAGS) $(INC_PARAMS) $< $(HISPID_DIR)/HiSpID_geometry.o $(STATIC_LIB) $(LFLAGS) -o $@

$(HISPID_DIR)/test_by_line_cache.x: tests/test_by_line_cache.c $(SRCD)/TP_Newton.c $(INC) $(STATIC_LIB)
	$(CC) $(CFLAGS) $(INC_PARAMS) $< $(filter-out $(OBJD)/TP_Newton.o,$(LIBOBJ)) $(LFLAGS) -o $@

test-solver-efficiency: $(HISPID_DIR)/test_by_line_cache.x $(HISPID_DIR)/test_krylov_memory.x
	$(HISPID_DIR)/test_by_line_cache.x
	$(HISPID_DIR)/test_krylov_memory.x

.PHONY: test-solver-efficiency

# Include exact optimization controls in the normal complete suite.
test-hispid: test-solver-efficiency

$(HISPID_DIR)/test_by_modal.x: tests/test_by_modal.c $(STATIC_LIB) $(INC)
	@mkdir -p $(HISPID_DIR)
	$(CC) $(CFLAGS) $(INC_PARAMS) $< $(STATIC_LIB) $(LFLAGS) -o $@

test-by-modal: $(HISPID_DIR)/test_by_modal.x
	$(HISPID_DIR)/test_by_modal.x

test-hispid: test-by-modal
.PHONY: test-by-modal

$(HISPID_DIR)/test_krylov.x: tests/test_krylov.c $(SRCD)/PunctureKrylov.c $(INCD)/PunctureKrylov.h
	@mkdir -p $(HISPID_DIR)
	$(CC) $(CFLAGS) $(INC_PARAMS) tests/test_krylov.c $(SRCD)/PunctureKrylov.c -lm -o $@

test-krylov: $(HISPID_DIR)/test_krylov.x
	$(HISPID_DIR)/test_krylov.x

test-hispid: test-krylov $(LIB)
.PHONY: test-krylov
