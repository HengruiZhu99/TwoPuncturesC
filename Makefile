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
