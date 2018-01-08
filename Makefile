#
# make         # compile library
# make forlanl # saves src files for LANL in private/src_forlanl (svn export first)
# make test    # run regression tests
#
# default compilers and flags are defined in src/CMakeLists.txt
# see src/Readme for expert build instructions
#
# Doug Wright, LLNL 8-24-2014

# default location to do build
BUILD_DIR = build
SRC_DIR = src

default: $(BUILD_DIR)
	cd $(BUILD_DIR); CXX=`which g++` CC=`which gcc` cmake ../src
	$(MAKE) -C $(BUILD_DIR) install

forlanl:
	cd $(SRC_DIR); CXX=`which g++` CC=`which gcc` cmake -DCOPYRIGHT:BOOL=ON -DFORLANL:BOOL=ON -DNEWXNORMAL:BOOL=ON .

#....create build directory
$(BUILD_DIR):
	mkdir $(BUILD_DIR)

#....run regression tests
test: default
	cd regr; ./run_test

#....commands handled by this make, should be filtered out from the following
FILTER_OUT = test forlanl default
#....pass command line targets to default build makefiles
MAKECMDGOALS := $(filter-out $(FILTER_OUT),$(MAKECMDGOALS))
$(MAKECMDGOALS):
	$(MAKE) -C $(BUILD_DIR) $(MAKECMDGOALS)

.PHONY: default forlanl test

MAKEFLAGS += --no-print-directory

# Notes
# -----
# March 23, 2017
# --------------
# By default, cmake searches for a list of known compiler names:
#   cc, gcc, cl, bcc, xlc, clang
# in that order, in PATH. On LC (borax), 
#   cc points to gcc 4.8.5 (chosen by default)
#   gcc points to gcc 4.9.3
# Because of incompatibilities between gcc/g++ version 4.8.5 and versions >=4.9
# leading to link time errors, we need to set the gcc/g++ compiler versions
# to >=4.9. To do this, prepend the cmake line by CXX=`which g++` CC=`which gcc`:
#   CXX=`which g++` CC=`which gcc` cmake ../src; make
