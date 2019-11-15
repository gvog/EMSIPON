#/// @file
#/// @author Georgios G. Vogiatzis
#/// @version 2.0 
#/// @brief A simple Makefile for compiling the source code. If you have Intel compilers installed
#///        it will use Intel compilers. Otherwise the GNU compilers will be used instead.

#/// @section LICENSE
#/// Copyright (c) 2012 Georgios Vogiatzis. 
#/// All rights reserved.
#/// This work is licensed under the terms of the MIT license.  
#/// For a copy, see <https://opensource.org/licenses/MIT>.


.SUFFIXES:  (.SUFFIXES) .c .cpp
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
## PREPROCESSOR                                ##
## Linux:
CPP = /lib/cpp -P -C -traditional -t -W
CPPFLAGS = #-DUSE_SSE
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

OBJDIR=obj
SRCDIR=src
INCLDIR:=-Isrc/include 


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# SHARED SOURCE
# Shared source files directories:
SHARDIR = src/include src
vpath %.h $(SHARDIR)

# LIBRARIES
LIBS = #-lstdc++

# LINKER
LD = ld


ifeq  ($(notdir $(shell which icc 2>&1)),icc)
CCOMP   = icc
CPPCOMP = icpc
RFLAGS = -O3 -xHost -prec-div -prec-sqrt -align -static -ip -ipo
RFLAGS = -g -O2 -fno-inline-functions
DFLAGS = -O0 -g -traceback -debug all
INCLDIR := $(INCLDIR) -I/usr/include/x86_64-linux-gnu/c++/4.8
else
CCOMP   = gcc 
CPPCOMP = g++ 
RFLAGS = -O3 -m64 -mtune=native -march=native -Wall
DFLAGS = -O0 -g -Wall 
endif

#CCFLAGS = $(DFLAGS)
CCFLAGS = $(RFLAGS) 

LDFLAGS = $(CCFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CPPCOMP) $(CCFLAGS) $(CPPFLAGS) $(INCLDIR) -c -o $@ $?

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CCOMP) $(CCFLAGS) $(INCLDIR) -c -o $@ $?

COMMON_OBJS =  $(OBJDIR)/rng.o\
               $(OBJDIR)/network.o\
               $(OBJDIR)/domain.o\
	       $(OBJDIR)/Auxiliary.o\
               $(OBJDIR)/b3D_integrator.o\
	       $(OBJDIR)/non_bonded_scheme_routines.o\
	       $(OBJDIR)/distributions.o\
	       $(OBJDIR)/grid.o\
	       $(OBJDIR)/hopping.o\
	       $(OBJDIR)/dump.o\
	       $(OBJDIR)/netmin.o

MAIN_OBJS = $(OBJDIR)/main.o

OBJECTS  = $(COMMON_OBJS) $(MAIN_OBJS)

# build
build: $(OBJECTS)
	$(CPPCOMP) $(LDFLAGS) $(LIBS) -o netmin.x $(OBJECTS)
# clean
clean :
	rm -f $(OBJDIR)/*.o $(PARSE_C) $(PARSE_H)
	rm -f */*~ *~ core *.out netmin.x
	rm -f dump_b3D.lammpstrj
	rm -f restart.data ss_lifetimes.txt events.txt
	cd ./doxygen/latex; rm -f *

# documentation
doc: doxygen/Doxyfile
	cd ./doxygen/; doxygen; cd ./latex; make; cp refman.pdf ../../manual.pdf


