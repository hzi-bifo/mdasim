#saigon (SunOS), gcc without openmp
#C = /usr/sfw/bin/g++
#CPP = /usr/sfw/bin/g++
#CFLAGS = -m64 -Wall -pedantic -O3 

#saigon (SunOS), cc with openmp
#C = cc
#CPP = CC
#CFLAGS = -fast -xarch=amd64 -xopenmp 

#IRIX, cc with openmp
#C = cc
#CPP = CC
#CFLAGS = -mp -O3 -64 -woffall

#C = gcc
#CPP = g++
#CFLAGS = -maix64 -Wall -O3 

#AIX, xlc_r with openmp
#C = xlc_r
#CPP = xlC_r
#CFLAGS = -q64 -O5 -qsmp=omp 

#AIX, xlc_r with openmp
#C = mpcc_r
#CPP = mpCC_r
#CPPFLAGS = -q64 -O5 -qsmp=omp -DMAXK=$(MAXK) -DFIXEDPOINTPROB=1
#LIBEXT = a
#LIBCREATE = -rcs -X64
#LIBTOOL = ar

#Linux, icc with openmp & pthread
#CPP = g++
CPP = mpic++
#CPPFLAGS = -m64 -fopenmp -O3 -s -static
MPICPPFLAGS = -pthread -lmpi_cxx -lmpi -lopen-rte -lopen-pal -ldl -lnsl -lutil -lm -ldl -Wl,--export-dynamic
CPPFLAGS = -m64 -openmp -pthread -O3 -DMAXK=$(MAXK) -DFIXEDPOINTPROB=1 $(MPICPPFLAGS)

#Linux, gcc with openmp
#CPP = g++
#CPPFLAGS = -m64 -fopenmp -O3 -s -static
#PROFILINGFLAG = -pg
#VALGRINDFLAG = -g
#DEBUGFLAG = -g
#CPPFLAGS = -m64 -fopenmp -O3 $(PROFILINGFLAG) $(DEBUGFLAG) -s
LIBEXT = a
LIBCREATE = rcs
LIBTOOL = ar 

PACKAGENAME = mdasim

TARGET = $(prefix)

OBJDIR = $(TARGET)/obj
BINDIR = $(TARGET)/bin

BASE = ./
SRCDIR = $(BASE)/src
INCLUDEDIR = $(BASE)/include

.PHONY: clean $(PACKAGENAME) lib

all : $(BINDIR)/$(PACKAGENAME)

$(OBJDIR): 
	mkdir -p $(OBJDIR)
$(BINDIR): 
	mkdir -p $(BINDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.C | $(OBJDIR)
	$(CPP) $(CPPFLAGS) -I$(INCLUDEDIR) -c $< -o $@ 

$(BINDIR)/%: $(OBJDIR)/%.o | $(BINDIR)
	$(CPP) $(CPPFLAGS) $< -l$(PACKAGENAME) -o $@

clean :
	rm -f obj/*.o *~ src/*~ include/*~ bin/* lib/*
	rm -f -r obj bin lib





