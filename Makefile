CPP = mpic++

MPICPPFLAGS = -pthread -lmpi_cxx -lmpi -lopen-rte -lopen-pal -ldl -lnsl -lutil -lm -ldl -Wl,--export-dynamic
CPPFLAGS = -m64 -openmp -pthread -O3 -DMAXK=$(MAXK) -DFIXEDPOINTPROB=1 $(MPICPPFLAGS)

PACKAGENAME = mdasim

prefix ?= .
TARGET = $(prefix)

BINDIR = $(TARGET)/bin

BASE = .
OBJDIR = $(BASE)/obj
SRCDIR = $(BASE)/src
INCLUDEDIR = $(BASE)/include

.PHONY: clean $(PACKAGENAME)

all : $(BINDIR)/$(PACKAGENAME)

$(OBJDIR):
	mkdir -p $(OBJDIR)
$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.C | $(OBJDIR)
	$(CPP) $(CPPFLAGS) -I $(INCLUDEDIR) -c $< -o $@

$(BINDIR)/%: $(OBJDIR)/%.o | $(BINDIR)
	$(CPP) $(CPPFLAGS) $< -o $@

clean :
	rm -f obj/*.o *~ src/*~ include/*~ bin/*
	rm -f -r obj bin
