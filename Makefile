#modied from htslib makefile
FLAGS=-O3

CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)


all: ngsPSMC

# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif

PACKAGE_VERSION  = 0.01

ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

version.h:
	echo '#define ngsPSMC_VERSION "$(PACKAGE_VERSION)"' > $@

.PHONY: misc clean test

-include $(OBJ:.o=.d)

ifdef HTSSRC
%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTS_INCDIR) -std=c99 $*.c
	$(CC) -MM $(CFLAGS)  -I$(HTS_INCDIR) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp >$*.d

ngsPSMC: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o ngsPSMC *.o $(HTS_LIBDIR) -lz -lm -lbz2 -llzma -lpthread
else
%.o: %.c
	$(CC) -c  $(CFLAGS) -std=c99 $*.c
	$(CC) -MM $(CFLAGS)  $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  $*.cpp
	$(CXX) -MM $(CXXFLAGS)  $*.cpp >$*.d

ngsPSMC: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o ngsPSMC *.o -lz -lpthread -lhts  -lbz2 -llzma
endif

clean:
	rm  -f *.o *.d ngsPSMC version.h *~

force:
