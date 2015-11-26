PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

CXX=icpc
CXXFLAGS=-Wall -O3 -parallel -ipo -std=c++11
CPPFLAGS=-I$(PREFIX)/include
LDLIBS=-L$(PREFIX)/lib -lboost_filesystem -lboost_system -lhdf5 -larmadillo -lz -lm

PROGRAMS=subsample_seer trimeth_subsample

all: $(PROGRAMS)

clean:
	$(RM) *.o ~* $(PROGRAMS)

install: all
	install -d $(BINDIR)
	install $(PROGRAMS) $(BINDIR)

.PHONY: all test clean install

