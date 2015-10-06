PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

#CXX=icpc
CXXFLAGS=-Wall -g -O0 -std=c++11
CPPFLAGS=-I$(PREFIX)/include
LDLIBS=-L$(PREFIX)/lib -larmadillo -lhdf5 -lm

PROGRAMS=subsample_seer

all: $(PROGRAMS)

clean:
	$(RM) *.o ~* $(PROGRAMS)

install: all
	install -d $(BINDIR)
	install $(PROGRAMS) $(BINDIR)

subsample_seer: subsample_seer.o
	$(LINK.cpp) $^ $(LDLIBS) -o $@

.PHONY: all test clean install

