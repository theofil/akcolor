# Build:  make
#   needs ROOT and Python
#   conda install -c conda-forge root fastjet

CXX = g++

FASTJET_CONFIG := $(shell which fastjet-config 2>/dev/null || \
    echo /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fastjet/3.4.1-5e99cb03f21d718e60d5b9b541840b9a/bin/fastjet-config)

CXXFLAGS = $(shell root-config --cflags) $(shell $(FASTJET_CONFIG) --cxxflags) -O2 -std=c++17
LDFLAGS  = $(shell root-config --libs)   $(shell $(FASTJET_CONFIG) --libs)

makefriends: makefriends.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

# nov26 campaign fork (see nov26.md); build with `make makefriends_nov26`
makefriends_nov26: makefriends_nov26.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

# fall26 campaign fork (see fall26.md); build with `make makefriends_fall26`
makefriends_fall26: makefriends_fall26.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f makefriends makefriends_nov26 makefriends_fall26
