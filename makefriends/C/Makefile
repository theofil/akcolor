# Build:  make
# No cmsenv needed — fastjet is resolved from CVMFS automatically.
#
# To decouple from CMSSW entirely, install via conda-forge instead:
#   conda install -c conda-forge root fastjet
# then open a fresh shell and make works as-is (root-config and
# fastjet-config will come from the conda environment).

CXX = g++

# Use fastjet-config from PATH if available (e.g. after cmsenv),
# otherwise fall back to the known CVMFS location.
FASTJET_CONFIG := $(shell which fastjet-config 2>/dev/null || \
    echo /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fastjet/3.4.1-5e99cb03f21d718e60d5b9b541840b9a/bin/fastjet-config)

CXXFLAGS = $(shell root-config --cflags) $(shell $(FASTJET_CONFIG) --cxxflags) -O2 -std=c++17
LDFLAGS  = $(shell root-config --libs)   $(shell $(FASTJET_CONFIG) --libs)

makefriends: makefriends.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f makefriends
