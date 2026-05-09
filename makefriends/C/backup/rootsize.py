#!/usr/bin/env python3
import sys
import ROOT

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <file.root>")
    sys.exit(1)

f = ROOT.TFile(sys.argv[1])
if not f or f.IsZombie():
    print(f"Error: cannot open {sys.argv[1]}")
    sys.exit(1)

for k in f.GetListOfKeys():
    print(f"{k.GetNbytes()/1e6:.3f} MB  (uncompressed: {k.GetObjlen()/1e6:.3f} MB)  {k.GetName()}")

