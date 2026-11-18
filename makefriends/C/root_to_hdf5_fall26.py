#!/usr/bin/env python3
"""Convert fall26 friend trees (makefriends_fall26) to HDF5.

    python3 root_to_hdf5_fall26.py friends/fall26/VBFZ_herwig.friend.root [...]
    python3 root_to_hdf5_fall26.py friends/fall26      # every *.friend.root in the folder

The nov26 converter (root_to_hdf5_nov26.py) with the three fall26 event branches added:
    (N,)  chunk, chunkEntry (int32): input chunk and entry in it; bosonCharge (int8): W charge
Everything else (2 jet slots, 80 constituents, gzip 4, an existing h5 is deleted first) is the
nov26 layout.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import root_to_hdf5_nov26 as base  # noqa: E402

FALL26_VARS = {'chunk': np.int32, 'chunkEntry': np.int32, 'bosonCharge': np.int8}
base.EVT_VARS.update(FALL26_VARS)
base.BRANCHES += list(FALL26_VARS)

if __name__ == '__main__':
    if not sys.argv[1:]:
        sys.exit(__doc__)
    base.main(sys.argv[1:])
