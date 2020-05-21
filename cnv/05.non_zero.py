#!/usr/bin/env python

import sys
import h5py
import numpy as np

f = h5py.File(sys.argv[1], "r")
values = f['counts']['values'][0]

print(f"File:  {sys.argv[1]}")
print(f"Nonzero: {np.count_nonzero(values)}")
print(f"Total: {len(values)}")
