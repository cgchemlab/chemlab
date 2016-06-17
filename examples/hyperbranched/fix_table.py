#!/usr/bin/env python

# Fix the tabulated potential by settings force for first and last
# element to the values of next or previous element.

import numpy as np
import sys

d = np.loadtxt(sys.argv[1])
# if first row has force eq 0.0 then take the force value from the next element
if d[0][2] == 0.0:
    d[0][2] = d[1][2]

if d[-1][2] == 0.0:
    d[-1][2] = d[-2][2]

np.savetxt(sys.argv[1], d)
