#!/usr/bin/env python
#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
#
#  This file is part of ChemLab.
#
#  ChemLab is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ChemLab is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import sys

d = np.loadtxt(sys.argv[1])
# if first row has force eq 0.0 then take the force value from the next element
if d[0][2] == 0.0:
    d[0][2] = d[1][2]

if d[-1][2] == 0.0:
    d[-1][2] = d[-2][2]

np.savetxt(sys.argv[1], d)
