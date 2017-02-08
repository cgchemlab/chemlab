#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import argparse
import h5py
import numpy as np

parser = argparse.ArgumentParser('Analyze structures in hyperbranched network')
parser.add_argument('h5', help='Input h5md file')
parser.add_argument('out', help='Output file')



def get_structures(h5, time_frame):
    bead_types = h5['/particles/atoms/species/value'][time_frame]
    bead_of_type_ra = len(bead_types[bead_types == 6])
    # Bead types signature
    res_order = [
        bead_types[i]*bead_types[i+1]*bead_types[i+2]*bead_types[i+3] for i in range(0, 4000, 4)
    ]
    res_type2name = {
        0: 'n_1t',    # [0,2,1,1]
        420: 'n_1l',  # [3,5,4,7] or [3,5,7,4]
        735: 'n_1d',  # [3,5,7,7]
        480: 'n_2t',  # [6,5,4,4]
        840: 'n_2l',  # [6,5,4,7] or [6,5,7,4]
        1470: 'n_2d'  # [6,5,7,7]
    }
    out = {k: 0 for k in res_type2name}
    for r in res_order:
        out[r] += 1
    return bead_of_type_ra, {res_type2name[k]: out[k] for k in res_type2name}

args = parser.parse_args()

h5 = h5py.File(args.h5, 'r')


structures_cr = {k: [v/1000.0] for k, v in get_structures(h5, 0)[1].items()}
cr_list = [0.0]
last_cr = 0.1
for time_frame in range(10002):
    if last_cr > 1.0:
        break
    conversion, structures = get_structures(h5, time_frame)
    if (conversion/1000.0) > last_cr:
        cr_list.append((conversion/1000.0))
        for k, v in structures.items():
            structures_cr[k].append(v/1000.0)
        last_cr += 0.1
header = 'cr ' + ' '.join(sorted(structures_cr.keys()))
data_out = []
for i in range(len(cr_list)):
    t = [cr_list[i]]
    for k in sorted(structures_cr.keys()):
        t.append(structures_cr[k][i])
    data_out.append(t)
np.savetxt(args.out, data_out, header=header)
