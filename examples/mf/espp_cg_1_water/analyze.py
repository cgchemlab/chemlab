#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
# 

import h5py
import sys
import networkx as nx
import numpy as np


h5 = h5py.File(sys.argv[1], 'r')

file_prefix = sys.argv[1].replace('.h5', '')

# [u'bonds_0', u'chem_bonds_0', u'chem_bonds_1']

cr = np.zeros(shape=(h5['/connectivity/chem_bonds_0/time'].shape[0], 2))
cr[:, 0] = h5['/connectivity/chem_bonds_0/time']

cr[:, 1] = [len([b for b in v if -1 not in b]) for v in h5['/connectivity/chem_bonds_0/value']]

np.savetxt('cr_{}.csv'.format(file_prefix), cr)

sys.exit(0)

components = np.zeros(shape=cr.shape)
components[:, 0] = cr[:, 0]

cv0 = h5['/connectivity/chem_bonds_0']
cv1 = h5['/connectivity/dynamic_bonds_1']
for t in range(components.shape[0]):
    g = nx.Graph()
    g.add_edges_from([b for b in cv0['value'][t] if -1 not in b])
    g.add_edges_from([b for b in cv1['value'][t] if -1 not in b])
    components[:, 1][t] = nx.number_connected_components(g)

np.savetxt('cmp_{}.csv'.format(file_prefix), components)
