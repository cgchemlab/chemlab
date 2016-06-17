#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import argparse
import h5py
import numpy as np
import networkx as nx

parser = argparse.ArgumentParser('Analyze structures in hyperbranched network')
parser.add_argument('h5', help='Input h5md file')
parser.add_argument('out_postfix', help='Output file postfix')

def get_structures(h5, time_frame):
    bead_types = h5['/particles/atoms/species/value'][time_frame]
    name2type_id = {t[1]: int(t[0]) for t in h5['/parameters/force_field/atomtypes']}
    bead_of_type_ra = len(bead_types[bead_types == name2type_id['RA']])
    # Bead types signature
    res_order = [
        bead_types[i]*bead_types[i+1]*bead_types[i+2]*bead_types[i+3] for i in range(0, 4000, 4)
    ]
    res_type2name = {
        ('MA', 'ML', 'MB', 'MB'): 'n_1t',    # [0,2,1,1]
        ('PA', 'PL', 'PB', 'RB'): 'n_1l',  # [3,5,4,7] or [3,5,7,4]
        ('PA', 'PL', 'RB', 'RB'): 'n_1d',  # [3,5,7,7]
        ('RA', 'PL', 'PB', 'PB'): 'n_2t',  # [6,5,4,4]
        ('RA', 'PL', 'RB', 'PB'): 'n_2l',  # [6,5,4,7] or [6,5,7,4]
        ('RA', 'PL', 'RB', 'RB'): 'n_2d'  # [6,5,7,7]
    }
    res_type2name = {reduce(lambda x, y: x*y, map(name2type_id.get, k)): v for k, v in res_type2name.items()}
    out = {k: 0 for k in res_type2name}
    for r in res_order:
        out[r] += 1
    return bead_of_type_ra, {res_type2name[k]: out[k] for k in res_type2name}


def CalPolym(h5, time_frame):
    cg_connections = set()
    cg_activated = set()
    bond_count = 0
    cl = h5['/connectivity/chem_bonds_0/value']
    size_of_cg = 4
    cg_beads = {
        x: range(x*size_of_cg+1, x*size_of_cg+size_of_cg+1)
        for x in range(1000)
    }
    at_cg_bead = {
        a: k for k, v in cg_beads.items() for a in v
    }
    
    for b1, b2 in cl[time_frame]:
        if b1 != -1 and b2 != -1:
            cg1, cg2 = map(at_cg_bead.get, (b1, b2))
            cg_connections.add(tuple(sorted((cg1, cg2))))
            cg_activated.add(cg1)
            cg_activated.add(cg2)
            bond_count += 1
    conversion = bond_count/1000.0
    g = nx.Graph()
    g.add_edges_from(cg_connections)
    graphs = list(nx.connected_component_subgraphs(g))
    return conversion, graphs



def main():
    args = parser.parse_args()

    h5 = h5py.File(args.h5, 'r')

    cl = h5['/connectivity/chem_bonds_0/value']
    
    convA = [0.0]
    Pw = [0.0]
    Pn = [0.0]
    PDI = [0.0]
    NumChain = [0]
    increament = 0.1
    for time in range(0, len(cl)):
        if increament > 1.0:
            break
        conv, add = CalPolym(h5, time)
        if conv > increament:
            Pw_temp = 0.0
            Pn_temp = 0.0
            convA.append(conv)
            NumChain.append(len(add))
            for i in range(0, len(add)):
                Pw_temp += pow(len(add[i]),2)
                Pn_temp += len(add[i])
            Pw_temp = Pw_temp/Pn_temp
            Pn_temp = Pn_temp/len(add)
            Pw.append(Pw_temp)
            Pn.append(Pn_temp)
            PDI.append(Pw_temp/Pn_temp)
            increament += 0.1
    
    out_file = '{}_Pw_Pn_PDI.csv'.format(args.out_postfix)
    data = np.column_stack((convA, Pw, Pn, PDI))
    np.savetxt(out_file, data, header='cr Pw Pn PDI')
    print('Writen file {}, size {}'.format(out_file, data.shape))
    
    frac = {k: [v/1000.0] for k, v in get_structures(h5, 0)[1].items()}
    cr_list = [0.0]
    last_cr = 0.1
    for time_frame in range(0, len(cl)):
        if last_cr > 1.0:
            break
        bond, structures = get_structures(h5, time_frame)
        if (bond/1000.0) > last_cr:
            cr_list.append((bond/1000.0))
            for k, v in structures.items():
                frac[k].append(v/1000.0)
            last_cr += 0.1
    header = 'cr ' + ' '.join(sorted(frac.keys()))
    data_out = []
    for i in range(len(cr_list)):
        t = [cr_list[i]]
        for k in sorted(frac.keys()):
            t.append(frac[k][i])
        data_out.append(t)
    np.savetxt('{}_struct.csv'.format(args.out_postfix), data_out, header=header)
    print('Writen file struct_{}.csv, size {}'.format(args.out_postfix, np.array(data_out).shape))

if __name__ == '__main__':
    main()
