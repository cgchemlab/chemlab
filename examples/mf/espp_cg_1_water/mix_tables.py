#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import argparse
import datetime
import numpy as np


def _args():
    parser = argparse.ArgumentParser('Mix table')
    parser.add_argument('--table1', '--tab1', dest='tab1', required=True, help='Table U1')
    parser.add_argument('--table2', '--tab2', dest='tab2', required=True, help='Table U2')
    parser.add_argument('--out_table', required=True, help='Output table')
    parser.add_argument('--scaling', type=float, default=0.5)
    parser.add_argument('--constant', type=float, default=0.0)
    parser.add_argument('--mix_type', type=int, default=0, choices=[0, 1],
                        help='coupling type, 0 for arithmetic, 1 for geometric')
    parser.add_argument('--output_type', default='espp', choices=['espp', 'gromacs'],
                        help='Define the output format')

    args = parser.parse_args()
    return args


def convertGromacsESPP(xvg):
    """Converts GROMACS .xvg table to ESPP r e f format."""
    eps = sig = c6 = c12 = 1.0
    # format r e f
    output_tab = np.zeros((xvg.shape[0], 3))
    r = xvg[:, 0]  # r column
    g = xvg[:, 3]
    dg = xvg[:, 4]
    h = xvg[:, 5]
    dh = xvg[:, 6]

    e = c6*g + c12*h
    f = c6*dg + c12*dh

    r /= sig
    e /= eps
    f = f*sig/eps

    output_tab[:, 0] = r
    output_tab[:, 1] = e
    output_tab[:, 2] = f
    return output_tab


def convertEsppToGromacs(espp):
    output_table = np.zeros((espp.shape[0], 7))

    output_table[:, 0] = espp[:, 0]
    output_table[:, 5] = espp[:, 1]
    output_table[:, 6] = espp[:, 2]

    return output_table


def mix_arithmetic(tab1, tab2, coupling):
    max_length = 0
    if tab1.shape[0] != tab2.shape[0]:
        if tab1.shape[0] <= tab2.shape[0]:
            out_tab = np.array(tab1)
        else:
            out_tab = np.array(tab2)
        if (tab1[:, 0][:out_tab.shape[0]] != tab2[:, 0][:out_tab.shape[0]]).all():
            raise RuntimeError('Both r columns should be the same')
    else:
        out_tab = np.array(tab1)
    max_length = out_tab.shape[0]
    if max_length == 0:
        raise RuntimeError('The length of output table is zero???')
    out_tab[:, 1] = coupling*tab1[:max_length, 1] + (1.0-coupling)*tab2[:max_length, 1]
    out_tab[:, 2] = coupling*tab1[:max_length, 2] + (1.0-coupling)*tab2[:max_length, 2]
    return out_tab


def mix_geometric(tab1, tab2, coupling, constant):
    max_length = 0
    if tab1.shape[0] != tab2.shape[0]:
        if tab1.shape[0] <= tab2.shape[0]:
            out_tab = np.array(tab1)
        else:
            out_tab = np.array(tab2)
        if (tab1[:, 0][:out_tab.shape[0]] != tab2[:, 0][:out_tab.shape[0]]).all():
            raise RuntimeError('Both r columns should be the same')
    else:
        out_tab = np.array(tab1)
    max_length = out_tab.shape[0]
    if max_length == 0:
        raise RuntimeError('The length of output table is zero???')
    e1 = tab1[:max_length, 1]
    e2 = tab2[:max_length, 1]
    f1 = tab1[:max_length, 2]
    f2 = tab2[:max_length, 2]
    out_tab[:, 1] = (pow(e1+constant, coupling) +
                     pow(e2+constant, (1-coupling)) - constant)
    out_tab[:, 2] = (coupling*pow(e1 + constant, coupling-1.0)*f1 +
                     (1.0-coupling)*pow(e2+constant, -coupling)*f2)
    return out_tab

def main():
    args = _args()

    if args.tab1.endswith('xvg'):
        u1 = convertGromacsESPP(np.loadtxt(args.tab1))
    else:
        u1 = np.loadtxt(args.tab1)
    if args.tab2.endswith('xvg'):
        u2 = convertGromacsESPP(np.loadtxt(args.tab2))
    else:
        u2 = np.loadtxt(args.tab2)

    if args.mix_type == 0:
        print('Performing arthmetic mixing')
        print('U_final = x*U1 + (1-x)*U2 where:')
        print('\tx={}\n\tU1={}\n\tU2={}'.format(args.scaling, args.tab1, args.tab2))
        mixed_table = mix_arithmetic(u1, u2, args.scaling)
    elif args.mix_type == 1:
        print('Performing geometric mixing')
        print('U_final = U1^x + U2^(1-x) - C where:')
        print('\tx={}\n\tU1={}\n\tU2={}\n\tC={}'.format(
            args.scaling, args.tab1, args.tab2, args.constant))
        mixed_table = mix_geometric(u1, u2, args.scaling, args.constant)
    else:
        raise RuntimeError('Unknown mixing type')

    if args.output_type == 'gromacs':
        mixed_table = convertEsppToGromacs(mixed_table)

    print('Saved {}'.format(args.out_table))
    np.savetxt(
        args.out_table,
        mixed_table,
        header='Mixed of {} and {} at {}'.format(
            args.tab1, args.tab2, datetime.datetime.now()),
        fmt='%2.9e')

if __name__ == '__main__':
    main()
