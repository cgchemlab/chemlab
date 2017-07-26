#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import argparse
import datetime
import numpy as np

parser = argparse.ArgumentParser('Mix table')
parser.add_argument('--tab1')
parser.add_argument('--tab2')
parser.add_argument('--output')
parser.add_argument('--scaling', type=float, default=0.5)
parser.add_argument('--constant', type=float, default=0.0)
parser.add_argument('--mix_type', type=int, default=0, choices=[0, 1],
                    help='coupling type, 0 for arithmetic, 1 for geometric')

args = parser.parse_args()


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

mono_tab = convertGromacsESPP(np.loadtxt(args.tab1))
poly_tab = convertGromacsESPP(np.loadtxt(args.tab2))
out_name = args.output

if args.mix_type == 0:
    mixed_table = mix_arithmetic(mono_tab, poly_tab, args.scaling)
elif args.mix_type == 1:
    mixed_table = mix_geometric(mono_tab, poly_tab, args.scaling, args.constant)
else:
    raise RuntimeError('Wrong mix type')
print('Saved {}'.format(out_name))
np.savetxt(
    out_name,
    mixed_table,
    header='Mixed of {} and {} at {}, method={} scaling={}'.format(
	args.tab1, args.tab2, datetime.datetime.now(), args.mix_type, args.scaling),
    fmt='%2.9e')
