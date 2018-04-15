#!/usr/bin/env python
#  Copyright (C) 2012,2013,2015(H),2016
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


import argparse
import math
import re


def convertTable(gro_in_file, esp_out_file, sigma=1.0, epsilon=1.0, c6=1.0, c12=1.0):
    """Convert GROMACS tabulated file into ESPResSo++ tabulated file (new file
    is created). First column of input file can be either distance or angle.
    For non-bonded files, c6 and c12 can be provided. Default value for sigma, epsilon,
    c6 and c12 is 1.0. Electrostatics are not taken into account (f and fd columns).

    Keyword arguments:
    gro_in_file -- the GROMACS tabulated file name (bonded, nonbonded, angle
    or dihedral).
    esp_out_file -- filename of the ESPResSo++ tabulated file to be written.
    sigma -- optional, depending on whether you want to convert units or not.
    epsilon -- optional, depending on whether you want to convert units or not.
    c6 -- optional
    c12 -- optional
    """
    # determine file type
    bonded, angle, dihedral = False, False, False

    re_bond = re.compile('.*_b[0-9]+.*')
    re_angle = re.compile('.*_a[0-9]+.*')
    re_dihedral = re.compile('.*_d[0-9]+.*')

    if re.match(re_bond, gro_in_file):
        bonded = True
    elif re.match(re_angle, gro_in_file):
        angle  = True
        bonded = True
    elif re.match(re_dihedral, gro_in_file):
        dihedral = True
        bonded = True

    fin = open(gro_in_file, 'r')
    fout = open(esp_out_file, 'w')

    if bonded: # bonded has 3 columns
        for line in fin:
            if line[0] == "#": # skip comment lines
                continue

            columns = line.split()
            r = float(columns[0])
            f = float(columns[1]) # energy
            fd= float(columns[2]) # force

            # convert units
            if angle or dihedral: # degrees to radians
                r = math.radians(r)
                fd=fd*180/math.pi
            else:
                r = r / sigma
            e = f / epsilon
            f = fd*sigma / epsilon

            if (not angle and not dihedral and r != 0) or \
                 (angle and r <= math.pi and r > 0) or \
                  (dihedral and r >= -math.pi and r <= math.pi):
                fout.write("%15.8g %15.8g %15.8g\n" % (r, e, f))

    else: # non-bonded has 7 columns
        for line in fin:
            if line.startswith('#'): # skip comment lines
                continue

            columns = line.split()
            r = float(columns[0])
            g = float(columns[3]) # dispersion
            gd= float(columns[4])
            h = float(columns[5]) # repulsion
            hd= float(columns[6])

            e = c6*g + c12*h
            f = c6*gd+ c12*hd

            # convert units
            r = r / sigma
            e = e / epsilon
            f = f*sigma / epsilon

            if r != 0: # skip 0
                fout.write("%15.8g %15.8g %15.8g\n" % (r, e, f))

    fin.close()
    fout.close()

def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_file')
    parser.add_argument('out_file')

    return parser

def main():
    args = _args().parse_args()

    convertTable(args.in_file, args.out_file)


if __name__ == '__main__':
    main()
