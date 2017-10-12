#!/usr/bin/env python
"""
Copyright (C) 2017 Jakub Krajniak <jkrajniak@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('out')

    return parser.parse_args()


def main():
    args = _args()
    cg_bonds = []
    ter = {k: (i, i+1, i+2) for k, i in enumerate(range(1001, 4001, 3), 1001)}

    for n in range(1, 1000, 3):
        cg_bonds.extend([
            (n, ter[1000+n][0]),
            (ter[1000+n][2], n+1),
            (n+1, ter[1000+n+1][0]),
            (ter[1000+n+1][2], n+2),
            (n+2, ter[1000+n+2][0]),
            (n, ter[1000+n+2][2])])

    with open(args.out, 'w') as outf:
        for l in cg_bonds:
            outf.write('{} {}\n'.format(*l))


if __name__ == '__main__':
    main()

