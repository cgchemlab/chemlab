#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import os
import shutil
import sys

table_pattern = sys.argv[1]
table_inc = int(sys.argv[2])

print('Table pattern: {}, table increment: {}'.format(table_pattern, table_inc))

table_files = sorted([x for x in os.listdir('.') if table_pattern in x])
table_files.reverse()

for f in table_files:
    table_type, table_idx = f.replace('.xvg', '').split('_')[1]
    src_file = f
    dst_file = 'table_{}{}.xvg'.format(table_type, int(table_idx) + table_inc)
    print('Move file {} -> {}'.format(src_file, dst_file))
    shutil.move(src_file, dst_file)
