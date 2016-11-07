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

# Defines hooks run in different part of the simulation toolbox.
# This is really advance option, use it if you know what are you doing.
# Every function will get in the paramters:
#  - system, integrator, args objects.
# It should return status at the end, True if success or False if failure.
# on failure, the simulation will be stoped.

# This file do not need to be system independent.


import espressopp


def hook_postsetup_reaction(system, integrator, topol, args, chemical_reaction):
    name2type = topol.atomsym_atomtype