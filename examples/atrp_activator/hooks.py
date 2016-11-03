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


def hook_init_reaction(system, integrator, topol, args):
    name2type = topol.atomsym_atomtype
    number_to_activate = 10

    max_pid = espressopp.analysis.MaxPID(system).compute()
    pid = 0
    activated = 0
    last_res_id = -1
    activated_monomer = False
    while pid <= max_pid and activated <= number_to_activate:
        if system.storage.particleExists(pid):
            p = system.storage.getParticle(pid)
            if p.type == name2type['MA']:
                if not activated_monomer:
                    new_property = topol.gt.atomtypes['FA']
                    # Transfer type to FA
                    p.type = name2type['FA']
                    p.state = 3
                    p.mass = new_property['mass']
                    activated_monomer = True
                else:
                    new_property = topol.gt.atomtypes['PA']
                    # Transfer type to PA
                    p.type = name2type['PA']
                    p.mass = new_property['mass']
            elif p.type == name2type['ML']:
                new_property = topol.gt.atomtypes['PL']
                p.type = name2type['PL']
                p.mass = new_property['mass']
            if last_res_id != p.res_id():
                last_res_id = p.res_id
                activated += 1
                activated_monomer = False
        pid += 1

    return activated == number_to_activate