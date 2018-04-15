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


import collections
import espressopp
import random
import numpy

random.seed(None)

class AngleDistribution:
    def __init__(self):
        self.histogram = numpy.zeros(100)
        self.steps = 0

angle_dist = AngleDistribution()

def hook_before_sim(system, integrator, ar, gt):
    global angle_dist
    analysis_angles = espressopp.analysis.AngleDistribution(system)
    analysis_angles.load_from_topology_manager(system.topology_manager)
    angle_dist.analysis_angles = analysis_angles

def hook_at_step(system, integrator, ar, gt, args, step):
    global angle_dist
    histogram = numpy.array(angle_dist.analysis_angles.compute(100))
    angle_dist.histogram += histogram
    angle_dist.steps += 1

def hook_end(system, integrator, ar, gt, args):
    global angle_dist
    theta = numpy.arange(0.0, numpy.pi, numpy.pi/100)
    numpy.savetxt('output_angle.csv', numpy.column_stack((theta, angle_dist.histogram)))

def hook_init_reaction(system, integrator, ar, topol, args):
    name2type = topol.atomsym_atomtype
    number_to_activate = 20

    max_pid = espressopp.analysis.MaxPID(system).compute()
    pid = 0
    activated = 0
    last_res_id = -1
    activated_monomer = False
    res_ids = range(1, 2001)
    res_id2pids = [(i, i+1) for i in range(1, 4001, 2)]
    res_id2pids = {i: v for i, v in enumerate(res_id2pids, 1)}
    print('Get {} res_ids'.format(len(res_ids)))
    res_ids = random.sample(res_ids, number_to_activate)
    print('Selected {} res_ids'.format(len(res_ids)))

    for res_id in res_ids:
        activated_monomer = False
        for pid in res_id2pids[res_id]:
            p = system.storage.getParticle(pid)
            if p.type == name2type['MA']:
                if not activated_monomer:
                    new_property = topol.gt.atomtypes['FA']
                    system.storage.modifyParticle(pid, 'type', name2type['FA'])
                    system.storage.modifyParticle(pid, 'state', 2)
                    system.storage.modifyParticle(pid, 'mass', new_property['mass'])
                    activated_monomer = True
                else:
                    new_property = topol.gt.atomtypes['PA']
                    p.type = name2type['PA']
                    p.mass = new_property['mass']
                    system.storage.modifyParticle(pid, 'type', p.type)
                    system.storage.modifyParticle(pid, 'mass', p.mass)
            else:
                print('Ignore {} with type {}'.format(pid, p.type))
    system.storage.decompose()
    print ('Activated {} monomers'.format(number_to_activate))

    # Double check
    type_list = []
    for res_id in res_ids:
        for pid in res_id2pids[res_id]:
            p = system.storage.getParticle(pid)
            type_list.append(p.type)
    print res_ids
    assert type_list.count(name2type['FA']) == number_to_activate
    assert type_list.count(name2type['PA']) == number_to_activate
    return True
