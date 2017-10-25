#  Copyright (C) 2017
#      Jakub Krajniak (jkrajniak at gmail.com)
#  Copyright (C) 2012,2013
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


import espressopp
import mpi4py.MPI as MPI
import sys


def PolymerMelt(num_chains, nsolvents, monomers_per_chain, box=(0, 0, 0), bondlen=0.97, rc=1.12246, skin=0.3,
                dt=0.005, epsilon=1.0, sigma=1.0, shift='auto', temperature=None, xyzfilename=None,
                monomer_type=0, solvent_type=1):
    if xyzfilename:
        pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf = espressopp.tools.readxyz(
            xyzfilename)
        box = (Lxf, Lyf, Lzf)
    else:
        if box[0] <= 0 or box[1] <= 0 or box[2] <= 0:
            print "WARNING: no valid box size specified, box size set to (100,100,100) !"
            box = (100, 100, 100)

    system = espressopp.System()
    system.rng = espressopp.esutil.RNG()
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.skin = skin
    nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
    system.storage = espressopp.storage.DomainDecomposition(
        system, nodeGrid, cellGrid)

    # LJ interaction
    vl = espressopp.VerletList(system, cutoff=rc)
    interaction = espressopp.interaction.VerletListLennardJones(vl)
    interaction.setPotential(
        type1=monomer_type,
        type2=monomer_type,
        potential=espressopp.interaction.LennardJones(epsilon, sigma, rc, shift))
    interaction.setPotential(
        type1=solvent_type,
        type2=solvent_type,
        potential=espressopp.interaction.LennardJones(epsilon, sigma, rc, shift))
    interaction.setPotential(
        type1=monomer_type,
        type2=solvent_type,
        potential=espressopp.interaction.LennardJones(epsilon, sigma, rc, shift))
    system.addInteraction(interaction, 'lj')

    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = dt
    if (temperature != None):
        thermostat = espressopp.integrator.LangevinThermostat(system)
        thermostat.gamma = 1.0
        thermostat.temperature = temperature
        integrator.addExtension(thermostat)

    mass = 1.0

    if xyzfilename:
        props = ['id', 'type', 'mass', 'pos', 'v']
        bondlist = espressopp.FixedPairList(system.storage)
        particles = []
        bonds = []
        for i in xrange(num_chains):
            for k in xrange(monomers_per_chain):
                idx = i * monomers_per_chain + k
                part = [pidf[idx], typef[idx], mass,
                        espressopp.Real3D(xposf[idx], yposf[idx], zposf[idx]),
                        espressopp.Real3D(xvelf[idx], yvelf[idx], zvelf[idx])]
                particles.append(part)
                if k > 0:
                    bonds.append((pidf[idx - 1], pidf[idx]))
        # Read solvent
        for i in xrange(0, nsolvents):
            idx = num_chains * monomers_per_chain + i
            part = [pidf[idx], typef[idx], mass,
                    espressopp.Real3D(xposf[idx], yposf[idx], zposf[idx]),
                    espressopp.Real3D(xvelf[idx], yvelf[idx], zvelf[idx])]
            particles.append(part)
        system.storage.addParticles(particles, *props)
        system.storage.decompose()
        vl.exclude(bonds)
        bondlist.addBonds(bonds)
    else:
        props = ['id', 'type', 'mass', 'pos', 'v']
        vel_zero = espressopp.Real3D(0.0, 0.0, 0.0)
        bondlist = espressopp.FixedPairList(system.storage)
        pid = 1
        particles = []
        bonds = []
        for i in xrange(num_chains):
            startpos = system.bc.getRandomPos()
            positions, b = espressopp.tools.topology.polymerRW(
                pid, startpos, monomers_per_chain, 1.2)
            for k in xrange(monomers_per_chain):
                part = [pid + k, monomer_type, mass, positions[k], vel_zero]
                particles.append(part)
            pid += monomers_per_chain
            bonds.extend(b)
        for i in xrange(1, nsolvents + 1):
            startpos = system.bc.getRandomPos()
            particles.append([pid, solvent_type, mass, startpos, vel_zero])
            pid += 1
        system.storage.addParticles(particles, *props)
        system.storage.decompose()
        vl.exclude(bonds)
        bondlist.addBonds(bonds)

    # FENE bonds
    potFENE=espressopp.interaction.FENELennardJones(K=30.0, r0=0.0, rMax=1.5)
    interFENE=espressopp.interaction.FixedPairListFENELennardJones(
        system, bondlist, potFENE)
    system.addInteraction(interFENE, 'fene')

    return system, integrator, vl
