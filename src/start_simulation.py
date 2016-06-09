#!/usr/bin/env python
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

import espressopp  # NOQA
import h5py
import math  # NOQA
try:
    import MPI
except ImportError:
    from mpi4py import MPI
import time
import logging
import random
import shutil

import chemlab
import tools

import os

h5md_group = 'atoms'

__doc__ = 'Run GROMACS-like simulation with chemical reactions'


def main():  #NOQA
    args = tools._args().parse_args()

    tools._args().save_to_file('{}params.out'.format(args.output_prefix), args)

    # GROMACS units, kJ/mol K
    kb = 0.0083144621
    mass_factor = 1.6605402

    if args.kb:
        kb = args.kb

    if args.mass_factor:
        mass_factor = args.mass_factor

    print('Using kB={} and mass-factor={}'.format(kb, mass_factor))

    if args.debug:
        for s in args.debug.split(','):
            print('Activating logger {}'.format(s))
            logging.getLogger(s.strip()).setLevel(logging.DEBUG)

    lj_cutoff = args.lj_cutoff
    cg_cutoff = args.cg_cutoff
    max_cutoff = max([lj_cutoff, cg_cutoff])
    dt = args.dt

    time0 = time.time()

    gt = chemlab.gromacs_topology.GromacsTopology(args.top)
    gt.read()

    input_conf = chemlab.files_io.GROFile(args.conf)
    input_conf.read()

    box = input_conf.box
    print('Setup simulation...')

    # Tune simulation parameter according to arguments
    integrator_step = min([args.int_step, args.trj_collect])
    sim_step = args.run / integrator_step

    if args.skin:
        skin = args.skin

    # Seed for RNG
    rng_seed = args.rng_seed
    if not args.rng_seed or args.rng_seed == -1:
        rng_seed = random.randint(10, 1000000)

    print('Skin: {}'.format(skin))
    print('RNG Seed: {}'.format(rng_seed))
    print('Boltzmann constant: {}'.format(kb))

    part_prop, particle_list = chemlab.gromacs_topology.gen_particle_list(input_conf, gt)
    NPart = len(particle_list)
    print('Reads {} particles with properties {}'.format(NPart, part_prop))

    particle_ids = [x[0] for x in particle_list]

    density = sum(x[3] for x in particle_list)*mass_factor/ (box[0] * box[1] * box[2])
    print('Density: {} kg/m^3'.format(density))
    print('Box: {} nm'.format(box))

    # Generate velocity.
    print('Generating velocities from Maxwell-Boltzmann distribution T={} ({})'.format(
        args.temperature, args.temperature*kb))
    vx, vy, vz = espressopp.tools.velocities.gaussian(
        args.temperature, len(particle_list), [x[3]*mass_factor for x in particle_list],
        kb=kb)
    part_prop.append('v')
    for i, p in enumerate(particle_list):
        p.append(espressopp.Real3D(vx[i], vy[i], vz[i]))

    system = espressopp.System()
    system.rng = espressopp.esutil.RNG(rng_seed)
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.skin = skin
    if args.node_grid:
        nodeGrid = map(int, args.node_grid.split(','))
    else:
        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    print('Number of nodes {}, node-grid: {}'.format(
        MPI.COMM_WORLD.size, nodeGrid))
    cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, max_cutoff, skin)

    print('Cell grid: {}'.format(cellGrid))

    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = dt
    system.integrator = integrator

    system.storage.addParticles(particle_list, *part_prop)
    system.storage.decompose()

# Dynamic exclude list, depends on the new create bonds as well.
    dynamic_exclusion_list = espressopp.DynamicExcludeList(integrator, gt.exclusions)
    print('Excluded pairs from LJ interaction: {}'.format(len(gt.exclusions)))

# Exclude all bonded interaction from the lennard jones
    verletlist = espressopp.VerletList(
        system,
        cutoff=max_cutoff,
        exclusionlist=dynamic_exclusion_list
        )

# define the potential, interaction_id = 0
    print('Bonds: {}'.format(len(gt.bonds)))
    print('Angles: {}'.format(len(gt.angles)))
    print('Dihedrals: {}'.format(len(gt.dihedrals)))

# Define the thermostat
    temperature = args.temperature*kb
    print('Temperature: {} ({}), gamma: {}'.format(args.temperature, temperature, args.thermostat_gamma))
    print('Thermostat: {}'.format(args.thermostat))
    if args.thermostat == 'lv':
        thermostat = espressopp.integrator.LangevinThermostat(system)
        thermostat.temperature = temperature
        thermostat.gamma = args.thermostat_gamma
    elif args.thermostat == 'vr':
        thermostat = espressopp.integrator.StochasticVelocityRescaling(system)
        thermostat.temperature = temperature
        thermostat.coupling = args.thermostat_gamma
    else:
        raise Exception('Wrong thermostat keyword: `{}`'.format(args.thermostat))
    integrator.addExtension(thermostat)

# Pressure coupling if needed,
    pressure_comp = espressopp.analysis.Pressure(system)
    pressure = 0.0
    if args.pressure:
        pressure = args.pressure * 0.060221374  # convert from bars to gromacs units kj/mol/nm^3
        if args.barostat == 'lv':
            print('Barostat: Langevin with P={}, gamma={}, mass={}'.format(
                pressure, 0.5, pow(10, 4)))
            barostat = espressopp.integrator.LangevinBarostat(system, system.rng, temperature)
            barostat.gammaP = args.barostat_gammaP
            barostat.mass = args.barostat_mass
            barostat.pressure = pressure
        elif args.barostat == 'br':
            print('Barostat: Berendsen with P={} and tau={}'.format(pressure, 0.5))
            barostat = espressopp.integrator.BerendsenBarostat(system, pressure_comp)
            barostat.tau = args.barostat_tau
            barostat.pressure = pressure
        else:
            raise Exception('Wrong barostat keyword: `{}`'.format(args.barostat))
        integrator.addExtension(barostat)

    print("Decomposing now ...")
    system.storage.decompose()

    # Conditional break of the reactions.
    cr_observs = None
    if args.maximum_conversion:
        pass
        #TODO(jakub): finish this part. Stop simulation whenever it's reach given conversion rate.

    print('Set topology manager')
    topology_manager = espressopp.integrator.TopologyManager(system)

    # Set potentials.
    cr_observs = chemlab.gromacs_topology.setNonbondedInteractions(
        system, gt, verletlist, lj_cutoff, cg_cutoff, tables=args.table_groups, cr_observs=cr_observs)
    dynamic_fpls, static_fpls = chemlab.gromacs_topology.set_bonded_interactions(system, gt)
    dynamic_ftls, static_ftls = chemlab.gromacs_topology.set_angle_interactions(system, gt)
    dynamic_fqls, static_fqls = chemlab.gromacs_topology.set_dihedral_interactions(system, gt)

    dynamic_fpairs, static_fpairs = chemlab.gromacs_topology.set_pair_interactions(system, gt, args)
    chemlab.gromacs_topology.set_coulomb_interactions(system, gt, args)

    print('Set Dynamic Exclusion lists.')
    for static_fpl in static_fpls:
        dynamic_exclusion_list.observe_tuple(static_fpl)
    for static_ftl in static_ftls:
        dynamic_exclusion_list.observe_triple(static_ftl)
    for static_fql in static_fqls:
        dynamic_exclusion_list.observe_quadruple(static_fql)

    for _, fpl in dynamic_fpls.items():
        dynamic_exclusion_list.observe_tuple(fpl)
    for _, ftl in dynamic_ftls.items():
        dynamic_exclusion_list.observe_triple(ftl)
    for _, fql in dynamic_fqls.items():
        dynamic_exclusion_list.observe_quadruple(fql)

    print('Set dynamic topology')
    for static_fpl in static_fpls:
        topology_manager.observe_tuple(static_fpl)
    for _, fpl in dynamic_fpls.items():
        topology_manager.observe_tuple(fpl)

    topology_manager.initialize_topology()

    for t, p in gt.bondparams.items():
        if p['func'] in dynamic_fpls:
            fpl = dynamic_fpls[p['func']]
            print('Register bonds for type: {}'.format(t))
            topology_manager.register_tuple(fpl, *t)

    for t, p in gt.angleparams.items():
        if p['func'] in dynamic_ftls:
            ftl = dynamic_ftls[p['func']]
            print('Register angles for type: {}'.format(t))
            topology_manager.register_triplet(ftl, *t)

    for t, p in gt.dihedralparams.items():
        if p['func'] in dynamic_fqls:
            fql = dynamic_fqls[p['func']]
            print('Register dihedral for type: {}'.format(t))
            topology_manager.register_quadruplet(fql, *t)

    integrator.addExtension(topology_manager)

    # Set chemical reactions, parser in reaction_parser.py
    fpls = []
    cr_interval = 0
    if args.reactions:
        if os.path.exists(args.reactions):
            print('Set chemical reactions from: {}'.format(args.reactions))
            reaction_config = chemlab.reaction_parser.parse_config(args.reactions)
            sc = chemlab.reaction_parser.SetupReactions(
                system, verletlist, gt, topology_manager, reaction_config)

            ar, fpls = sc.setup_reactions()
            output_reaction_config = '{}_{}_{}'.format(args.output_prefix, rng_seed, args.reactions)
            print('Save copy of reaction config to: {}'.format(output_reaction_config))
            shutil.copyfile(args.reactions, output_reaction_config)
            cr_interval = sc.ar_interval
            integrator_step = min(cr_interval, integrator_step)
            sim_step = args.run / integrator_step
            args.topol_collect = cr_interval
    else:
        cr_interval = integrator_step

    cr_interval = min([integrator_step, cr_interval])

    for f in fpls:
        topology_manager.observe_tuple(f)
        dynamic_exclusion_list.observe_tuple(f)
    # Define SystemMonitor that will store data from observables into a .csv file.
    energy_file = '{}_energy_{}.csv'.format(args.output_prefix, rng_seed)
    print('Energy saved to: {}'.format(energy_file))
    system_analysis = espressopp.analysis.SystemMonitor(
        system,
        integrator,
        espressopp.analysis.SystemMonitorOutputCSV(energy_file))
    temp_comp = espressopp.analysis.Temperature(system)
    system_analysis.add_observable('T', temp_comp)
    system_analysis.add_observable(
        'Ekin', espressopp.analysis.KineticEnergy(
            system, temp_comp))
    for label, interaction in sorted(system.getAllInteractions().items()):
        print('System analysis: adding {}'.format(label))
        system_analysis.add_observable(
            label, espressopp.analysis.PotentialEnergy(system, interaction))
    for (cr_type, _), obs in cr_observs.items():
        system_analysis.add_observable(
            'cr_{}'.format(cr_type), obs)
    for fidx, f in enumerate(fpls):
        system_analysis.add_observable(
            'count_{}'.format(fidx), espressopp.analysis.NFixedPairListEntries(system, f))
    ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, cr_interval)
    integrator.addExtension(ext_analysis)
    print('Configured system analysis')

    print('Configure H5MD trajectory writer')
    h5md_output_file = '{}_{}_traj.h5'.format(args.output_prefix, rng_seed)
    traj_file = espressopp.io.DumpH5MD(
        system,
        h5md_output_file,
        group_name='atoms',
        static_box=False,
        author='XXX',
        email='xxx',
        store_species=args.store_species,
        store_res_id=True,
        store_charge=False,
        store_state=args.store_state,
        store_lambda=args.store_lambda,
        store_force=args.store_force,
        chunk_size=int(NPart/MPI.COMM_WORLD.size))

    print('Set topology writer')
    dump_topol = espressopp.io.DumpTopology(system, integrator, traj_file)
    for i, f in enumerate(fpls):
        dump_topol.observe_tuple(f, 'chem_bonds_{}'.format(i))

    for i, f in dynamic_fpls.items():
        dump_topol.observe_tuple(f, 'dynamic_bonds_{}'.format(i))

    bcount = 0
    for static_fpl in static_fpls:
        dump_topol.add_static_tuple(static_fpl, 'bonds_{}'.format(bcount))
        bcount += 1
    for dynamic_fpl in dynamic_fpls.values():
        dump_topol.add_static_tuple(dynamic_fpl, 'bonds_{}'.format(bcount))
        bcount += 1

    if args.topol_collect > 0:
        print('Collect topology: {}'.format(args.topol_collect))
        ext_dump = espressopp.integrator.ExtAnalyze(dump_topol, args.topol_collect)
        integrator.addExtension(ext_dump)

    k_trj_collect = int(math.ceil(args.trj_collect/float(integrator_step)))
    k_trj_flush = 10 if 10 < k_trj_collect else k_trj_collect
    print('Store trajectory every {} steps'.format(args.trj_collect))

    if args.start_ar >= 0:
        k_enable_reactions = int(math.ceil(args.start_ar/float(integrator_step)))
        print('Enable chemical reactions at {} step'.format(args.start_ar))
    else:
        k_enable_reactions = -1

    print('Reset total velocity')
    total_velocity = espressopp.analysis.TotalVelocity(system)
    total_velocity.reset()

    if args.max_force > -1:
        cap_force = espressopp.integrator.CapForce(system, args.max_force)
        integrator.addExtension(cap_force)
        print('Cap force to {}'.format(args.max_force))

    print('Mapping Type name  type id')
    for at_sym in gt.used_atomtypes:
        print('        {:9}  {:8}'.format(at_sym, gt.atomsym_atomtype[at_sym]))

    print('Running {} steps'.format(sim_step*integrator_step))
    print('Temperature: {} ({} K)'.format(args.temperature*kb, args.temperature))
    system_analysis.dump()

    for k in range(sim_step):
        system_analysis.info()
        if k % k_trj_collect == 0:
            traj_file.dump(k*integrator_step, k*integrator_step*args.dt)
        if k % k_trj_flush == 0:
            dump_topol.update()
            traj_file.flush()   # Write HDF5 to disk.
        if k_enable_reactions == k:
            print('Enabling chemical reactions')
            integrator.addExtension(ar)
        integrator.run(integrator_step)

    system_analysis.info()
    traj_file.dump(sim_step*integrator_step, sim_step*integrator_step*args.dt)
    dump_topol.dump()
    dump_topol.update()
    traj_file.flush()
    traj_file.close()

    # Write some parameters of the simulation.
    h5 = h5py.File(h5md_output_file, 'r+')
    if 'parameters' not in h5:
        h5.create_group('parameters')
    g_params = h5['/parameters']
    sim_params = {
        'thermostat': args.thermostat,
        'thermostat_gamma': args.thermostat_gamma,
        'temperature': args.temperature,
        'kb': kb,
        'barostat': args.barostat if args.pressure else 'no',
        'pressure': pressure,
        'total_steps': sim_step*integrator_step,
        'total_time': sim_step*integrator_step*args.dt
    }
    for k, v in sim_params.items():
        g_params.attrs[k] = v
    tools.save_forcefield(h5, gt)
    h5.close()

    # Saves output file.
    output_gro_file = '{}_{}_confout.gro'.format(args.output_prefix, rng_seed)
    dump_gro = espressopp.io.DumpGRO(
        system, integrator, filename=output_gro_file,
        unfolded=True, append=False)
    dump_gro.dump()
    print('Wrote end configuration to: {}'.format(output_gro_file))

    print('finished!')
    print('total time: {}'.format(time.time()-time0))
    print('Some timers:')
    print topology_manager.get_timers()

    print traj_file.getTimers()
    espressopp.tools.analyse.final_info(system, integrator, verletlist, time0, time.time())


if __name__ == '__main__':
    main()
