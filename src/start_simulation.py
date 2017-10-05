#!/usr/bin/env python
#  Copyright (C) 2016-2017
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
import numpy

try:
    import MPI
except ImportError:
    from mpi4py import MPI
import collections
import cPickle
import time
import logging
import random
import re
import shutil

import chemlab
import tools
import app_args

import os

h5md_group = 'atoms'

__doc__ = 'Run GROMACS-like simulation with chemical reactions'


def main():  # NOQA
    args = app_args._args().parse_args()

    app_args._args().save_to_file('{}params.out'.format(args.output_prefix), args)

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
            name_filter = s.split(':')
            print('Activating logger {}'.format(name_filter[0]))
            logging.getLogger(name_filter[0].strip()).setLevel(logging.DEBUG)
            if len(name_filter) == 2:
                log_filter = tools.RegexpFilter(name_filter[1])
                logging.getLogger(name_filter[0].strip()).addFilter(log_filter)

    if args.check_topology:
        logging.getLogger('TopologyManager').setLevel(logging.WARN)

    lj_cutoff = args.lj_cutoff
    cg_cutoff = args.cg_cutoff
    max_cutoff = max([lj_cutoff, cg_cutoff])
    dt = args.dt
    print('LJ cutoff: {} Tabulated cutoff: {} time-step: {}'.format(
        lj_cutoff, cg_cutoff, dt))

    time0 = time.time()

    # Check if exclusion list file is available
    has_exclusions = args.exclusion_list is not None and os.path.exists(args.exclusion_list)

    gt = chemlab.gromacs_topology.GromacsTopology(args.top, generate_exclusions=not has_exclusions)
    gt.read()

    input_conf = chemlab.files_io.GROFile(args.conf)
    input_conf.read()

    box = input_conf.box
    print('Setup simulation...')

    # Tune simulation parameter according to arguments
    integrator_step = args.int_step
    if args.trj_collect > 0:
        integrator_step = min([args.int_step, args.trj_collect])
    sim_step = args.run / integrator_step

    if args.skin == 'auto':
        skin = 0.16
    else:
        skin = float(args.skin)

    # Seed for RNG
    rng_seed = args.rng_seed
    if not args.rng_seed or args.rng_seed == -1:
        rng_seed = random.randint(10, 1000000)
        args.rng_seed = rng_seed

    global_file_prefix = '{}_{}'.format(args.output_prefix, rng_seed)

    print('Skin: {}'.format(skin))
    print('RNG Seed: {}'.format(rng_seed))
    print('Boltzmann constant: {}'.format(kb))

    part_prop, particle_list = chemlab.gromacs_topology.gen_particle_list(input_conf, gt)
    NPart = len(particle_list)
    print('Reads {} particles with properties {}'.format(NPart, part_prop))

    particle_ids = [x[0] for x in particle_list]

    density = sum(x[3] for x in particle_list) * mass_factor / (box[0] * box[1] * box[2])
    print('Density: {} kg/m^3'.format(density))
    print('Box: {} nm'.format(box))

    # Generate velocity.
    if args.temperature is None:
        raise RuntimeError('Temperature not defined!')
    temperature = args.temperature * kb
    if args.gen_velocity:
        print('Generating velocities from Maxwell-Boltzmann distribution T={} ({})'.format(
            args.temperature, args.temperature * kb))
        vx, vy, vz = espressopp.tools.velocities.gaussian(
            args.temperature,
            len(particle_list),
            [x[3] * mass_factor for x in particle_list],
            kb=kb)
        part_prop.append('v')
        for i, p in enumerate(particle_list):
            p.append(espressopp.Real3D(vx[i], vy[i], vz[i]))

    system = espressopp.System()
    system.rng = espressopp.esutil.RNG(rng_seed)

    system.skin = skin
    if args.node_grid:
        nodeGrid = map(int, args.node_grid.split(','))
    else:
        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    print('Number of nodes {}, node-grid: {}'.format(
        MPI.COMM_WORLD.size, nodeGrid))
    cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, max_cutoff, skin)

    print('Cell grid: {}'.format(cellGrid))

    system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = dt
    system.integrator = integrator

    system.storage.addParticles(particle_list, *part_prop)

    system.storage.decompose()

    # Dynamic exclude list, depends on the new create bonds as well.
    if has_exclusions:
        exclusion_file = open(args.exclusion_list, 'r')
        exclusions = [map(int, x.split()) for x in exclusion_file.readlines()]
        print('Read exclusion list from {} (total: {})'.format(args.exclusion_list, len(exclusions)))
        if len(exclusions) == 0 and gt.bonds and any([x > 0 for x in gt.gt.moleculetype.values()]):
            raise RuntimeError('Exclusion list in {} is empty'.format(args.exclusion_list))
        gt.exclusions = exclusions

    topol_file = os.path.basename(args.top)
    output_filename = 'exclusion_{}.list'.format(topol_file.split('.')[0])
    out_file = open(output_filename, 'w')
    out_file.writelines('\n'.join(['{} {}'.format(*d) for d in sorted(gt.exclusions)]))
    out_file.close()

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

    print("Decomposing now ...")
    system.storage.decompose()

    # Conditional break of the reactions.
    cr_observs = None

    print('Set topology manager')
    topology_manager = espressopp.integrator.TopologyManager(system)

    # Hooks
    hook_init_reaction = lambda *_, **__: True
    hook_postsetup_reaction = lambda *_, **__: True
    hook_at_step = lambda *_, **__: True
    hook_before_sim = lambda *_, **__: True
    if os.path.exists('hooks.py'):
        print('Found hooks.py')
        locals = {}
        execfile('hooks.py', globals(), locals)
        hook_init_reaction = locals.get('hook_init_reaction', hook_init_reaction)
        hook_postsetup_reaction = locals.get('hook_postsetup_reaction', hook_postsetup_reaction)
        hook_at_step = locals.get('hook_at_step', hook_at_step)
        hook_before_sim = locals.get('hook_before_sim', hook_before_sim)

    # Set chemical reactions, parser in reaction_parser.py
    chem_dynamic_types = set()
    chem_dynamic_bond_types = set()
    separate_fpls = set()
    chem_fpls = []
    reactions = []
    extensions_integrator = []
    cr_interval = 0
    has_reaction = False
    sc = None
    ar = None
    if args.reactions is not None and os.path.exists(args.reactions):
        print('Set chemical reactions from: {}'.format(args.reactions))
        reaction_config = chemlab.reaction_parser.parse_config(args.reactions)
        sc = chemlab.reaction_setup.SetupReactions(
            system,
            verletlist,
            gt,
            topology_manager,
            reaction_config,
            args)
        ar, chem_fpls, reactions, extensions_integrator = sc.setup_reactions()
        chem_dynamic_types = sc.dynamic_types
        chem_dynamic_bond_types = sc.observed_bondtypes
        separate_fpls = sc.separate_fpls

        if cr_observs is None:
            cr_observs = {}
        if sc.cr_observs is not None:
            cr_observs.update(sc.cr_observs)

        output_reaction_config = '{}_{}_{}'.format(args.output_prefix, rng_seed, args.reactions)
        print('Save copy of reaction config to: {}'.format(output_reaction_config))
        shutil.copyfile(args.reactions, output_reaction_config)

        cr_interval = sc.ar_interval
        integrator_step = min(cr_interval, integrator_step)
        print('Change integrator step to {}'.format(integrator_step))
        sim_step = args.run / integrator_step
        print('Change topology collect interval to {}'.format(cr_interval))
        args.topol_collect = min([cr_interval, args.topol_collect])
        has_reaction = True
        hook_postsetup_reaction(system, integrator, gt, args, ar)
    else:
        cr_interval = integrator_step

    print('Dynamic type ids: {}'.format(chem_dynamic_types))

    cr_interval = min([integrator_step, cr_interval])

    maximum_conversion = []
    eq_run = 0
    if args.maximum_conversion:
        if cr_observs is None:
            cr_observs = {}
        maximum_conversion = tools.get_maximum_conversion(args, system, chem_fpls, gt, cr_observs)
        if args.eq_steps > 0:
            eq_run = int(args.eq_steps / sim_step)

    if args.t_hybrid_bond > 0:
        list_dynamic_resolution = espressopp.integrator.FixedListDynamicResolution(system)
        for fpl in chem_fpls:
            list_dynamic_resolution.register_pair_list(fpl.fpl, 1.0 / args.t_hybrid_bond)
        integrator.addExtension(list_dynamic_resolution)

    system.storage.decompose()

    # Set potentials.
    cr_observs, particle_pair_scales = chemlab.gromacs_topology.set_nonbonded_interactions(
        system, gt, verletlist, lj_cutoff, cg_cutoff, tables=args.table_groups, cr_observs=cr_observs)
    dynamic_fpls, static_fpls, registered_fpls = chemlab.gromacs_topology.set_bonded_interactions(
        system, gt, chem_dynamic_types, chem_dynamic_bond_types, separate_fpls)
    dynamic_ftls, static_ftls = chemlab.gromacs_topology.set_angle_interactions(
        system, gt, chem_dynamic_types, chem_dynamic_bond_types)
    dynamic_fqls, static_fqls = chemlab.gromacs_topology.set_dihedral_interactions(
        system, gt, chem_dynamic_types, chem_dynamic_bond_types)

    dynamic_fpairs, static_fpairs = chemlab.gromacs_topology.set_pair_interactions(
        system, gt, args, chem_dynamic_types)
    chemlab.gromacs_topology.set_coulomb_interactions(system, gt, args)

    if args.thermal_groups:
        thermal_groups = map(gt.atomsym_atomtype.get, args.thermal_groups.split(','))
    elif args.table_groups:
        thermal_groups = map(gt.atomsym_atomtype.get, args.table_groups.split(','))
        print('Thermal groups: {} ({})'.format(args.table_groups, thermal_groups))
    else:
        thermal_groups = []

    # Add cap force
    if args.max_force > -1:
        cap_force = espressopp.integrator.CapForce(system, args.max_force)
        integrator.addExtension(cap_force)
        print('Cap force to {}'.format(args.max_force))

    # Define the thermostat
    print('Temperature: {} ({}), gamma: {}'.format(args.temperature, temperature, args.thermostat_gamma))
    print('Thermostat: {}'.format(args.thermostat))
    thermostat = None
    if args.thermostat == 'lv':
        thermostat = espressopp.integrator.LangevinThermostat(system)
        thermostat.temperature = temperature
        thermostat.gamma = args.thermostat_gamma
        if has_reaction and sc and sc.use_thermal_group:
            print('Running thermostat on thermal groups: {}'.format(thermal_groups))
            thermostat.add_valid_types(thermal_groups)
    elif args.thermostat == 'vr':
        thermostat = espressopp.integrator.StochasticVelocityRescaling(system)
        thermostat.temperature = temperature
        thermostat.coupling = args.thermostat_gamma
    elif args.thermostat == 'iso':
        thermostat = espressopp.integrator.Isokinetic(system)
        thermostat.temperature = temperature
        thermostat.coupling = int(args.thermostat_gamma)
    elif args.thermostat == 'no':
        print('No thermostat selected, runing NVE simulation?')
    else:
        raise Exception('Wrong thermostat keyword: `{}`'.format(args.thermostat))
    if thermostat is not None:
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
    # Observe tuples, any new bond here trigger new angles, dihedrals.
    for static_fpl in static_fpls:
        print('Observe tuple {}'.format(static_fpl))
        topology_manager.observe_tuple(static_fpl)
    for _, fpl in dynamic_fpls.items():
        topology_manager.observe_tuple(fpl)

    topology_manager.initialize_topology()

    # for t, p in gt.bondparams.items():
    #     fpl = dynamic_fpls.get((p['func'], True), dynamic_fpls.get((p['func'], False)))
    #     if fpl:
    #         print('Register bonds for type: {}'.format(t))
    #         topology_manager.register_tuple(fpl, *t)

    # Any new bond will trigger update here.
    for t, p in gt.angleparams.items():
        ftl = dynamic_ftls.get((p['func'], True), dynamic_ftls.get((p['func'], False)))
        if ftl:
            print('Register angles for type: {} ({})'.format(t, ftl))
            topology_manager.register_triplet(ftl, *t)

    for t, p in gt.dihedralparams.items():
        fql = dynamic_fqls.get((p['func'], True), dynamic_fqls.get((p['func'], False)))
        if fql:
            print('Register dihedral for type: {} ({})'.format(t, fql))
            topology_manager.register_quadruplet(fql, *t)

    integrator.addExtension(topology_manager)

    # Add fpls defined by chemical reactions to exclude lists and topology_manager.
    for f in chem_fpls:
        print("TopologyManager: Observe tuple {}".format(f.fpl))
        topology_manager.observe_tuple(f.fpl)
        dynamic_exclusion_list.observe_tuple(f.fpl)

    # Register chemistry tuples in topology_manager
    for def_f in chem_fpls:
        for t in def_f.type_list:
            print('Register chem_fpl for type: {} ({})'.format(t, def_f.fpl))
            topology_manager.register_tuple(def_f.fpl, *t)

    for ptypes, fpl in registered_fpls:
        print('Register fpls for type: {}-{} ({})'.format(ptypes[0], ptypes[1], fpl))
        topology_manager.register_tuple(fpl, *ptypes)
        dynamic_exclusion_list.observe_tuple(fpl)
    # All data defined in topology manager, we can rebuild reactions
    if has_reaction:
        sc.rebuild_fixed_pair_lists()

    # Define SystemMonitor that will store data from observables into a .csv file.
    energy_file = '{}_energy_{}.csv'.format(args.output_prefix, rng_seed)
    print('Energy saved to: {}'.format(energy_file))
    system_analysis = espressopp.analysis.SystemMonitor(
        system,
        integrator,
        espressopp.analysis.SystemMonitorOutputCSV(energy_file))
    temp_comp = espressopp.analysis.Temperature(system)
    if has_reaction and sc and sc.use_thermal_group:
        for t_id in thermal_groups:
            temp_comp.add_type(t_id)

    system_monitor_filter = None
    if args.system_monitor_filter:
        system_monitor_filter = args.system_monitor_filter.split(',')

    system_analysis.add_observable('T', temp_comp)
    system_analysis.add_observable(
        'Ekin', espressopp.analysis.KineticEnergy(
            system, temp_comp))
    if args.store_pressure:
        system_analysis.add_observable('P', pressure_comp)
    for label, interaction in sorted(system.getAllInteractions().items()):
        print('System analysis: adding {}'.format(label))
        show_in_system_info = True
        if system_monitor_filter:
            show_in_system_info = False
            for v in system_monitor_filter:
                if v in label:
                    show_in_system_info = True
                    break
        system_analysis.add_observable(
            label, espressopp.analysis.PotentialEnergy(system, interaction), show_in_system_info)
    for (cr_type, _, ts), obs in cr_observs.items():
        if not isinstance(cr_type, collections.Iterable):
            cr_type = [cr_type]
        if ts is None:
            system_analysis.add_observable(
                'cr_{}'.format('_'.join(map(str, cr_type))), obs)
        else:
            system_analysis.add_observable(
                'cr_{}_{}'.format('_'.join(map(str, cr_type)), ts), obs)

    for fidx, f in enumerate(chem_fpls):
        system_analysis.add_observable(
            'count_{}'.format(fidx), espressopp.analysis.NFixedPairListEntries(system, f.fpl))

    if args.t_hybrid_bond > 0:
        for fpl_idx, fpl in enumerate(chem_fpls):
            system_analysis.add_observable(
                'res_fpl_{}'.format(fpl_idx), espressopp.analysis.ResolutionFixedPairList(system, fpl.fpl))

    # system_analysis.add_observable('Fmax', espressopp.analysis.MaxForce(system))

    # This is a bit expensive
    if args.count_tuples:
        bcount = 0
        for static_fpl in static_fpls:
            system_analysis.add_observable(
                'bcount_{}'.format(bcount), espressopp.analysis.NFixedPairListEntries(system, static_fpl))
            bcount += 1
        for _, fpl in dynamic_fpls.items():
            system_analysis.add_observable(
                'bcount_{}'.format(bcount), espressopp.analysis.NFixedPairListEntries(system, fpl))
            bcount += 1

        for ptypes, fpls in registered_fpls:
            system_analysis.add_observable(
                'bcount_{}-{}'.format(*ptypes), espressopp.analysis.NFixedPairListEntries(system, fpls))

        bcount = 0
        for static_ftl in static_ftls:
            system_analysis.add_observable(
                'acount_{}'.format(bcount), espressopp.analysis.NFixedTripleListEntries(system, static_ftl))
            bcount += 1
        for _, ftl in dynamic_ftls.items():
            system_analysis.add_observable(
                'acount_{}'.format(bcount), espressopp.analysis.NFixedTripleListEntries(system, ftl))
            bcount += 1
        bcount = 0
        for static_fql in static_fqls:
            system_analysis.add_observable(
                'qcount_{}'.format(bcount), espressopp.analysis.NFixedQuadrupleListEntries(system, static_fql))
            bcount += 1
        for _, fql in dynamic_fqls.items():
            system_analysis.add_observable(
                'qcount_{}'.format(bcount), espressopp.analysis.NFixedQuadrupleListEntries(system, fql))
            bcount += 1
        # Add counter on the exclude list pairs
        system_analysis.add_observable('vl_excl', espressopp.analysis.NExcludeListEntries(system, verletlist))

        # Observe ParticlePairsScale
        for pps_idx, pps in enumerate(particle_pair_scales, 1):
            system_analysis.add_observable(
                'pair_scale_{}'.format(pps_idx), espressopp.analysis.NParticlePairScalingEntries(system, pps))

    if args.count_types:
        for at_sym in args.count_types.split(','):
            print('Observer {:9} ({})'.format(at_sym, gt.atomsym_atomtype[at_sym]))
            obs_type_id = gt.atomsym_atomtype[at_sym]
            chem_conver_obs = espressopp.analysis.ChemicalConversion(system, obs_type_id)
            system_analysis.add_observable('num_type_{}_{}'.format(at_sym, obs_type_id), chem_conver_obs)

    if args.count_types_state is not None:
        types_state = args.count_types_state.split(',')
        for ts in types_state:
            type_name, state = ts.split(':')
            type_id = gt.atomsym_atomtype[type_name]
            state = int(state)
            system_analysis.add_observable(
                'st_{}_{}'.format(type_name, state),
                espressopp.analysis.ChemicalConversionTypeState(system, type_id, state))

    if args.count_fix_distances:
        for fd_idx, fd in enumerate(sc.fix_distances):
            system_analysis.add_observable('fd_{}'.format(fd_idx), espressopp.analysis.NumFixDistances(system, fd))

    cr_interval = min([cr_interval, args.energy_collect])
    if args.energy_collect > 0:
        ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, min([cr_interval, args.energy_collect]))
        integrator.addExtension(ext_analysis)
        print('Configured system analysis, collect data every {} steps'.format(min([cr_interval, args.energy_collect])))

    print('Configure H5MD trajectory writer')
    NPart = espressopp.analysis.NPart(system).compute()
    h5md_output_file = '{}_{}_traj.h5'.format(args.output_prefix, rng_seed)
    traj_file = espressopp.io.DumpH5MD(
        system,
        h5md_output_file,
        group_name='atoms',
        static_box=False if args.pressure else True,
        author='XXX',
        email='xxx',
        store_species=args.store_species,
        store_res_id=args.store_res_id,
        store_charge=args.store_charge,
        store_position=args.store_position,
        store_state=args.store_state,
        store_lambda=args.store_lambda,
        store_force=args.store_force,
        store_velocity=args.store_velocity,
        store_mass=args.store_mass,
        is_single_prec=args.store_single_precision,
        chunk_size=256)  # int(NPart/MPI.COMM_WORLD.size))

    print('Set topology writer')
    dump_topol = espressopp.io.DumpTopology(system, integrator, traj_file)
    for i, f in enumerate(chem_fpls):
        dump_topol.observe_tuple(f.fpl, 'chem_bonds_{}'.format(i))

    bcount = acount = qcount = 0
    for (i, observe_tuple), f in dynamic_fpls.items():
        if observe_tuple and args.store_angdih:
            print('DumpTopol: observe dynamic_bonds_{}'.format(bcount))
            dump_topol.observe_tuple(f, 'dynamic_bonds_{}'.format(bcount))
        else:
            print('DumpTopol: save static list from bonds_{}'.format(bcount))
            dump_topol.add_static_tuple(f, 'bonds_{}'.format(bcount))
        bcount += 1

    for ptypes, fpl in registered_fpls:
        print('DumpTopol: observe fpls for type: {}-{} ({})'.format(ptypes[0], ptypes[1], fpl))
        dump_topol.observe_tuple(fpl, 'dynamic_bonds_{}_{}'.format(ptypes[0], ptypes[1]))

    for (i, observe_triple), f in dynamic_ftls.items():
        if observe_triple and args.store_angdih:
            print('DumpTopol: observe dynamic_angles_{}'.format(acount))
            dump_topol.observe_triple(f, 'dynamic_angles_{}'.format(acount))
        else:
            print('DumpTopol: save static list from angles_{}'.format(acount))
            dump_topol.add_static_triple(f, 'angles_{}'.format(acount))
        acount += 1

    for (i, observe_quadruple), f in dynamic_fqls.items():
        if observe_quadruple and args.store_angdih:
            print('DumpTopol: observe dynamic_dihedrals_{}'.format(qcount))
            dump_topol.observe_quadruple(f, 'dynamic_dihedrals_{}'.format(qcount))
        else:
            print('DumpTopol: save static list from dihedrals_{}'.format(qcount))
            dump_topol.add_static_quadruple(f, 'dihedrals_{}'.format(qcount))
        qcount += 1

    for static_fpl in static_fpls:
        print('DumpTopol: store bonds_{}'.format(bcount))
        dump_topol.add_static_tuple(static_fpl, 'bonds_{}'.format(bcount))
        bcount += 1

    for static_ftl in static_ftls:
        print('DumpTopol: store angles_{}'.format(acount))
        dump_topol.add_static_triple(static_ftl, 'angles_{}'.format(acount))
        acount += 1

    for static_fql in static_fqls:
        print('DumpTopol: store dihedrals_{}'.format(qcount))
        dump_topol.add_static_quadruple(static_fql, 'dihedrals_{}'.format(qcount))
        qcount += 1

    if args.start_ar >= 0 and has_reaction:
        k_enable_reactions = int(math.ceil(args.start_ar / float(integrator_step)))
        print('Enable chemical reactions at {} step'.format(args.start_ar))
    else:
        k_enable_reactions = -1
    save_traj_topology = args.save_before_reaction if k_enable_reactions > 2 else True

    if args.topol_collect > 0 and save_traj_topology:
        print('Collect topology every {} steps'.format(args.topol_collect))
        ext_dump = espressopp.integrator.ExtAnalyze(dump_topol, args.topol_collect)
        integrator.addExtension(ext_dump)
        dump_topol.dump()
        dump_topol.update()

    trj_collect = min([args.trj_collect, cr_interval]) if cr_interval > 0 else args.trj_collect
    k_trj_collect = int(math.ceil(trj_collect / float(integrator_step)))
    if args.trj_flush is None:
        k_trj_flush = 25 if 25 < 10 * k_trj_collect else 10 * k_trj_collect
    else:
        k_trj_flush = int(math.ceil(args.trj_flush / float(integrator_step)))
    print('Collect trajectory every {} steps'.format(trj_collect))
    print('Collect energy data every {} steps'.format(cr_interval))
    print('Flush trajectory and topology to disk every {} steps'.format(k_trj_flush * integrator_step))

    if args.stop_ar >= 0 and has_reaction:
        k_stop_reactions = int(math.ceil(args.stop_ar / float(integrator_step)))
        print('Disable reactions at {} step'.format(args.stop_ar))
    else:
        k_stop_reactions = -1

    if args.rate_arrhenius:
        print(('Warning! Rate will change based on the arrhenius law. Keep in mind that it will'
               ' change rate for every reactions defined in configuration file based on the'
               ' global value of the total energy divided by the number of created bonds.'))

    print('Reset total velocity')
    total_velocity = espressopp.analysis.CMVelocity(system)
    total_velocity.reset()

    if args.gro_trj_collect:
        dump_gro_trj_fname = '{}_{}_traj.gro'.format(args.output_prefix, rng_seed)
        dump_gro_trj = espressopp.io.DumpGRO(
            system,
            integrator,
            filename=dump_gro_trj_fname,
            unfolded=True,
            append=True)
        ext_dump_gro = espressopp.integrator.ExtAnalyze(dump_gro_trj, args.gro_trj_collect)
        integrator.addExtension(ext_dump_gro)
        print('Set gro trajectory saver, save every {} steps'.format(args.gro_trj_collect))
        print('Warning, this will slow down simulation.')
        print('File saved to {}'.format(dump_gro_trj_fname))

    print('{:9}    {:8}'.format('Type name', 'type id'))
    for at_sym, type_id in sorted(gt.atomsym_atomtype.items(), key=lambda x: x[1]):
        print('{:9}    {:8}'.format(at_sym, type_id))

    print('Running {} steps'.format(sim_step * integrator_step))
    print('Temperature: {} ({} K)'.format(args.temperature * kb, args.temperature))
    print('Number of particles: {}'.format(NPart))
    system_analysis.dump()

    stop_simulation = False
    reactions_enabled = False
    energy0 = 0.0
    bonds0 = 0.0

    rate_file = None
    if args.rate_arrhenius:
        rate_file = open('{}_{}_new_rates.csv'.format(args.output_prefix, rng_seed), 'w')

    if args.skin == 'auto':
        print('Tunning skin parameter.')
        skin = espressopp.tools.decomp.tuneSkin(
            system, integrator, minSkin=0.1, maxSkin=1.5, precision=0.0001, printInfo=True)
        print('Found skin: {}'.format(skin))
        integrator.step = 0

    totalTime = time.time()
    integratorLoop = 0.0

    hook_before_sim(system, integrator, ar, gt)

    for k in range(sim_step):
        system_analysis.info()
        if save_traj_topology and k_trj_collect > 0 and k % k_trj_collect == 0:
            traj_file.dump(k * integrator_step, k * integrator_step * args.dt)
        if save_traj_topology and k_trj_flush > 0 and k % k_trj_flush == 0:
            dump_topol.update()
            traj_file.flush()  # Write HDF5 to disk.
        if k_enable_reactions == k:
            print('Enabling chemical reactions')
            integrator.addExtension(ar)
            if extensions_integrator:
                for ext in extensions_integrator:
                    integrator.addExtension(ext)
            reactions_enabled = True
            # Saves coordinate output file.
            output_gro_file = '{}_{}_before_reaction_confout.gro'.format(args.output_prefix, args.rng_seed)
            input_conf.update_position(system, unfolded=True)
            input_conf.write(output_gro_file, force=True)
            print('Save configuration before start of the reaction, filename: {}'.format(output_gro_file))
            if sc.exclusions_list:
                dynamic_exclusion_list.exclude(sc.exclusions_list)
                print('Add {} new exclusions from restrict reactions'.format(len(sc.exclusions_list)))

            if not hook_init_reaction(system, integrator, ar, gt, args):
                raise RuntimeError('hook_init_reaction return False')
            if not save_traj_topology:
                print('Enabling saving topology and trajectory')
                save_traj_topology = True
                ext_dump = espressopp.integrator.ExtAnalyze(dump_topol, args.topol_collect)
                integrator.addExtension(ext_dump)
                dump_topol.dump()
                dump_topol.update()

        if reactions_enabled:
            for obs, stop_value in maximum_conversion:
                val = obs.compute()
                if val >= stop_value:
                    print('Reaches {} of the conversion => Stop simulation'.format(val))
                    stop_simulation = True
            if stop_simulation:
                if eq_run == 0:
                    break
                else:
                    eq_run -= 1
            # Support for arrhenius law.
            if args.rate_arrhenius:
                bonds0 = sum(f.fpl.totalSize() for f in chem_fpls)  # TODO(jakub): this is terrible.
                energy0 = system_analysis.potential_energy

            if k_stop_reactions == k:
                ar.disconnect()

        loopTimer = time.time()
        integrator.run(integrator_step)
        integratorLoop += (time.time() - loopTimer)

        hook_at_step(system, integrator, ar, gt, args, k * integrator_step)

        if args.rate_arrhenius and reactions_enabled:
            bonds1 = sum(f.fpl.totalSize() for f in chem_fpls)  # TODO(jakub): this is terrible.
            delta_bonds = bonds1 - bonds0
            if delta_bonds > 0:
                energy_delta = (system_analysis.potential_energy - energy0) / float(delta_bonds)
                new_rate = math.exp(-energy_delta / temperature)
                print('{}\tChange reaction rate, delta_E={}, new_k={}, delta_bonds={}'.format(k * integrator_step,
                                                                                              energy_delta, new_rate,
                                                                                              delta_bonds))
                rate_file.write('{} {:e}\n'.format(k * integrator_step, new_rate))
                for r in reactions:
                    r.rate = new_rate
    totalTime = time.time() - totalTime
    ##### END of main integrator loop ###########

    system_analysis.info()
    traj_file.dump(sim_step * integrator_step, sim_step * integrator_step * args.dt)
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
        'total_steps': sim_step * integrator_step,
        'total_time': sim_step * integrator_step * args.dt,
        'integrator_step': integrator_step,
        'start_reaction': args.start_ar,
        'topology_collect': args.topol_collect,
        'trajectory_collect': trj_collect
    }
    for k, v in sim_params.items():
        g_params.attrs[k] = v
    tools.save_forcefield(h5, gt)
    h5.close()
    print('Closing HDF5 {} '.format(h5md_output_file))

    # Save topology
    max_pid = espressopp.analysis.MaxPID(system).compute()
    output_topol_filename = '{}_{}_output_topol.top'.format(args.output_prefix, args.rng_seed)
    out_topol = chemlab.files_io.GROMACSTopologyFile(output_topol_filename)
    out_topol.atomtypes = gt.topol.atomtypes
    out_topol.bondtypes = gt.topol.bondtypes
    out_topol.angletypes = gt.topol.angletypes
    out_topol.dihedraltypes = gt.topol.dihedraltypes
    out_topol.nonbond_params = gt.topol.nonbond_params
    out_topol.atomstate = gt.topol.atomstate
    out_topol.defaults = gt.topol.defaults
    out_topol.system_name = gt.topol.system_name
    out_topol.moleculetype = {'MOL': 3}
    out_topol.molecules = [('MOL', 1)]
    out_topol.content = None

    valid_type_ids = None
    if args.table_groups:
        valid_type_ids = map(gt.atomsym_atomtype.get, args.table_groups.split(','))

    for at_pid in xrange(1, int(max_pid) + 1):
        p = system.storage.getParticle(at_pid)
        if p:
            if valid_type_ids and p.type not in valid_type_ids:
                continue
            if at_pid in gt.atoms:
                at_data = gt.atoms[at_pid]
                mol_name = at_data['molecule_name']
                if mol_name not in out_topol.molecules_data:
                    out_topol.molecules_data[mol_name] = {'atoms': {}}

                topo_atom = out_topol.molecules_data[mol_name]['atoms'].setdefault(at_pid, chemlab.files_io.TopoAtom())
                topo_atom.atom_id = at_pid
                topo_atom.name = at_data['name']
                topo_atom.atom_type = gt.atomtype_atomsym[p.type]
                topo_atom.atom_type_id = p.type
                topo_atom.mass = p.mass
                topo_atom.chain_idx = p.res_id
                topo_atom.chain_name = at_data['chain_name']
                topo_atom.charge = p.q
                topo_atom.cgnr = at_pid
                out_topol.atoms[at_pid] = topo_atom
            else:
                if 'MOLX' not in out_topol.molecules_data:
                    out_topol.molecules_data['MOLX'] = {'atoms': {}}

                topo_atom = chemlab.files_io.TopoAtom(
                    atom_id=at_pid,
                    atom_type=gt.atomtype_atomsym[p.type],
                    chain_idx=p.res_id,
                    chain_name='CH{}'.format(p.res_id),
                    name='X{}'.format(p.type),
                    cgnr=at_pid,
                    charge=p.q,
                    mass=p.mass)
                topo_atom.atom_type_id = p.type
                out_topol.molecules_data['MOLX']['atoms'][at_pid] = topo_atom
                out_topol.atoms[at_pid] = out_topol.molecules_data['MOLX']['atoms'][at_pid]
                # Extend also input_conf
                input_conf.atoms[at_pid] = chemlab.files_io.Atom(
                    atom_id=at_pid,
                    name='X{}'.format(p.type),
                    chain_idx=p.res_id,
                    chain_name='C{}'.format(p.type),
                    position=numpy.zeros(3),
                    velocity=numpy.zeros(3))

    with open('{}_{}_bonds.dat'.format(args.output_prefix, args.rng_seed), 'w') as of:
        bond_lists = []
        for fpl in static_fpls:
            fpl_params = fpl.params
            for p in fpl.getAllBonds():
                bond_lists.append([p[0], p[1], fpl_params[0]] + list(fpl_params[1]) + ['; static'])
        for func, fpl in dynamic_fpls.items():
            for p in fpl.getAllBonds():
                t0, t1 = out_topol.atoms[p[0]].atom_type_id, out_topol.atoms[p[1]].atom_type_id
                namet0, namet1 = out_topol.atoms[p[0]].atom_type, out_topol.atoms[p[1]].atom_type
                fpl_params = fpl.params[t0][t1]
                if fpl_params:
                    bond_lists.append([p[0], p[1], fpl_params['func']] + fpl_params['params'] + ['; dynamic'])
                else:
                    bond_lists.append([p[0], p[1]] + ['; MISSING params type: {}-{} dynamic'.format(namet0, namet1)])
        for def_f in chem_fpls:
            for p in def_f.fpl.getAllBonds():
                t0, t1 = out_topol.atoms[p[0]].atom_type_id, out_topol.atoms[p[1]].atom_type_id
                namet0, namet1 = out_topol.atoms[p[0]].atom_type, out_topol.atoms[p[1]].atom_type
                params = None
                if namet0 in out_topol.bondtypes:
                    if namet1 in out_topol.bondtypes[namet0]:
                        params = out_topol.bondtypes[namet0][namet1]
                elif namet1 in out_topol.bondtypes:
                    if namet0 in out_topol.bondtypes[namet1]:
                        params = out_topol.bondtypes[namet1][namet0]
                if params:
                    bond_lists.append([p[0], p[1], '{} {} ; chem {}-{}'.format(
                        params['func'], ' '.join(params['params']), namet0, namet1)])
                else:
                    bond_lists.append([p[0], p[1], '; chem MISSING params type: {}-{}'.format(namet0, namet1)])

        for ptypes, fpl in registered_fpls:
            for b in fpl.getAllBonds():
                fpl_params = fpl.params
                bond_lists.append([b[0], b[1], fpl_params[0]] + list(fpl_params[1]) +
                                  ['; special tuple types: {}-{}'.format(ptypes[0], ptypes[1])])

        for b in bond_lists:
            of.write('{}\n'.format(' '.join(map(str, b))))
            out_topol.new_data['bonds'][(b[0], b[1])] = b[2:]

    with open('{}_{}_angles.dat'.format(args.output_prefix, args.rng_seed), 'w') as of:
        angle_lists = []
        for ftl in static_ftls:
            ftl_params = ftl.params
            for p in ftl.getAllTriples():
                angle_lists.append(list(p) + [ftl_params[0]] + list(ftl_params[1]) + ['; static'])
        for func, ftl in dynamic_ftls.items():
            for p in ftl.getAllTriples():
                t0, t1, t2 = (out_topol.atoms[p[0]].atom_type_id, out_topol.atoms[p[1]].atom_type_id,
                              out_topol.atoms[p[2]].atom_type_id)
                namet0, namet1, namet2 = (
                    out_topol.atoms[p[0]].atom_type,
                    out_topol.atoms[p[1]].atom_type,
                    out_topol.atoms[p[2]].atom_type)
                ftl_params = ftl.params[t0][t1][t2]
                if ftl_params:
                    angle_lists.append(list(p) + [ftl_params['func']] + list(ftl_params['params']) + ['; dynamic'])
                else:
                    angle_lists.append(
                        list(p) + ['; MISSING params type: {}-{}-{} dynamic'.format(namet0, namet1, namet2)])
        for a in angle_lists:
            of.write('{}\n'.format(' '.join(map(str, a))))
            out_topol.new_data['angles'][tuple(a[:3])] = a[3:]

    with open('{}_{}_dihedrals.dat'.format(args.output_prefix, args.rng_seed), 'w') as of:
        dih_lists = []
        for fql in static_fqls:
            fql_params = fql.params
            for p in fql.getAllQuadruples():
                dih_lists.append(list(p) + [fql_params[0]] + list(fql_params[1]) + ['; static'])
        for func, fql in dynamic_fqls.items():
            for p in fql.getAllQuadruples():
                t0, t1, t2, t3 = (out_topol.atoms[p[0]].atom_type_id, out_topol.atoms[p[1]].atom_type_id,
                                  out_topol.atoms[p[2]].atom_type_id, out_topol.atoms[p[3]].atom_type_id)
                namet0, namet1, namet2, namet3 = (
                    out_topol.atoms[p[0]].atom_type,
                    out_topol.atoms[p[1]].atom_type,
                    out_topol.atoms[p[2]].atom_type,
                    out_topol.atoms[p[3]].atom_type)
                fql_params = fql.params[t0][t1][t2][t3]
                if fql_params:
                    dih_lists.append(list(p) + [fql_params['func']] + list(fql_params['params']) + ['; dynamic'])
                else:
                    dih_lists.append(list(p) + ['; MISSING params type: {}-{}-{}-{} dynamic'.format(
                        namet0, namet1, namet2, namet3)])
        for d in dih_lists:
            of.write('{}\n'.format(' '.join(map(str, d))))
            out_topol.new_data['dihedrals'][tuple(d[:3])] = d[3:]
    print('Write output topology: {}'.format(output_topol_filename))
    out_topol.write(output_topol_filename)
    # End save tuples.

    with open('{}_{}_benchmark.csv'.format(args.output_prefix, args.rng_seed), 'a+') as benchmark_file:
        benchmark_file.write('{} {} {} {}\n'.format(MPI.COMM_WORLD.size, NPart, totalTime, integratorLoop))

    if args.rate_arrhenius:
        print('Changes in reaction rates written to {}'.format(rate_file.name))
        rate_file.close()

    topology_manager.save_topology('{}_{}_topology.dat'.format(args.output_prefix, args.rng_seed))
    topology_manager.save_res_topology('{}_{}_res_topology.dat'.format(args.output_prefix, args.rng_seed))
    topology_manager.save_residues('{}_{}_residue_list.dat'.format(args.output_prefix, args.rng_seed))

    # Saves coordinate output file.
    output_gro_file = '{}_{}_confout.gro'.format(args.output_prefix, rng_seed)
    input_conf.update_position(system, unfolded=False)
    input_conf.write(output_gro_file, force=True)
    print('Wrote end configuration to: {}'.format(output_gro_file))

    output_whole_gro = '{}_{}_whole_confout.gro'.format(args.output_prefix, rng_seed)
    dump_gro = espressopp.io.DumpGRO(system, integrator, filename=output_whole_gro)
    dump_gro.dump()
    print('Wrote whole configuration to: {}'.format(output_whole_gro))

    # Save fix distances, this is temporary
    if sc is not None and sc.fix_distances:
        all_fix_distances = []
        for fd in sc.fix_distances:
            all_fix_distances.extend(fd.get_all_triplets())
        with open('{}_all_fix_distances.pck'.format(global_file_prefix), 'wb') as out_fd:
            cPickle.dump(all_fix_distances, out_fd)

    if sc is not None:
        ar.save_reaction_counters('{}_reaction_counters'.format(global_file_prefix))
        with open('{}_reaction_counters'.format(global_file_prefix), 'a') as outf:
            outf.write('\n\nReaction index\n')
            for ridx in sorted(sc.reaction_index):
                outf.write('{} {}\n'.format(ridx, sc.reaction_index[ridx]))
        print('Saved reactioncounters: {}_reaction_counters'.format(global_file_prefix))

        ar.save_intra_inter_counter('{}_intra_inter_counters'.format(global_file_prefix))
        print('Saved intra/inter reactions counter: {}_intra_inter_counters'.format(global_file_prefix))

    total_time = time.time() - time0

    topol_timers = collections.defaultdict(list)
    for kv in topology_manager.get_timers():
        for k, v in kv:
            topol_timers[k].append(v)
    for k in topol_timers:
        if len(topol_timers[k]) > 0:
            topol_timers[k] = sum(topol_timers[k]) / float(len(topol_timers[k]))
    print('Topology manager timers:')
    for k, v in topol_timers.items():
        print('\t{}: {}'.format(k, v))

    print('DumpH5MD timers:')
    traj_timers = collections.defaultdict(list)
    for kv in traj_file.getTimers():
        for k, v in kv.items():
            traj_timers[k].append(v)
    for k in traj_timers:
        if len(traj_timers[k]) > 0:
            traj_timers[k] = sum(traj_timers[k]) / float(len(traj_timers[k]))
    for k, v in traj_timers.items():
        print('\t{}: {}'.format(k, v))

    with open('{}_{}_benchmark.pck'.format(args.output_prefix, args.rng_seed), 'wb') as benchmark_file:
        benchmark_data = {}
        benchmark_data['traj_timers'] = traj_timers
        benchmark_data['topol_timers'] = topol_timers
        benchmark_data['integrator_timers'] = tools.get_integrator_timers(integrator.getTimers(), system)
        benchmark_data['extension_timers'] = {}
        for ext_id in range(integrator.getNumberOfExtensions()):
            ext = integrator.getExtension(ext_id)
            ext_timers = tools.average_timers(ext.get_timers())
            if ext_timers:
                ext_name = str(ext)
                ext_name = '{}_{}'.format(re.match(r'\<(.*) object', ext_name).groups()[0], id(ext))
                benchmark_data['extension_timers'][ext_name] = tools.average_timers(ext.get_timers())
        benchmark_data['verlet_list'] = tools.average_timers(verletlist.get_timers())
        cPickle.dump(benchmark_data, benchmark_file)

    print('Final time analysis (per CPUs - {}) [s]:'.format(MPI.COMM_WORLD.size))
    espressopp.tools.analyse.final_info(system, integrator, verletlist, time0, time.time())

    print('Total time: {}'.format(total_time))
    print('Finished! Thanks!')


if __name__ == '__main__':
    main()
