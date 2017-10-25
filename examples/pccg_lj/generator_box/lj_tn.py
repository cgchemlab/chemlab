import espressopp
import argparse
import os
import mpi4py.MPI as MPI
import time

import polymer_melt


def _args():
    parser = argparse.ArgumentParser('Eq box for control chain growth')
    parser.add_argument('--conf', help='Input coordinates')
    parser.add_argument('--monomer_size', default=2, type=int)
    parser.add_argument('--nmonomers', default=2000, type=int)
    parser.add_argument('--nsolvents', default=11200, type=int)
    parser.add_argument('--density', default=0.74, type=float)

    return parser


def main():
    args = _args().parse_args()

    # Settings
    Npart = args.nmonomers * args.monomer_size + args.nsolvents
    CPUs = MPI.COMM_WORLD.size
    chain_size = args.monomer_size
    rho = args.density
    L = pow(Npart / rho, 1.0 / 3.0)
    box = (L, L, L)
    r_cutoff = pow(2, 1.0/6.0)
    skin = 0.3
    temperature = 0.5
    dt = 0.001
    epsilon = 1.0
    sigma = 1.0

    monomer_type = 0
    solvent_type = 1

    types = [monomer_type, solvent_type]
    matrix_types = [(types[i], types[j]) for i in range(len(types)) for j in range(i, len(types))]

    warmup_cutoff = pow(2.0, 1.0 / 6.0)
    warmup_nloops = 200
    warmup_isteps = 100
    total_warmup_steps = warmup_nloops * warmup_isteps
    epsilon_start = 0.001
    epsilon_end = 1.0
    epsilon_delta = (epsilon_end - epsilon_start) / warmup_nloops
    capradius = 0.6
    equil_nloops = 100
    equil_isteps = 100
    print(r_cutoff)

    if args.conf:
        conf = args.conf if os.path.exists(args.conf) else None
    else:
        conf = None

    system, integrator, vl = polymer_melt.PolymerMelt(
        args.nmonomers,
        args.nsolvents,
        chain_size,
        box=box,
        rc=r_cutoff,
        skin=skin,
        dt=dt,
        epsilon=epsilon,
        sigma=sigma,
        xyzfilename=conf,
        temperature=temperature,
        monomer_type=monomer_type,
        solvent_type=solvent_type)


    if conf is None:
        # First warm up so the defined potential should be removed
        system.removeInteractionByName('lj')

        LJpot = espressopp.interaction.LennardJonesCapped(
            epsilon=epsilon_start, sigma=sigma, cutoff=warmup_cutoff, caprad=capradius, shift='auto')
        interaction = espressopp.interaction.VerletListLennardJonesCapped(vl)
        for t1, t2 in matrix_types:
            interaction.setPotential(type1=t1, type2=t2, potential=LJpot)
        system.addInteraction(interaction, 'lj-cap')

        #cap_force = espressopp.integrator.CapForce(system, 1000)
        #integrator.addExtension(cap_force)

        espressopp.tools.analyse.info(system, integrator)
        print('warm-up')
        for step in range(warmup_nloops):
            integrator.run(warmup_isteps)
            LJpot.epsilon += epsilon_delta
            for t1, t2 in matrix_types:
                interaction.setPotential(type1=t1, type2=t2, potential=LJpot)
            espressopp.tools.analyse.info(system, integrator)

        system.removeInteractionByName('lj-cap')
        #cap_force.disconnect()

        interaction = espressopp.interaction.VerletListLennardJones(vl)
        ljpot = espressopp.interaction.LennardJones(
            epsilon=epsilon, sigma=sigma, cutoff=r_cutoff, shift='auto')
        for t1, t2 in matrix_types:
            interaction.setPotential(type1=t1, type2=t2, potential=ljpot)

        system.addInteraction(interaction, 'lj')
        # system.storage.cellAdjust()

        integrator.resetTimers()
        integrator.step = 0
        filename = 'lennard_jones_fluid_N_{}.xyz'.format(Npart)
        espressopp.tools.writexyz(
            filename, system, velocities=True, unfolded=True)

    dump_gro = espressopp.io.DumpGRO(system, integrator, 'conf.gro', unfolded=True, append=True)
    ext_dump_gro = espressopp.integrator.ExtAnalyze(dump_gro, 100)
    integrator.addExtension(ext_dump_gro)

    print('simulation')
    espressopp.tools.analyse.info(system, integrator)
    time0 = time.time()
    for step in range(equil_nloops):
        integrator.run(equil_isteps)
        espressopp.tools.analyse.info(system, integrator)
    dtime = time.time() - time0

    from md_libs import files_io
    gro = files_io.GROFile('confout.gro')
    gro.box = box
    pid2rid = {i: (i+1)/2 for i in range(1, args.nmonomers*chain_size+1)}
    pid2rid.update({i: i for i in range(args.nmonomers*chain_size+1, Npart+1)})
    pid2name = {i: 'A' for i in range(1, args.nmonomers*chain_size+1)}
    pid2name.update({i: 'B' for i in range(args.nmonomers*chain_size+1, Npart+1)})
    pid2res = {i: 'MON' for i in range(1, args.nmonomers*chain_size+1)}
    pid2res.update({i: 'SOL' for i in range(args.nmonomers*chain_size+1, Npart+1)})
    for pid in range(1, Npart+1):
        p = system.storage.getParticle(pid)
        gro.atoms[pid] = files_io.Atom(
            atom_id=pid,
            name=pid2name[pid],
            chain_name=pid2res[pid],
            chain_idx=pid2rid[pid],
            position=p.pos)
    gro.write('confout.gro', force=True)

if __name__ == '__main__':
    main()
