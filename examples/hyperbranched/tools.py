"""
Copyright (C) 2015-2016 Jakub Krajniak <jkrajniak@gmail.com>

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
import ast
import numpy as np
import random

__doc__ = "Tool functions."


class MyArgParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(MyArgParser, self).__init__(*args, **kwargs)

    def convert_arg_line_to_args(self, line):
        for arg in line.split():
            t = arg.strip()
            if not t:
                continue
            if t.startswith('#'):
                break
            if not t.startswith('--'):
                t = '--{}'.format(t)
            yield t

    @staticmethod
    def save_to_file(output_file, namespace):
        """Saves arguments to file so it can be read again.

        Args:
            output_file: The string with the name of output file.
            namespace: The namespace with arguements.
        """
        with open(output_file, "w") as of:
            keys = sorted(namespace.__dict__.keys())
            for k in keys:
                v = namespace.__dict__[k]
                if v is not None:
                    of.write('{}={}\n'.format(k, v))


def dump_topol(file_name, topol, system, particle_ids, bonds, angles, dihedrals, pairs):
    # Get current atom set.
    for atid in particle_ids:
        p = system.storage.getParticle(atid)
        atom_type_params = topol.gt.atomtypes[topol.atomtype_atomsym[p.type]]
        topo_atom = espressopp.tools.chemlab.files_io.TopoAtom()
        topo_atom.atom_id = atid
        topo_atom.atom_type = atom_type_params['name']
        topo_atom.chain_name = topol.gt.molecules.keys()[0]
        topo_atom.name = 'T{}'.format(p.type)
        topo_atom.mass = atom_type_params['mass']
        topo_atom.charge = atom_type_params['charge']
        topo_atom.chain_idx = p.res_id
        topol.topol.atoms[atid] = topo_atom

    for fpl in bonds:
        for bp in fpl.getBonds():
            for b12 in bp:
                topol.topol.new_data['bonds'][b12] = []

    for ftl in angles:
        for bt in ftl.getTriples():
            for b123 in bt:
                topol.topol.angles[b123] = []

    for fql in dihedrals:
        for bq in fql.getQuadruples():
            for b1234 in bq:
                topol.topol.dihedrals[b1234] = []

    for fpr in pairs:
        for bl in fpr.getBonds():
            for b12 in bl:
                topol.topol.pairs[b12] = []

    topol.topol.write(file_name)

def save_forcefield(h5, gt):
    """Saves force-field to H5MD file under the /parameters/forcefield group."""
    if 'force_field' not in h5['/parameters']:
        g_ff = h5['/parameters'].create_group('force_field')
    g_ff = h5['/parameters/force_field']
    atomtypes = []
    for at_sym, atd in gt.topol.atomtypes.items():
        atomtypes.append([
            gt.atomsym_atomtype[at_sym], at_sym, atd['mass'],
            atd['charge'], atd['epsilon'], atd['sigma'],
            atd['type']])
    atomtypes = np.array(atomtypes)
    d_atomtypes = g_ff.create_dataset('atomtypes', data=atomtypes)
    d_atomtypes.attrs['keys'] = [
        'type_id', 'name', 'mass', 'charge', 'epsilon', 'sigma', 'type']


def _args():
    parser = MyArgParser(description='Runs classical MD simulation',
                         fromfile_prefix_chars='@')
    parser.add_argument('--conf', required=True, help='Input .gro coordinate file')
    parser.add_argument('--top', '--topology', required=True, help='Topology file',
                        dest='top')
    parser.add_argument('--node_grid')
    parser.add_argument('--skin', type=float, default=0.16,
                        help='Skin value for Verlet list')
    parser.add_argument('--run', type=int, default=10000,
                        help='Number of simulation steps')
    parser.add_argument('--int_step', default=1000, type=int, help='Steps in integrator')
    parser.add_argument('--rng_seed', type=int, help='Seed for RNG', required=False,
                        default=random.randint(1000, 10000))
    parser.add_argument('--output_prefix',
                        default='sim', type=str,
                        help='Prefix for output files')
    parser.add_argument('--output_file',
                        default='trjout.h5', type=str,
                        help='Name of output trajectory file')
    parser.add_argument('--thermostat',
                        default='lv',
                        choices=('lv', 'vr'),
                        help='Thermostat to use, lv: Langevine, vr: Stochastic velocity rescale')
    parser.add_argument('--barostat', default='lv', choices=('lv', 'br'),
                        help='Barostat to use, lv: Langevine, br: Berendsen')
    parser.add_argument('--barostat_tau', default=5.0, type=float,
                        help='Tau parameter for Berendsen barostat')
    parser.add_argument('--barostat_mass', default=50.0, type=float,
                        help='Mass parameter for Langevin barostat')
    parser.add_argument('--barostat_gammaP', default=1.0, type=float,
                        help='gammaP parameter for Langevin barostat')
    parser.add_argument('--thermostat_gamma', type=float, default=5.0,
                        help='Thermostat coupling constant')
    parser.add_argument('--temperature', default=458.0, type=float, help='Temperature')
    parser.add_argument('--pressure', help='Pressure', type=float)
    parser.add_argument('--trj_collect', default=1000, type=int,
                        help='Collect trajectory every (step)')
    parser.add_argument('--energy_collect', default=1000, type=int,
                        help='Collect energy every (step)')
    parser.add_argument('--topol_collect', default=1000, type=int,
                        help='Collect topology every (step)')
    parser.add_argument('--dt', default=0.001, type=float,
                        help='Integrator time step')
    parser.add_argument('--lj_cutoff', default=1.2, type=float,
                        help='Cutoff of atomistic non-bonded interactions')
    parser.add_argument('--cg_cutoff', default=1.4, type=float,
                        help='Cuoff of coarse-grained non-bonded interactions')
    parser.add_argument('--coulomb_epsilon1', default=1.0, type=float,
                        help='Epsilon_1 for coulomb interactions')
    parser.add_argument('--coulomb_epsilon2', default=80.0, type=float,
                        help='Epsilon_2 for coulomb interactions')
    parser.add_argument('--coulomb_kappa', default=0.0, type=float,
                        help='Kappa paramter for coulomb interactions')
    parser.add_argument('--coulomb_cutoff', default=0.9, type=float,
                        help='Coulomb cut-off')
    parser.add_argument('--reactions', default='reaction.cfg',
                        help='Configuration file with chemical reactions')
    parser.add_argument('--debug', default=None, help='Turn on logging mechanism')
    parser.add_argument('--start_ar', default=0, type=int, help='When to start chemical reactions')
    parser.add_argument('--store_species', default=False, type=ast.literal_eval,
                        help='Store particle types')
    parser.add_argument('--store_state', default=True, type=ast.literal_eval,
                        help='Store chemical state')
    parser.add_argument('--store_lambda', default=False, type=ast.literal_eval,
                        help='Store lambda parameter')

    return parser
