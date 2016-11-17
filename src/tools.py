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


def save_forcefield(h5, gt):
    """Saves force-field to H5MD file under the /parameters/forcefield group."""
    if 'force_field' not in h5['/parameters']:
        h5['/parameters'].create_group('force_field')
    g_ff = h5['/parameters/force_field']
    atomtypes = []
    for at_sym in gt.used_atomtypes:
        atd = gt.topol.atomtypes[at_sym]
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
    general_options = parser.add_argument_group('General options')
    general_options.add_argument('--conf', required=True, help='Input .gro coordinate file')
    general_options.add_argument('--top', '--topology', required=True, help='Topology file',
                        dest='top')
    general_options.add_argument('--node_grid')
    general_options.add_argument('--skin', default=0.16, help='Skin value for Verlet list')
    general_options.add_argument('--output_prefix',
                        default='sim', type=str,
                        help='Prefix for output files')
    general_options.add_argument('--output_file',
                        default='trjout.h5', type=str,
                        help='Name of output trajectory file')
    general_options.add_argument('--trj_collect', default=1000, type=int,
                        help='Collect trajectory every (step)')
    general_options.add_argument('--energy_collect', default=1000, type=int,
                        help='Collect energy every (step)')
    general_options.add_argument('--topol_collect', default=1000, type=int,
                        help='Collect topology every (step)')
    general_options.add_argument('--reactions', default=None,
                        help='Configuration file with chemical reactions')
    general_options.add_argument('--debug', default=None, help='Turn on logging mechanism')
    general_options.add_argument('--start_ar', default=0, type=int, help='When to start chemical reactions')
    general_options.add_argument('--stop_ar', default=-1, type=int, help='When to stop chemical reactions')
    general_options.add_argument('--maximum_conversion',
                        default=None,
                        help=('The comma separated list of conditions on which '
                              'the simulation will stop. (format: atom type symbol:max number:total number)'))
    general_options.add_argument('--eq_steps', default=0, help=('Run simulation after conversion reached for n-steps'), type=int)
    general_options.add_argument('--table_groups', default=None,
                        help='The list of atom type names that should be simulated with tabulated potential.')
    general_options.add_argument('--max_force', default=-1, type=float,
                        help='Maximum force in the system.')
    general_options.add_argument('--rate_arrhenius', default=False, help='Change rate based on the Arrhenius equation.',
                        type=ast.literal_eval)
    general_options.add_argument('--exclusion_list', default=None, help='Read exclusion list from external file')
    general_options.add_argument('--benchmark_data', default=None, help='Store time measurement in the file')

    args_simulation_options = parser.add_argument_group('Simulation parameters')
    args_simulation_options.add_argument('--kb', type=float, default=0.0083144621, help='Boltzmann constant, default kJ/mol')
    args_simulation_options.add_argument('--mass_factor', type=float, default=1.6605402, help='Mass units, default a.u.')
    args_simulation_options.add_argument('--run', type=int, default=10000,
                        help='Number of simulation steps')
    args_simulation_options.add_argument('--int_step', default=1000, type=int, help='Steps in integrator')
    args_simulation_options.add_argument('--rng_seed', type=int, help='Seed for RNG', required=False,
                        default=random.randint(1000, 10000))
    args_simulation_options.add_argument('--thermostat',
                        default='lv',
                        choices=('lv', 'vr', 'iso'),
                        help='Thermostat to use, lv: Langevine, vr: Stochastic velocity rescale')
    args_simulation_options.add_argument('--barostat', default='lv', choices=('lv', 'br'),
                        help='Barostat to use, lv: Langevine, br: Berendsen')
    args_simulation_options.add_argument('--barostat_tau', default=5.0, type=float,
                        help='Tau parameter for Berendsen barostat')
    args_simulation_options.add_argument('--barostat_mass', default=50.0, type=float,
                        help='Mass parameter for Langevin barostat')
    args_simulation_options.add_argument('--barostat_gammaP', default=1.0, type=float,
                        help='gammaP parameter for Langevin barostat')
    args_simulation_options.add_argument('--thermostat_gamma', type=float, default=5.0,
                        help='Thermostat coupling constant')
    args_simulation_options.add_argument('--temperature', default=458.0, type=float, help='Temperature')
    args_simulation_options.add_argument('--pressure', help='Pressure', type=float)
    args_simulation_options.add_argument('--dt', default=0.001, type=float,
                        help='Integrator time step')
    args_simulation_options.add_argument('--lj_cutoff', default=1.2, type=float,
                        help='Cutoff of atomistic non-bonded interactions')
    args_simulation_options.add_argument('--cg_cutoff', default=1.4, type=float,
                        help='Cuoff of coarse-grained non-bonded interactions')
    args_simulation_options.add_argument('--coulomb_epsilon1', default=1.0, type=float,
                        help='Epsilon_1 for coulomb interactions')
    args_simulation_options.add_argument('--coulomb_epsilon2', default=80.0, type=float,
                        help='Epsilon_2 for coulomb interactions')
    args_simulation_options.add_argument('--coulomb_kappa', default=0.0, type=float,
                        help='Kappa paramter for coulomb interactions')
    args_simulation_options.add_argument('--coulomb_cutoff', default=0.9, type=float,
                        help='Coulomb cut-off')

    args_storing_options = parser.add_argument_group('H5MD storage')
    args_storing_options.add_argument('--store_species', default=True, type=ast.literal_eval,
                        help='Store particle types')
    args_storing_options.add_argument('--store_state', default=True, type=ast.literal_eval,
                        help='Store chemical state')
    args_storing_options.add_argument('--store_position', default=True, type=ast.literal_eval,
                        help='Store positions')
    args_storing_options.add_argument('--store_lambda', default=False, type=ast.literal_eval,
                        help='Store lambda parameter')
    args_storing_options.add_argument('--store_force', default=False, type=ast.literal_eval,
                        help='Store forces')
    args_storing_options.add_argument('--store_velocity', default=False, type=ast.literal_eval,
                        help='Store velocity')
    args_storing_options.add_argument('--store_charge', default=False, type=ast.literal_eval,
                        help='Store charge')
    args_storing_options.add_argument('--store_mass', default=True, type=ast.literal_eval, help='Store mass')
    args_storing_options.add_argument('--store_res_id', default=True, type=ast.literal_eval, help='Store res_id')
    args_storing_options.add_argument('--store_pressure', default=False, type=ast.literal_eval,
                        help='Compute and store pressure')
    args_storing_options.add_argument('--store_single_precision', default=True, type=ast.literal_eval,
                        help='Write data in single precision format')

    args_counters = parser.add_argument_group('Counters')
    args_counters.add_argument('--count_types', default=None, help='List of particle types to count; eq. A,B')
    args_counters.add_argument('--count_tuples', default=False, type=ast.literal_eval, help='Count tuples')
    args_counters.add_argument('--count_types_state', default=None, help='List of particle types, state; eq. A:3,B:4')

    args_hybrid_bonds = parser.add_argument_group('Hybrid bonds')
    args_hybrid_bonds.add_argument('--t_hybrid_bond', default=0, type=int)
    args_hybrid_bonds.add_argument('--t_hybrid_angle', default=0, type=int)
    args_hybrid_bonds.add_argument('--t_hybrid_dihedral', default=0, type=int)

    return parser
