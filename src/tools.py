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
import collections
import logging
import numpy as np
import random
import re
import espressopp

__doc__ = "Tool functions."


class MyArgParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(MyArgParser, self).__init__(*args, **kwargs)

    def convert_arg_line_to_args(self, line):
        for arg in line.split():
            t = arg.strip()
            if not t:
                continue
            if t.startswith('#') or t.startswith(';'):
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


class RegexpFilter(logging.Filter):
    """This filter allows only messages that matches regexp to pass to user."""

    def __init__(self, regexp, name=''):
        self.regexp = re.compile(regexp)
        super(RegexpFilter, self).__init__(name)

    def filter(self, record):
        return self.regexp.match(record.msg) or self.regexp.match(record.funcName)


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
    general_options.add_argument('--trj_flush', default=None, type=int,
                                 help='Flush trajectory to storage')
    general_options.add_argument('--gro_trj_collect', default=None, type=int,
                                 help='Save trajectory in .gro file')
    general_options.add_argument('--energy_collect', default=1000, type=int,
                                 help='Collect energy every (step)')
    general_options.add_argument('--topol_collect', default=1000, type=int,
                                 help='Collect topology every (step)')
    general_options.add_argument('--reactions', default=None,
                                 help='Configuration file with chemical reactions')
    general_options.add_argument('--debug', default=None, help='Turn on logging mechanism')
    general_options.add_argument(
        '--check_topology',
        default=False,
        type=ast.literal_eval,
        help='If set to true then show the message generated by TopologyManager')
    general_options.add_argument('--start_ar', default=0, type=int, help='When to start chemical reactions')
    general_options.add_argument('--stop_ar', default=-1, type=int, help='When to stop chemical reactions')

    general_options.add_argument('--table_groups', default=None,
                                 help='The list of atom type names that should be simulated with tabulated potential.')
    general_options.add_argument('--thermal_groups', default=None,
                                 help='The list of atom type names that should be thermalized.')
    general_options.add_argument('--max_force', default=-1, type=float,
                                 help='Maximum force in the system.')
    general_options.add_argument('--rate_arrhenius', default=False, help='Change rate based on the Arrhenius equation.',
                                 type=ast.literal_eval)
    general_options.add_argument('--exclusion_list', default=None, help='Read exclusion list from external file')
    general_options.add_argument('--benchmark_data', default=None, help='Store time measurement in the file')
    general_options.add_argument('--system_monitor_filter', default=None,
                                 help='Print all (empty) or only selected elements in SystemMonitor')
    general_options.add_argument('--gen_velocity', default=True, type=ast.literal_eval,
                                 help='Generate velocity?')

    args_simulation_options = parser.add_argument_group('Simulation parameters')
    args_simulation_options.add_argument('--kb', type=float, default=0.0083144621,
                                         help='Boltzmann constant, default kJ/mol')
    args_simulation_options.add_argument('--mass_factor', type=float, default=1.6605402,
                                         help='Mass units, default a.u.')
    args_simulation_options.add_argument('--run', type=int, default=10000,
                                         help='Number of simulation steps')
    args_simulation_options.add_argument('--int_step', default=1000, type=int, help='Steps in integrator')
    args_simulation_options.add_argument('--rng_seed', type=int, help='Seed for RNG', required=False,
                                         default=random.randint(1000, 10000))
    args_simulation_options.add_argument('--thermostat',
                                         default='lv',
                                         choices=('lv', 'vr', 'iso', 'no'),
                                         help=('Thermostat to use, '
                                               'lv: Langevine, vr: Stochastic velocity rescale '
                                               'no: no thermostat'))
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
    args_storing_options.add_argument('--store_angdih', default=False, type=ast.literal_eval,
                                      help='Store angles, dihedrals')

    maximum_conversion_options = parser.add_argument_group('Maximum conversion')
    maximum_conversion_options.add_argument('--maximum_conversion',
                                            default=None,
                                            help=('The comma separated list of conditions on which '
                                                  'the simulation will stop. (format: '
                                                  'atom type symbol(state):max number:total number)'))
    maximum_conversion_options.add_argument('--eq_steps', default=0,
                                            help=('Run simulation after conversion reached for n-steps'), type=int)
    maximum_conversion_options.add_argument('--keep_simulation', default=False,
                                            help='Keep simulation running until the conversion is reached')

    args_counters = parser.add_argument_group('Counters')
    args_counters.add_argument('--count_types', default=None, help='List of particle types to count; eq. A,B')
    args_counters.add_argument('--count_tuples', default=False, type=ast.literal_eval, help='Count tuples')
    args_counters.add_argument('--count_types_state', default=None, help='List of particle types, state; eq. A:3,B:4')
    args_counters.add_argument('--count_fix_distances', default=False, type=ast.literal_eval, help='Count size of fix distances')

    args_hybrid_bonds = parser.add_argument_group('Hybrid bonds')
    args_hybrid_bonds.add_argument('--t_hybrid_bond', default=0, type=int)
    args_hybrid_bonds.add_argument('--t_hybrid_angle', default=0, type=int)
    args_hybrid_bonds.add_argument('--t_hybrid_dihedral', default=0, type=int)

    return parser


def get_integrator_timers(alltimers, system):
    skip_timers = ['timeRun']
    nprocs = len(alltimers)
    timers = {k: 0.0 for k, _ in alltimers[0]}
    alltimers = [{k: float(v) for k, v in ntimer if k not in skip_timers} for ntimer in alltimers]
    for ntimer in alltimers:
        for k, v in ntimer.items():
            timers[k] += v
    total_t = 0.0
    for k in timers:
        total_t += timers[k]
        timers[k] /= nprocs

    for k, v in sorted(timers.items()):
        if k.startswith('f'):
            lbl = system.getNameOfInteraction(int(k.replace('f', '')))
            timers[lbl] = v
            del timers[k]

    return timers


def average_timers(timer_list):
    avg_timers = collections.defaultdict(list)
    for cpu_list in timer_list:
        for k, v in cpu_list:
            avg_timers[k].append(v)
    for k in avg_timers:
        if len(avg_timers[k]) > 0:
            avg_timers[k] = sum(avg_timers[k]) / float(len(avg_timers[k]))

    return avg_timers

def get_maximum_conversion(args, system, chem_fpls, gt, cr_observs=None):
    if cr_observs is None:
        cr_observs = {}

    maximum_conversion = []

    re_type_state = re.compile(r'(?P<type>[A-Za-z0-9-]+)\(?(?P<state>\d?)\)?')
    re_max_conversion = re.compile(r'(?P<type>[A-Za-z0-9-]+)\(?(?P<state>\d?)\)?:(?P<maxnum>\d+):(?P<tot>\d+)')
    for o in args.maximum_conversion.split(','):
        type_symbols, max_number, tot_number = o.split(':')
        max_number = int(max_number)
        tot_number = int(tot_number)

        if '-' in type_symbols:
            vals = re_type_state.match(type_symbols)
            vals = vals.groupdict()
            type_symbol = vals['type']
            type_sym_1, type_sym_2 = type_symbol.split('-')
            type_id_1 = gt.used_atomsym_atomtype[type_sym_1]
            type_id_2 = gt.used_atomsym_atomtype[type_sym_2]
            for fpl_def in chem_fpls:
                if ((type_id_1, type_id_2) in fpl_def.type_list or
                        (type_id_2, type_id_1) in fpl_def.type_list):
                    obs_fpl = espressopp.analysis.NFixedPairListEntries(system, fpl_def.fpl)
                    maximum_conversion.append((obs_fpl, max_number))
                    break
        elif '+' in type_symbols:
            type_symbols = type_symbols.split('+')
            obs = espressopp.analysis.ChemicalConversionTypeState(system, total_count=tot_number)
            stop_value = float(max_number) / tot_number
            cr_types = []
            for type_symbol in type_symbols:
                vals = re_type_state.match(type_symbol)
                vals = vals.groupdict()
                if vals['state']:
                    type_state = int(vals['state'])
                else:
                    type_state = None
                type_name = vals['type']
                type_id_symbol = gt.used_atomsym_atomtype[type_name]
                obs.count_type(type_id_symbol, type_state)
                cr_types.append(type_name)
            cr_observs[(tuple(cr_types), tot_number, None)] = obs
            maximum_conversion.append((obs, stop_value))
        else:
            vals = re_type_state.match(type_symbols)
            vals = vals.groupdict()
            type_symbol = vals['type']
            if vals['state']:
                type_state = int(vals['state'])
            else:
                type_state = None
            stop_value = float(max_number) / tot_number
            type_id_symbol = gt.used_atomsym_atomtype[type_symbol]
            if (type_id_symbol, tot_number, type_state) not in cr_observs:
                if type_state is None:
                    obs = espressopp.analysis.ChemicalConversion(system, type_id_symbol, tot_number)
                else:
                    obs = espressopp.analysis.ChemicalConversionTypeState(system, type_id_symbol, type_state, tot_number)
                cr_observs[(type_id_symbol, tot_number, type_state)] = obs
            else:
                obs = cr_observs[(type_id_symbol, tot_number, type_state)]
            maximum_conversion.append((obs, stop_value))

    return maximum_conversion
