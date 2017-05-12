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


import collections
import numpy as np
import re
import espressopp

__doc__ = "Tool functions."


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
