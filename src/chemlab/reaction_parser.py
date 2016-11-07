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
import espressopp
import ConfigParser
import random
import re
import warnings

__doc__ = """This is a reaction parser."""

# Constants
REACTION_NORMAL = 'normal'
REACTION_DISSOCATION = 'diss'
REACTION_EXCHANGE = 'exchange'
REACTION_RESTRICT = 'restricted'
EXT_POSTPROCESS = 'PP'
EXT_INTEGRATOR = 'Integrator'


def parse_equation(input_string):
    """Parse chemical equation and returns properties extracted from the input string

    The reaction is in the form some reactants -> some products
    where the reactants are splited by '+' and the same with products.

    Args:
        input_string: The input string.

    Returns:
        The dictionary with the extracted data.
    """
    re_reactant = re.compile(r'(?P<name>\w+)\((?P<min>\d+),\s*(?P<max>\d+)\)')
    re_product = re.compile(r'(?P<name>\w+)\((?P<delta>[0-9-]+)\)')

    reactants, products = input_string.split('->')

    reactant_list = {}

    mol_a, mol_b = map(re_reactant.match, map(str.strip, reactants.split('+')))

    reactant_list['type_1'] = mol_a.groupdict()
    reactant_list['type_2'] = mol_b.groupdict()

    products = [re_product.match(x).groupdict() for x in map(str.strip, products.split(':'))]
    reactant_list['type_1']['delta'] = products[0]['delta']
    reactant_list['type_2']['delta'] = products[1]['delta']
    reactant_list['type_1']['new_type'] = products[0]['name']
    reactant_list['type_2']['new_type'] = products[1]['name']

    return reactant_list, REACTION_NORMAL


def parse_reverse_equation(input_string):
    """Parse chemical reaction: dissociation

    Equation:
        A[min,max):B[min,max) -> A(deltaA) + B(deltaB)

    Args:
        input_string: The input string with equation to parse.
    """
    re_reactant = re.compile(r'(?P<name>\w+)\((?P<min>\d+),\s*(?P<max>\d+)\)')
    re_product = re.compile(r'(?P<name>\w+)\((?P<delta>[0-9-]+)\)')

    reactant_list = {}

    reactants, products = map(str.strip, input_string.split('->'))
    mol_a, mol_b = map(re_reactant.match, map(str.strip, reactants.split(':')))

    reactant_list['type_1'] = mol_a.groupdict()
    reactant_list['type_2'] = mol_b.groupdict()
    products = [re_product.match(x).groupdict() for x in map(str.strip, products.split('+'))]
    reactant_list['type_1']['delta'] = products[0]['delta']
    reactant_list['type_2']['delta'] = products[1]['delta']
    reactant_list['type_1']['new_type'] = products[0]['name']
    reactant_list['type_2']['new_type'] = products[1]['name']

    return reactant_list, REACTION_DISSOCATION


def parse_exchange_equation(input_string):
    """Parse chemical reaction: exchange

    Equation:
        A[min,max):B[min,max) + C[min,max) -> A(deltaA):C(deltaC) + B(deltaB)

    Args:
        input_string: The input string with equation to parse.
    """
    re_reactant = re.compile(r'(?P<name>\w+)\((?P<min>\d+),\s*(?P<max>\d+)\)')
    re_product = re.compile(r'(?P<new_type>\w+)\((?P<delta>[0-9-]+)\)')

    reactant_list = {}

    reactants, products = map(str.strip, input_string.split('->'))
    part_a, part_b = [map(str.strip, x.split(':')) for x in reactants.split('+')]
    mol_a, mol_b = [x.groupdict() for x in map(re_reactant.match, part_a)]
    mol_c = re_reactant.match(part_b[0]).groupdict()

    product_a, product_b = [map(str.strip, x.split(':')) for x in products.split('+')]
    prod_a, prod_b = [x.groupdict() for x in map(re_product.match, product_a)]
    prod_c = re_product.match(product_b[0]).groupdict()

    reactant_list['type_1'] = mol_a
    reactant_list['type_2'] = mol_b
    reactant_list['type_3'] = mol_c
    reactant_list['type_1'].update(prod_a)
    reactant_list['type_2'].update(prod_c)
    reactant_list['type_3'].update(prod_b)

    return reactant_list, REACTION_EXCHANGE


def process_reaction(reaction):
    """Process a single reaction section."""
    reaction = dict(reaction)

    group = reaction['group']
    data = {
        'rate': float(reaction['rate']),
        'intramolecular': eval(reaction.get('intramolecular', 'False')),
        'intraresidual': eval(reaction.get('intraresidual', 'False')),
        'virtual': eval(reaction.get('virtual', 'False')),
    }

    reaction_parsers = [parse_equation, parse_reverse_equation, parse_exchange_equation]
    reaction_type = None
    for reaction_parser in reaction_parsers:
        try:
            data['reactant_list'], reaction_type = reaction_parser(reaction['reaction'])
        except:
            continue

    if reaction_type is None:
        raise RuntimeError('Could not parse reaction equation: {}'.format(reaction['reaction']))
    print('Reaction_type: {}: {}'.format(reaction['reaction'], reaction_type))

    data['reaction_type'] = reaction_type

    if reaction_type == REACTION_NORMAL:
        if 'min_cutoff' in reaction:
            data['min_cutoff'] = float(reaction['min_cutoff'])

        # Support for smooth cut-off
        if 'sigma' in reaction and 'eq_distance' in reaction:
            data['sigma'] = float(reaction['sigma'])
            data['eq_distance'] = float(reaction['eq_distance'])
        elif 'cutoff' in reaction:
            data['cutoff'] = float(reaction['cutoff'])
        else:
            raise RuntimeError('Please define cutoff of the reaction')
    elif reaction_type == REACTION_DISSOCATION:
        if 'diss_rate' in reaction:
            if not data['reverse']:
                raise Exception('Invalid keyword `diss_rate` for non-dissociation reaction')
            data['diss_rate'] = float(reaction['diss_rate'])

    if 'active' in reaction:
        data['active'] = eval(reaction['active'])
    else:
        data['active'] = True

    if 'connectivity_map' in reaction:
        data['reaction_type'] = REACTION_RESTRICT

    return (group, data)


def process_general(cfg):
    """Process general section."""
    cfg = dict(cfg)
    if cfg.get('bond_limit'):
        warnings.warn('Bond limit not supported anymore!')

    return {
        'interval': int(cfg['interval']),
        'nearest': bool(cfg.get('nearest', False)),
        'pair_distances_filename': cfg.get('pair_distances_filename')
    }


def process_group(cfg):
    """Process group section."""
    cfg = dict(cfg)
    group_cfg = {
        'reaction_list': [],
        'connectivity_map': cfg.get('connectivity_map'),
        'extensions': {}
    }
    if 'extensions' in cfg:
        group_cfg['extensions'] = {s.strip(): None for s in cfg['extensions'].split(',')}

    if 'potential' in cfg:
        group_cfg['potential'] = cfg['potential']
        group_cfg['potential_options'] = dict([s.split('=') for s in cfg['potential_options'].split(',')])

    if 'eq_length' in cfg:
        group_cfg['eq_length'] = float(cfg['eq_length'])
        group_cfg['final_type'] = cfg['final_type']
        group_cfg['alpha'] = float(cfg['alpha'])

    return group_cfg


def process_extension(cfg):
    """Process extension entry."""
    cfg = dict(cfg)
    ret = {'class': cfg['ext_type']}
    del cfg['ext_type']
    ret['options'] = cfg
    return ret


def parse_config(input_file):
    parser = ConfigParser.SafeConfigParser()
    parser.read(input_file)

    config = {'general': None, 'reactions': {}}
    extensions = {}

    for s in parser.sections():
        if s == 'general':
            config['general'] = process_general(parser.items(s))
        elif s.startswith('group_'):
            group_name = s.replace('group_', '').strip()
            group_opt = None
            if group_name not in config['reactions']:
                group_opt = process_group(parser.items(s))
                config['reactions'][group_name] = group_opt
                for ext in group_opt['extensions']:
                    group_opt['extensions'][ext] = extensions[ext]
        elif s.startswith('ext_'):
            name = s.replace('ext_', '').strip()
            properties = process_extension(parser.items(s))
            if name in extensions:
                raise RuntimeError('Name of extension already exists')
            extensions[name] = properties
        elif s.startswith('reaction_'):
            group_name, data = process_reaction(parser.items(s))
            if group_name not in config['reactions']:
                raise RuntimeError(
                    'Wrong order, first reaction groups and then referring reactions')
            config['reactions'][group_name]['reaction_list'].append(data)
    return config
