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

import espressopp
import ConfigParser
import random
import re
import warnings

__doc__ = """This is a reaction parser."""


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

    return reactant_list


def parse_reverse_equation(input_string):
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

    return reactant_list


def process_reaction(reaction):
    """Process a single reaction section."""
    reaction = dict(reaction)

    group = reaction['group']
    data = {
        'rate': float(reaction['rate']),
        'intramolecular': eval(reaction.get('intramolecular', 'False')),
        'intraresidual': eval(reaction.get('intraresidual', 'False')),
        'virtual': eval(reaction.get('virtual', 'False'))
    }

    try:
        data['reactant_list'] = parse_equation(reaction['reaction'])
        data['reverse'] = False
    except:
        data['reactant_list'] = parse_reverse_equation(reaction['reaction'])
        data['reverse'] = True

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

    if 'diss_rate' in reaction:
        if not data['reverse']:
            raise Exception('Invalid keyword `diss_rate` for non-dissociation reaction')
        data['diss_rate'] = float(reaction['diss_rate'])

    if 'active' in reaction:
        data['active'] = eval(reaction['active'])
    else:
        data['active'] = True

    return (group, data)


def process_general(cfg):
    """Process general section."""
    cfg = dict(cfg)
    if cfg.get('bond_limit'):
        warnings.warn('Bond limit not supported anymore!')

    return {
        'interval': int(cfg['interval']),
        'nearest': bool(cfg.get('nearest', False)),
    }


def process_group(cfg):
    """Process group section."""
    cfg = dict(cfg)
    return_gr = {'potential': cfg['potential'],
                 'potential_options': dict(
                     [s.split('=') for s in cfg['potential_options'].split(',')]),
                 'reaction_list': [],
                 'connectivity_map': cfg.get('connectivity_map'),
                 'extensions': {}
                 }
    if 'extensions' in cfg:
        return_gr['extensions'] = {s.strip(): None for s in cfg['extensions'].split(',')}
    return return_gr


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


class SetupReactions:
    """Main class to setup the reactions

    Args:
        system: The espressopp.System object.
        vl: The VerletList that is used to find the reacted pairs.
        topology_manager: The espressopp.integrator.TopologyManager object.
        config: The config file.
    """

    def __init__(self, system, vl, topol, topol_manager, config):
        self.system = system
        self.vl = vl
        self.topol = topol
        self.tm = topol_manager
        self.cfg = config
        self.name2type = topol.atomsym_atomtype
        self.dynamic_types = set()  # Stores the particle types that will change during the reactions.

        self.use_thermal_group = False
        self.fix_distance = None
        self.cr_observs = None  # Observs conversion types.

        self.exclusions_list = [] # For restric reactions, the exclusion lists has to be extended.

    def _setup_reaction(self, chem_reaction, fpl):
        """Setup single reaction.

        Args:
            chem_reaction: Dictionary with definition of the reactions
            fpl: The FixedPairList.

        Returns:
            The espressopp.integrator.Reaction object.
        """
        rl = chem_reaction['reactant_list']
        if not chem_reaction['active']:
            return None
        # Select reaction class.
        if chem_reaction['reverse']:
            r_class = espressopp.integrator.DissociationReaction
        elif chem_reaction.get('connectivity_map'):
            r_class = espressopp.integrator.RestrictReaction
        else:
            r_class = espressopp.integrator.Reaction

        rt1 = rl['type_1']['name']
        rt2 = rl['type_2']['name']
        r = r_class(
            type_1=self.name2type[rl['type_1']['name']],
            type_2=self.name2type[rl['type_2']['name']],
            delta_1=int(rl['type_1']['delta']),
            delta_2=int(rl['type_2']['delta']),
            min_state_1=int(rl['type_1']['min']),
            max_state_1=int(rl['type_1']['max']),
            min_state_2=int(rl['type_2']['min']),
            max_state_2=int(rl['type_2']['max']),
            rate=float(chem_reaction['rate']),
            fpl=fpl,
            cutoff=float(chem_reaction.get('cutoff', 0.0)))

        self.dynamic_types.add(self.name2type[rl['type_1']['name']])
        self.dynamic_types.add(self.name2type[rl['type_2']['name']])

        print('Setup reaction: {}({})-{}({})'.format(
            rt1, self.name2type[rt1], rt2, self.name2type[rt2]))
        if not chem_reaction['reverse']:
            r.intramolecular = bool(chem_reaction['intramolecular'])
            r.intraresidual = bool(chem_reaction['intraresidual'])
            r.is_virtual = bool(chem_reaction['virtual'])

            if 'sigma' in chem_reaction:
                r.set_reaction_cutoff(espressopp.integrator.ReactionCutoffRandom(
                    chem_reaction['eq_distance'], chem_reaction['sigma'], seed=random.randint(100, 100000)))

        if 'min_cutoff' in chem_reaction:
            r.get_reaction_cutoff().min_cutoff = float(chem_reaction['min_cutoff'])
        if 'diss_rate' in chem_reaction:
            r.diss_rate = float(chem_reaction['diss_rate'])
        if 'active' in chem_reaction:
            r.active = chem_reaction['active']

        if chem_reaction.get('connectivity_map'):
            print('Reading connectivity map {}, reaction will be restricted'.format(
                chem_reaction['connectivity_map']))
            connectivity_map = open(chem_reaction['connectivity_map'])
            ex_list = set()
            for l in connectivity_map.readlines():
                b1, b2 = map(int, l.strip().split())
                ex_list.add(tuple(sorted(b1, b2)))
            for b1, b2 in ex_list:
                r.define_connection(b1, b2)
            self.exclusions_list.extend(list(ex_list))
            print('Restricted to {} connections'.format(len(ex_list)))

        # Change type if necessary.
        if (rl['type_1']['name'] != rl['type_1']['new_type'] or
                    rl['type_2']['name'] != rl['type_2']['new_type']):
            r_pp = espressopp.integrator.PostProcessChangeProperty()
            t1_old = self.name2type[rl['type_1']['name']]
            t1_new = self.name2type[rl['type_1']['new_type']]
            if t1_old != t1_new:
                self.dynamic_types.add(t1_old)
                self.dynamic_types.add(t1_new)
                print('Reaction: {}-{}, change type {}->{}'.format(rt1, rt2, t1_old, t1_new))
                new_property = self.topol.gt.atomtypes[rl['type_1']['new_type']]
                r_pp.add_change_property(
                    t1_old,
                    espressopp.ParticleProperties(
                        t1_new, new_property['mass'],
                        new_property['charge']))

            t2_old = self.name2type[rl['type_2']['name']]
            t2_new = self.name2type[rl['type_2']['new_type']]
            if t2_old != t2_new:
                self.dynamic_types.add(t2_old)
                self.dynamic_types.add(t2_new)
                print('Reaction: {}-{}, change type {}->{}'.format(rt1, rt2, t2_old, t2_new))
                new_property = self.topol.gt.atomtypes[rl['type_2']['new_type']]
                r_pp.add_change_property(
                    t2_old,
                    espressopp.ParticleProperties(
                        t2_new, new_property['mass'],
                        new_property['charge']))

            r.add_postprocess(r_pp)
        return r

    def _prepare_group_postprocess(self, cfg):
        """Prepare extensions for the reactions."""
        pps = []

        # Belows are the sub-functions to create PostProcess objects.
        def _cfg_post_process_change_neighbour(cfg):
            """Setup PostProcessChangeNeighbourProperties"""
            pp = espressopp.integrator.PostProcessChangeNeighboursProperty(self.tm)
            type_transfers = [
                x.split('->') for x in cfg['type_transfers'].split(',')]
            for old_type, new_type in type_transfers:
                old_type, nb_level = old_type.split(':')
                nb_level = int(nb_level)
                if old_type != new_type:
                    print('Change property {}->{} nb={} and {}'.format(
                        old_type, new_type, nb_level, nb_level + 1))
                    t1_old = self.name2type[old_type]
                    t1_new = self.name2type[new_type]
                    self.dynamic_types.add(t1_old)
                    self.dynamic_types.add(t1_new)
                    new_property = self.topol.gt.atomtypes[new_type]
                    pp.add_change_property(
                        t1_old,
                        espressopp.ParticleProperties(
                            t1_new, new_property['mass'], new_property['charge']),
                        nb_level
                    )
                    pp.add_change_property(
                        t1_old,
                        espressopp.ParticleProperties(
                            t1_new, new_property['mass'], new_property['charge']),
                        nb_level + 1
                    )
            return pp

        def _cfg_post_process_remove_neighbour_bonds(cfg):
            """Setup PostProcessRemoveNeighbourBonds"""
            pp = espressopp.integrator.PostProcessRemoveNeighbourBond(self.tm)
            bond_types = [
                x.split('->') for x in cfg['bonds_to_remove'].split(',')
                ]
            # bonds_to_remove=opls_220->opls_220:opls_154:1,opls_268->opls_268:opls_270:1
            for anchor_type, pairs_to_remove in bond_types:
                anchor_type_id = self.topol.used_atomsym_atomtype[anchor_type]
                type_name1, type_name2, nb_level = pairs_to_remove.split(':')
                print('Remove bond anchored to {} at distance {} between {}-{}'.format(
                    anchor_type, nb_level, type_name1, type_name2
                ))
                nb_level = int(nb_level)
                type_pid1 = self.topol.used_atomsym_atomtype[type_name1]
                type_pid2 = self.topol.used_atomsym_atomtype[type_name2]
                pp.add_bond_to_remove(anchor_type_id, nb_level, type_pid1, type_pid2)
            return pp

        def _cfg_post_process_freeze_region(cfg):
            """Setup freeze region."""
            directions = cfg.get('directions', '-x,x,-y,y,-z,z').split(',')
            target_type = cfg['target_type']
            target_type_id = self.topol.atomsym_atomtype[target_type]
            final_type_id = max(self.topol.atomsym_atomtype.values()) + 1
            self.topol.atomsym_atomtype['FREEZE_{}'.format(final_type_id)] = final_type_id
            boxL = self.system.bc.boxL
            if cfg.get('width_type', 'static') == 'ratio':
                width = float(cfg['width'])*boxL
            else:
                width = espressopp.Real3D(float(cfg['width']))

            dir_to_region = {
                '-x': (espressopp.Real3D(0.0), espressopp.Real3D(width[0], boxL[1], boxL[2])),
                '-y': (espressopp.Real3D(0.0), espressopp.Real3D(boxL[0], width[1], boxL[2])),
                '-z': (espressopp.Real3D(0.0), espressopp.Real3D(boxL[0], boxL[1], width[2])),
                'x': (espressopp.Real3D(boxL[0] - width[0], 0, 0), boxL),
                'y': (espressopp.Real3D(0, boxL[1] - width[1], 0), boxL),
                'z': (espressopp.Real3D(0, 0, boxL[2] - width[2]), boxL)}

            for d in directions:
                print('Define region {}: {}-{} with type: {}'.format(
                    d,
                    dir_to_region[d][0],
                    dir_to_region[d][1],
                    target_type_id))
                particle_region = espressopp.ParticleRegion(
                    self.system.storage,
                    self.system.integrator,
                    dir_to_region[d][0],
                    dir_to_region[d][1])
                particle_region.add_type_id(target_type_id)
                change_in_region = espressopp.integrator.ChangeInRegion(
                    self.system, particle_region)
                change_in_region.set_particle_properties(
                    target_type_id, espressopp.ParticleProperties(final_type_id))
                change_in_region.set_flags(target_type_id, True, True)
                self.system.integrator.addExtension(change_in_region)

        def _cfg_post_process_release_molecule(cfg):
            """Setup release molecules."""
            host_type = cfg['host_type']
            target_type = cfg['target_type']
            eq_length = float(cfg['eq_length'])
            alpha = float(cfg['alpha'])
            init_res = float(cfg['init_res'])
            final_type = cfg.get('final_type', target_type)

            # Generate dummy molecules
            max_pid = max(self.topol.atoms)
            dummy_type_id = max(self.topol.atomsym_atomtype.values()) + 1
            self.topol.atomsym_atomtype['DUMMY_{}'.format(dummy_type_id)] = dummy_type_id
            host_pids = sorted([x for x, v in self.topol.atoms.items() if v['type'] == host_type])
            target_type_id = self.topol.atomsym_atomtype[target_type]
            target_properties = self.topol.gt.atomtypes[target_type]
            print('Generate {} of dummy particles (type: {}) linked to {}'.format(
                len(host_pids), dummy_type_id, host_type))
            particle_list = []
            fix_list = []
            for idx, host_pid in enumerate(host_pids):
                host_p = self.system.storage.getParticle(host_pid)
                dummy_pos = host_p.pos + espressopp.Real3D(eq_length, 0.0, 0.0)
                new_pid = max_pid + idx + 1
                fix_list.append((host_pid, new_pid, eq_length))
                particle_list.append((new_pid, dummy_type_id, dummy_pos, target_properties['mass'], new_pid, init_res))
            props = ['id', 'type', 'pos', 'mass', 'res_id', 'lambda_adr']
            self.system.storage.addParticles(particle_list, *props)
            self.system.storage.decompose()

            fix_distance = espressopp.integrator.FixDistances(
                self.system,
                fix_list,
                self.topol.atomsym_atomtype[host_type],
                dummy_type_id)
            self.fix_distance = fix_distance

            fxd_post_process = espressopp.integrator.PostProcessChangeProperty()
            fxd_post_process.add_change_property(
                dummy_type_id,
                espressopp.ParticleProperties(
                    target_type_id,
                    target_properties['mass'],
                    0.0
                ))
            fix_distance.add_postprocess(fxd_post_process)
            self.system.integrator.addExtension(fix_distance)

            basic_dynamic_res = espressopp.integrator.BasicDynamicResolution(self.system, {target_type_id: alpha})
            # If the final_type != target_type then we have to change the type of molecules after resolution reaches
            # 1.0
            final_type_id = target_type_id
            if target_type != final_type:
                final_type_id = self.topol.atomsym_atomtype[final_type]
                final_properties = self.topol.gt.atomtypes[final_type]
                final_particle_properties = espressopp.ParticleProperties(
                            final_type_id,
                            final_properties['mass'],
                            final_properties['charge'],
                            1.0)
                basic_dynamic_res.add_postprocess(
                    espressopp.integrator.PostProcessChangeProperty(
                        target_type_id, final_particle_properties))
                print('Change property of final type {}->{} whenever resolution reaches 1.0'.format(
                    target_type_id, final_type_id))

            self.system.integrator.addExtension(basic_dynamic_res)

            # Because of the dummy particle, we have to use thermal group for thermostat to not
            # thermoset the dummy particle
            self.use_thermal_group = True

            # Observ progress of generating that molecule by checking the total number of target type_id
            if self.cr_observs is None:
                self.cr_observs = {}

            self.cr_observs[(final_type_id, len(particle_list))] = espressopp.analysis.ChemicalConversion(
                self.system, final_type_id, len(particle_list))

            return None

        class_to_cfg = {
            'ChangeNeighboursProperty': _cfg_post_process_change_neighbour,
            'RemoveNeighboursBonds': _cfg_post_process_remove_neighbour_bonds,
            'ReleaseMolecule': _cfg_post_process_release_molecule,
            'FreezeRegion': _cfg_post_process_freeze_region
        }
        for pp_cfg in cfg.values():
            cfg_setup = class_to_cfg[pp_cfg['class']]
            post_process_obj = cfg_setup(pp_cfg['options'])
            if post_process_obj:
                pps.append(post_process_obj)

        return pps

    def setup_reactions(self):
        """Setup reactions.

        Returns:
            The espressopp.integrator.ChemicalReaction extension and the list
            of fixed pair lists with new bonds.

        """
        self.ar_interval = int(self.cfg['general']['interval'])
        ar = espressopp.integrator.ChemicalReaction(
            self.system,
            self.vl,
            self.system.storage,
            self.tm,
            self.ar_interval)
        ar.nearest_mode = self.cfg['general']['nearest']

        fpls = []
        reactions = []

        for group_name, reaction_group in self.cfg['reactions'].items():
            print('Setting reaction group {}'.format(group_name))

            # Setting the interaction for the pairs created by this reaction group.
            fpl = espressopp.FixedPairList(self.system.storage)
            fpls.append(fpl)
            pot_class = eval('espressopp.interaction.{}'.format(reaction_group['potential']))
            # Convert if it's possible, values for float
            pot_options = {}
            for k, v in reaction_group['potential_options'].items():
                try:
                    pot_options[k] = float(v)
                except ValueError:
                    pot_options[k] = v
            print('Setting potential for bond with class {}, options {}'.format(
                reaction_group['potential'], reaction_group['potential_options']))
            potential = pot_class(**pot_options)
            interaction = eval('espressopp.interaction.FixedPairList{}'.format(
                reaction_group['potential']))(self.system, fpl, potential)
            fpl.interaction = interaction
            self.system.addInteraction(interaction, 'fpl_{}'.format(group_name))

            # Setting the post process extensions.
            extensions = self._prepare_group_postprocess(reaction_group['extensions'])

            print('Setting chemical reactions in group')
            for chem_reaction in reaction_group['reaction_list']:
                # Pass connectivity map from group level to reaction level
                chem_reaction['connectivity_map'] = reaction_group['connectivity_map']
                r = self._setup_reaction(chem_reaction, fpl)
                if r is not None:
                    for pp in extensions:
                        r.add_postprocess(pp)
                    ar.add_reaction(r)
                    reactions.append(r)

        return ar, fpls, reactions
