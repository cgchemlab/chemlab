#  Copyright (C) 2016,2017
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

__doc__ = """Main module for setup reactions."""

import collections
import espressopp
import random
from reaction_parser import REACTION_DISSOCATION, REACTION_NORMAL, REACTION_EXCHANGE, EXT_INTEGRATOR, \
    EXT_POSTPROCESS
import reaction_post_process


class SetupReactions:
    """Main class to setup the reactions

    Args:
        system: The espressopp.System object.
        vl: The VerletList that is used to find the reacted pairs.
        topology_manager: The espressopp.integrator.TopologyManager object.
        config: The config file.
    """

    def __init__(self, system, vl, topol, topol_manager, config, args):
        self.missing_fpls = []
        self.system = system
        self.vl = vl
        self.topol = topol
        self.tm = topol_manager
        self.cfg = config
        self.args = args
        self.name2type = topol.atomsym_atomtype
        self.dynamic_types = set()  # Stores the particle types that will change during the reactions.

        self.fix_distances = []
        self.cr_observs = {}  # Observs conversion types.

        # Bond types that will change and has to be observed by dump topology
        self.observed_bondtypes = set()

        # Bond types that has to be placed in separate FixedPairLists
        self.separate_fpls = set()

        # Link type tuple (sorted) with FixedPairList
        self.type2fpl = {}

        self.exclusions_list = []  # For restrict reactions, the exclusion lists has to be extended.

        self.post_process_setup = reaction_post_process.PostProcessSetup(system, topol, topol_manager, args)
        self.post_process_setup.dynamic_types = self.dynamic_types
        self.post_process_setup.observed_bondtypes = self.observed_bondtypes
        self.post_process_setup.cr_observs = self.cr_observs
        self.post_process_setup.fix_distances = self.fix_distances
        self.post_process_setup.separate_fpls = self.separate_fpls
        self.post_process_setup.type2fpl = self.type2fpl

    @property
    def use_thermal_group(self):
        return self.post_process_setup.use_thermal_group

    def _setup_reaction_normal(self, chem_reaction, fpl):
        rl = chem_reaction['reactant_list']

        if chem_reaction.get('connectivity_map'):
            r_class = espressopp.integrator.RestrictReaction
        else:
            r_class = espressopp.integrator.Reaction

        rt1 = rl['type_1']['name']
        rt2 = rl['type_2']['name']
        reaction = r_class(
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
            cutoff=float(chem_reaction['cutoff']))

        self.dynamic_types.add(self.name2type[rl['type_1']['name']])
        self.dynamic_types.add(self.name2type[rl['type_2']['name']])

        print('Setup reaction: {}({})-{}({})'.format(
            rt1, self.name2type[rt1], rt2, self.name2type[rt2]))

        print("Intramolecular bonds: {}".format(chem_reaction['intramolecular']))
        reaction.intramolecular = bool(chem_reaction['intramolecular'])

        reaction.intraresidual = bool(chem_reaction['intraresidual'])
        reaction.is_virtual = bool(chem_reaction['virtual'])

        if 'sigma' in chem_reaction:
            reaction.set_reaction_cutoff(espressopp.integrator.ReactionCutoffRandom(
                chem_reaction['eq_distance'], chem_reaction['sigma'], seed=random.randint(100, 100000)))

        if 'min_cutoff' in chem_reaction:
            reaction.get_reaction_cutoff().min_cutoff = float(chem_reaction['min_cutoff'])
        if 'active' in chem_reaction:
            reaction.active = chem_reaction['active']

        if chem_reaction.get('connectivity_map'):
            print('Reading connectivity map {}, reaction will be restricted'.format(
                chem_reaction['connectivity_map']))
            connectivity_map = open(chem_reaction['connectivity_map'])
            ex_list = set()
            for l in connectivity_map.readlines():
                b1, b2 = map(int, l.strip().split())
                ex_list.add(tuple(sorted(b1, b2)))
            for b1, b2 in ex_list:
                reaction.define_connection(b1, b2)
            self.exclusions_list.extend(list(ex_list))
            print('Restricted to {} connections'.format(len(ex_list)))
            if chem_reaction.get('reaction_type') == REACTION_DISSOCATION:
                reaction.revert = True

        t1_old = self.name2type[rl['type_1']['name']]
        t1_new = self.name2type[rl['type_1']['new_type']]
        t2_old = self.name2type[rl['type_2']['name']]
        t2_new = self.name2type[rl['type_2']['new_type']]
        # Change type if necessary.
        if (rl['type_1']['name'] != rl['type_1']['new_type'] or
                    rl['type_2']['name'] != rl['type_2']['new_type']):
            r_pp = espressopp.integrator.PostProcessChangeProperty()
            if t1_old != t1_new:
                self.dynamic_types.add(t1_old)
                self.dynamic_types.add(t1_new)
                print('Reaction: {}-{}, change type {}->{}'.format(rt1, rt2, t1_old, t1_new))
                new_property = self.topol.gt.atomtypes[rl['type_1']['new_type']]
                r_pp.add_change_property(
                    t1_old,
                    espressopp.integrator.TopologyParticleProperties(
                        type=t1_new, mass=new_property['mass'],
                        q=new_property['charge']))
            if t2_old != t2_new:
                self.dynamic_types.add(t2_old)
                self.dynamic_types.add(t2_new)
                print('Reaction: {}-{}, change type {}->{}'.format(rt1, rt2, t2_old, t2_new))
                new_property = self.topol.gt.atomtypes[rl['type_2']['new_type']]
                r_pp.add_change_property(
                    t2_old,
                    espressopp.integrator.TopologyParticleProperties(
                        type=t2_new, mass=new_property['mass'],
                        q=new_property['charge']))


            reaction.add_postprocess(r_pp)
        return reaction, [(t1_old, t2_old), (t1_new, t2_new)]

    def _setup_reaction_exchange(self, chem_reaction, fpl):
        """
        Setup exchange reaction of the form:
            A[min,max):B[min,max) + C[min,max) -> A(deltaA):C(deltaC) + B(deltaB)
            ^^^^^^^^^^ ^^^^^^^^^     ^^^^^^^^    ^^^^^^     ^^^^^^^     ^^^^^
               type_1   type_2        type_3     type_1     type_3      type_2
        Args:
            chem_reaction: Definition of reaction object from the reaction_parser.
            fpl: The fixed pair list.

        Returns:
            espressopp.integrator.Reaction object
        """

        rl = chem_reaction['reactant_list']
        if chem_reaction.get('connectivity_map'):
            raise RuntimeError('connectivity_map not supported by exchange reaction')
        r_class = espressopp.integrator.Reaction
        rt1 = rl['type_1']
        rt2 = rl['type_2']
        rt3 = rl['type_3']
        reaction = r_class(
            type_1=self.name2type[rt1['name']],
            type_2=self.name2type[rt3['name']],
            delta_1=int(rt1['delta']),
            delta_2=int(rt3['delta']),
            min_state_1=int(rt1['min']),
            max_state_1=int(rt1['max']),
            min_state_2=int(rt3['min']),
            max_state_2=int(rt3['max']),
            rate=float(chem_reaction['rate']),
            fpl=fpl,
            cutoff=float(chem_reaction['cutoff']))
        reaction.is_virtual = True  # We don't mean to make a bond, only to catch the event.

        print('Exchange reaction: {}-{}, restricted when {} in state [{},{})'.format(
            rt1['name'], rt3['name'], rt2['name'], int(rt2['min']), int(rt2['max'])))

        reaction.add_constraint(espressopp.integrator.ReactionConstraintNeighbourState(
            self.name2type[rt2['name']], int(rt2['min']), int(rt2['max'])), 'type_1')

        self.dynamic_types.add(self.name2type[rl['type_1']['name']])
        self.dynamic_types.add(self.name2type[rl['type_2']['name']])
        self.dynamic_types.add(self.name2type[rl['type_3']['name']])
        print('Setup exchange reaction: {}({})-{}({})'.format(
            rt1['name'], self.name2type[rt1['name']], rt3['name'], self.name2type[rt3['name']]))

        reaction.intraresidual = bool(chem_reaction['intraresidual'])
        if 'sigma' in chem_reaction:
            reaction.set_reaction_cutoff(espressopp.integrator.ReactionCutoffRandom(
                chem_reaction['eq_distance'], chem_reaction['sigma'], seed=random.randint(100, 100000)))

        if 'min_cutoff' in chem_reaction:
            reaction.get_reaction_cutoff().min_cutoff = float(chem_reaction['min_cutoff'])

        t1_old = self.name2type[rt1['name']]
        t1_new = self.name2type[rt1['new_type']]
        t2_old = self.name2type[rt2['name']]
        t2_new = self.name2type[rt2['new_type']]
        # Change type if necessary.
        if t1_old != t1_new:
            r_pp = espressopp.integrator.PostProcessChangePropertyByTopologyManager(self.tm)
            self.dynamic_types.add(t1_old)
            self.dynamic_types.add(t1_new)
            print('Exchange reaction: {}-{}, change type {}->{} by TM'.format(rt1['name'], rt3['name'], t1_old, t1_new))
            new_property = self.topol.gt.atomtypes[rt1['new_type']]
            r_pp.add_change_property(
                t1_old,
                espressopp.integrator.TopologyParticleProperties(
                    type=t1_new, mass=new_property['mass'],
                    q=new_property['charge']))
            reaction.add_postprocess(r_pp, 'type_1')

        self.dynamic_types.add(t2_old)
        self.dynamic_types.add(t2_new)
        print('Exchange reaction: {}-{}, change type {}->{} NB'.format(rt1['name'], rt2['name'], t2_old, t2_new))

        new_property = self.topol.gt.atomtypes[rl['type_2']['new_type']]
        print('Change type {}->{} property={}'.format(t2_old, t2_new, new_property))

        # This property change only when E is in certain state min,max
        tm_particle_properties=espressopp.integrator.TopologyParticleProperties(
            type=t2_new, mass=new_property['mass'], q=new_property['charge'], incr_state=int(rt2['delta']))
        tm_particle_properties.set_min_max_state(int(rt2['min']), int(rt2['max']))
        pp = espressopp.integrator.PostProcessChangeNeighboursProperty(self.tm)
        pp.add_change_property(t2_old, tm_particle_properties, 1)
        reaction.add_postprocess(pp, 'type_1')

        return reaction, [(t1_old, t2_old), (t1_new, t2_new)]

    def _setup_reaction_dissocation(self, chem_reaction):
        """Setup dissociation reaction."""
        rl = chem_reaction['reactant_list']

        r_class = espressopp.integrator.DissociationReaction

        rt1 = rl['type_1']['name']
        rt2 = rl['type_2']['name']

        # Get the correct fpl.
        t1, t2 = self.name2type[rl['type_1']['name']], self.name2type[rl['type_2']['name']]
        fpl = self.type2fpl.get((t1, t2))

        reaction = r_class(
            type_1=t1,
            type_2=t2,
            delta_1=int(rl['type_1']['delta']),
            delta_2=int(rl['type_2']['delta']),
            min_state_1=int(rl['type_1']['min']),
            max_state_1=int(rl['type_1']['max']),
            min_state_2=int(rl['type_2']['min']),
            max_state_2=int(rl['type_2']['max']),
            rate=float(chem_reaction['rate']),
            fpl=fpl,
            cutoff=float(chem_reaction['cutoff']))

        if not fpl:
            self.missing_fpls.append((reaction, (t1, t2)))


        if 'diss_rate' in chem_reaction:
            reaction.diss_rate = chem_reaction['diss_rate']

        self.dynamic_types.add(self.name2type[rl['type_1']['name']])
        self.dynamic_types.add(self.name2type[rl['type_2']['name']])

        if 'active' in chem_reaction:
            reaction.active = chem_reaction['active']

        if 'min_cutoff' in chem_reaction:
            reaction.get_reaction_cutoff().min_cutoff = float(chem_reaction['min_cutoff'])

        if 'sigma' in chem_reaction:
            reaction.set_reaction_cutoff(espressopp.integrator.ReactionCutoffRandom(
                chem_reaction['eq_distance'], chem_reaction['sigma'], seed=random.randint(100, 100000)))

        if 'alpha' not in chem_reaction:
            raise RuntimeError('alpha parameter for DissociationReaction not found')
        alpha = float(chem_reaction['alpha'])

        reaction.is_virtual = bool(chem_reaction['virtual'])
        print('Setup dissociation reaction of types: {}({})-{}({})'.format(
            rt1, self.name2type[rt1], rt2, self.name2type[rt2]))

        t1_old = self.name2type[rl['type_1']['name']]
        t1_new = self.name2type[rl['type_1']['new_type']]
        t2_old = self.name2type[rl['type_2']['name']]
        t2_new = self.name2type[rl['type_2']['new_type']]

        self.observed_bondtypes.add((t1_old, t2_old))

        # Change type if necessary.
        r_pp = espressopp.integrator.PostProcessChangeProperty()
        self.dynamic_types.add(t1_old)
        self.dynamic_types.add(t2_old)
        basic_dynamic_res = espressopp.integrator.BasicDynamicResolution(
            self.system, {t1_old: alpha, t2_old: alpha})
        if t1_old != t1_new:
            self.dynamic_types.add(t1_new)
            new_property = self.topol.gt.atomtypes[rl['type_1']['new_type']]
            tpp = espressopp.integrator.TopologyParticleProperties(
                type=t1_new,
                mass=new_property['mass'],
                lambda_adr=1.0,
                q=new_property['charge'])
            basic_dynamic_res.add_postprocess(
                espressopp.integrator.PostProcessChangeProperty(t1_old, tpp))

        r_pp.add_change_property(
            t1_old,
            espressopp.integrator.TopologyParticleProperties(lambda_adr=0.0))

        if t2_old != t2_new:
            self.dynamic_types.add(t2_new)
            new_property = self.topol.gt.atomtypes[rl['type_2']['new_type']]
            tpp = espressopp.integrator.TopologyParticleProperties(
                    type=t2_new,
                    mass=new_property['mass'],
                    lambda_adr=1.0,
                    q=new_property['charge'])
            basic_dynamic_res.add_postprocess(
                espressopp.integrator.PostProcessChangeProperty(t2_old, tpp))
        r_pp.add_change_property(
            t2_old,
            espressopp.integrator.TopologyParticleProperties(lambda_adr=0.0))

        reaction.add_postprocess(r_pp)
        self.system.integrator.addExtension(basic_dynamic_res)

        return reaction, [(t1_old, t2_old), (t1_new, t2_new)]


    def _setup_reaction(self, chem_reaction, fpl):
        """Setup single reaction.

        Args:
            chem_reaction: Dictionary with definition of the reactions
            fpl: The FixedPairList.

        Returns:
            The espressopp.integrator.Reaction object.
        """
        if not chem_reaction['active']:
            return None, None

        # Select reaction class.
        reaction_type2class = {
            REACTION_NORMAL: self._setup_reaction_normal,
            REACTION_EXCHANGE: self._setup_reaction_exchange,
            REACTION_DISSOCATION: self._setup_reaction_dissocation
        }

        if chem_reaction['reaction_type'] not in reaction_type2class:
            raise RuntimeError('reaction_type {} not supported'.format(chem_reaction['reaction_type']))
        return reaction_type2class[chem_reaction['reaction_type']](chem_reaction, fpl)

    def _prepare_group_postprocess(self, cfg):
        """Prepare extensions for the reactions.

        Args:
            cfg: The dictionary with configuration.

        Returns:
            The list of triplets, object, if PostProcess then which particle will be involved and extension type.
            The extension type can be: 'PP' - PostProcess, 'Integrator' - Extension to integrator.
            If either postprocess nor extension to integrator is added this triplet has to be (None, None, None)
        """

        list_of_extensions = collections.defaultdict(list)

        for ext_name, pp_cfg in cfg.items():
            cfg_setup = self.post_process_setup.setup_post_process(pp_cfg['class'])
            post_process_obj = cfg_setup(pp_cfg['options'])
            if post_process_obj:
                if isinstance(post_process_obj, list):
                    list_of_extensions[ext_name].extend(post_process_obj)
                else:
                    list_of_extensions[ext_name].append(post_process_obj)

        return list_of_extensions

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
        if self.cfg['general']['pair_distances_filename']:
            ar.pair_distances_filename = self.cfg['general']['pair_distances_filename']
        if self.cfg['general']['max_per_interval'] > 0:
            ar.max_per_interval = self.cfg['general']['max_per_interval']

        fpls = []
        reactions = []
        extensions_to_integrator = []

        fpl_def = collections.namedtuple('FPLDef', ['fpl', 'type_list'])

        self.reaction_index = {}
        reaction_idx = 0

        for group_name, reaction_group in self.cfg['reactions'].items():
            print('Setting reaction group {}'.format(group_name))

            # Setting the interaction for the pairs created by this reaction group.
            if self.args.t_hybrid_bond > 0:
                fpl = espressopp.FixedPairListLambda(self.system.storage, 0.0)
                interaction_class = eval('espressopp.interaction.FixedPairListLambda{}'.format(
                    reaction_group['potential']))
            else:
                fpl = espressopp.FixedPairList(self.system.storage)
                interaction_class = eval('espressopp.interaction.FixedPairList{}'.format(
                    reaction_group['potential']))

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

            interaction = interaction_class(self.system, fpl, potential)
            fpl.interaction = interaction
            self.system.addInteraction(interaction, 'chem_fpl_{}'.format(group_name))

            # Setting the post process extensions.
            group_extensions = self._prepare_group_postprocess(reaction_group['extensions'])
            extensions_to_integrator = []
            extensions_to_reactions = collections.defaultdict(list)
            for k, v in group_extensions.items():
                for x in v:
                    if x.ext is not None:
                        if x.ext_type == EXT_INTEGRATOR:
                            extensions_to_integrator.append(x.ext)
                        elif x.ext_type == EXT_POSTPROCESS:
                            extensions_to_reactions[k].append(x)
                        else:
                            raise RuntimeError('Wrong ext_type={}'.format(x.ext_type))

            # Process the reactions.
            reaction_type_list = []
            print('Setting chemical reactions in group')
            for chem_reaction in reaction_group['reaction_list']:
                # Pass connectivity map from group level to reaction level
                chem_reaction['connectivity_map'] = reaction_group['connectivity_map']
                # Dissociation reaction will be done in second run.
                if chem_reaction['reaction_type'] == REACTION_DISSOCATION:
                    continue
                r, reaction_types = self._setup_reaction(chem_reaction, fpl)
                if r is not None:
                    reaction_type_list.extend(reaction_types)
                    for ext_name, extensions in extensions_to_reactions.items():
                        if ext_name in chem_reaction['exclude_extensions']:
                            print('Skip extension: {} ({})'.format(ext_name, chem_reaction['equation']))
                        else:
                            for extension in extensions:
                                print('Add extension {} to {} ({}): {}'.format(
                                    ext_name, chem_reaction['equation'], extension.pp_type, extension.ext))
                                if extension.pp_type:
                                    r.add_postprocess(extension.ext, extension.pp_type)
                                else:
                                    r.add_postprocess(extension.ext)
                    ar.add_reaction(r)
                    self.reaction_index[reaction_idx] = chem_reaction['equation']
                    reactions.append(r)
                    reaction_idx += 1
                    for t1, t2 in reaction_types:
                        self.type2fpl[(t1, t2)] = fpl
                        self.type2fpl[(t2, t1)] = fpl

            # Now process dissociation reactions.
            # Pitfal, It can only sees fixed pair lists in given group
            for chem_reaction in reaction_group['reaction_list']:
                if chem_reaction['reaction_type'] != REACTION_DISSOCATION:  # only dissociation
                    continue
                r, reaction_types = self._setup_reaction_dissocation(chem_reaction)
                if r is not None:
                    #reaction_type_list.extend(reaction_types)
                    self.separate_fpls.add(tuple(reaction_types[0]))
                    for ext_name, extensions in extensions_to_reactions.items():
                        if ext_name in chem_reaction['exclude_extensions']:
                            print('Skip extension: {} ({})'.format(ext_name, chem_reaction['equation']))
                        else:
                            for extension in extensions:
                                print('Add extension {} to {} ({}): {}'.format(
                                    ext_name, chem_reaction['equation'], extension.pp_type, extension.ext))
                                if extension.pp_type:
                                    r.add_postprocess(extension.ext, extension.pp_type)
                                else:
                                    r.add_postprocess(extension.ext)
                    ar.add_reaction(r)
                    self.reaction_index[reaction_idx] = chem_reaction['equation']
                    reactions.append(r)
                    reaction_idx += 1

            fpls.append(fpl_def(fpl, set(reaction_type_list)))

        return ar, fpls, reactions, extensions_to_integrator


    def rebuild_fixed_pair_lists(self):
        """If the tuple for dissociation reaction is not found then the TopologyManager has to be used"""
        for r, (t1, t2) in self.missing_fpls:
            fpl = self.tm.get_fixed_pair_list(t1, t2)
            if fpl:
                r.set_fixed_pair_list(fpl)
                print fpl.getAllBonds()
            else:
                raise RuntimeError('Fixed pair list for {}-{} not found'.format(t1, t2))
