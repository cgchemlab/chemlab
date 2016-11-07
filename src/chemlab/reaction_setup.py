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
        self.system = system
        self.vl = vl
        self.topol = topol
        self.tm = topol_manager
        self.cfg = config
        self.args = args
        self.name2type = topol.atomsym_atomtype
        self.dynamic_types = set()  # Stores the particle types that will change during the reactions.

        self.use_thermal_group = False
        self.fix_distance = None
        self.cr_observs = None  # Observs conversion types.

        # Bond types that will change and has to be observed by dump topology
        self.observed_bondtypes = set()

        self.exclusions_list = []  # For restrict reactions, the exclusion lists has to be extended.

        self.post_process_setup = reaction_post_process.PostProcessSetup(system, topol, topol_manager)
        self.post_process_setup.dynamic_types = self.dynamic_types
        self.post_process_setup.observed_bondtypes = self.observed_bondtypes
        self.post_process_setup.cr_observs = self.cr_observs

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
            cutoff=float(chem_reaction.get('cutoff', 0.0)))

        self.dynamic_types.add(self.name2type[rl['type_1']['name']])
        self.dynamic_types.add(self.name2type[rl['type_2']['name']])

        print('Setup reaction: {}({})-{}({})'.format(
            rt1, self.name2type[rt1], rt2, self.name2type[rt2]))
        if 'intramolecular' in chem_reaction:
            print('Warning, tag intramolecular not used anymore!')

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

            reaction.add_postprocess(r_pp)

        return reaction

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
            cutoff=float(chem_reaction.get('cutoff', 0.0)))
        reaction.is_virtual = True  # We don't mean to make a bond, only to catch the event.

        self.dynamic_types.add(self.name2type[rl['type_1']['name']])
        self.dynamic_types.add(self.name2type[rl['type_2']['name']])
        self.dynamic_types.add(self.name2type[rl['type_3']['name']])
        print('Setup reaction: {}({})-{}({})'.format(
            rt1['name'], self.name2type[rt1['name']], rt3['name'], self.name2type[rt3['name']]))

        reaction.intraresidual = bool(chem_reaction['intraresidual'])
        if 'sigma' in chem_reaction:
            reaction.set_reaction_cutoff(espressopp.integrator.ReactionCutoffRandom(
                chem_reaction['eq_distance'], chem_reaction['sigma'], seed=random.randint(100, 100000)))

        if 'min_cutoff' in chem_reaction:
            reaction.get_reaction_cutoff().min_cutoff = float(chem_reaction['min_cutoff'])

        # Change type if necessary.
        if rt1['name'] != rt1['new_type'] or rt2['name'] != rt2['new_type']:
            r_pp = espressopp.integrator.PostProcessChangeProperty()
            t1_old = self.name2type[rt1['name']]
            t1_new = self.name2type[rt1['new_type']]
            if t1_old != t1_new:
                self.dynamic_types.add(t1_old)
                self.dynamic_types.add(t1_new)
                print('Reaction: {}-{}, change type {}->{}'.format(rt1['name'], rt2['name'], t1_old, t1_new))
                new_property = self.topol.gt.atomtypes[rt1['new_type']]
                r_pp.add_change_property(
                    t1_old,
                    espressopp.ParticleProperties(
                        t1_new, new_property['mass'],
                        new_property['charge']))

            t2_old = self.name2type[rt2['name']]
            t2_new = self.name2type[rt2['new_type']]
            if t2_old != t2_new:
                self.dynamic_types.add(t2_old)
                self.dynamic_types.add(t2_new)
                print('Reaction: {}-{}, change type {}->{}'.format(rt1['name'], rt2['name'], t2_old, t2_new))
                new_property = self.topol.gt.atomtypes[rl['type_2']['new_type']]
                r_pp.add_change_property(
                    t2_old,
                    espressopp.ParticleProperties(
                        t2_new, new_property['mass'],
                        new_property['charge']))

            reaction.add_postprocess(r_pp)
        return reaction


    def _setup_reaction(self, chem_reaction, fpl):
        """Setup single reaction.

        Args:
            chem_reaction: Dictionary with definition of the reactions
            fpl: The FixedPairList.

        Returns:
            The espressopp.integrator.Reaction object.
        """
        if not chem_reaction['active']:
            return None
        # Select reaction class.
        reaction_type2class = {
            REACTION_NORMAL: self._setup_reaction_normal,
            REACTION_EXCHANGE: self._setup_reaction_exchange,
        }

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

        list_of_extensions = []

        for pp_cfg in cfg.values():
            cfg_setup = self.post_process_setup.setup_post_process(pp_cfg['class'])
            post_process_obj = cfg_setup(pp_cfg['options'])
            if post_process_obj:
                if isinstance(post_process_obj, list):
                    list_of_extensions.extend(post_process_obj)
                else:
                    list_of_extensions.append(post_process_obj)

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

        fpls = []
        reactions = []
        extensions_to_integrator = []

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

            interaction = interaction_class(self.system, fpl, potential)
            fpl.interaction = interaction
            self.system.addInteraction(interaction, 'fpl_{}'.format(group_name))

            # Setting the post process extensions.
            extensions = self._prepare_group_postprocess(reaction_group['extensions'])
            extensions_to_integrator = [x.ext for x in extensions
                                        if x.ext_type == EXT_INTEGRATOR and x.ext is not None]
            extensions_to_reactions = [x for x in extensions
                                       if x.ext_type == EXT_POSTPROCESS and x.ext is not None]

            print('Setting chemical reactions in group')
            for chem_reaction in reaction_group['reaction_list']:
                # Pass connectivity map from group level to reaction level
                chem_reaction['connectivity_map'] = reaction_group['connectivity_map']
                r = self._setup_reaction(chem_reaction, fpl)
                if r is not None:
                    for extension in extensions_to_reactions:
                        if extension.pp_type:
                            r.add_postprocess(extension.ext, extension.pp_type)
                        else:
                            r.add_postprocess(extension.ext)
                    ar.add_reaction(r)
                    reactions.append(r)

        return ar, fpls, reactions, extensions_to_integrator
