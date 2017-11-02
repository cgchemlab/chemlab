#  Copyright (C) 2016,2017
#      Jakub Krajniak (jkrajniak at gmail.com)
#      Zidan Zhang (zidan.zhang at kuleuven.be)
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

import cPickle
import espressopp
import os
import re

from reaction_parser import EXT_INTEGRATOR, EXT_POSTPROCESS

PP_TYPE_1 = 'type_1'
PP_TYPE_2 = 'type_2'
PP_BOTH = 'both'

# Belows are the sub-functions to create PostProcess objects.
output_triplet = collections.namedtuple('Extension', ['ext', 'pp_type', 'ext_type'])


class PostProcessSetup(object):
    def __init__(self, system, topol, topol_manager, args):
        self.fix_distances = []
        self.system = system
        self.tm = topol_manager
        self.name2type = topol.atomsym_atomtype
        self.topol = topol
        self.dynamic_types = None
        self.observed_bondtypes = None
        self.cr_observs = {}
        self.use_thermal_group = False
        self.simulation_args = args

    def setup_post_process(self, post_process_type):
        """Process post-processing options"""

        # The entry points
        pp_type_to_cfg = {
            'ChangeNeighboursProperty': self._setup_post_process_change_neighbour,
            'RemoveNeighboursBonds': self._setup_post_process_remove_neighbour_bonds,
            'ReleaseMolecule': self._setup_post_process_release_molecule,
            'JoinMolecule': self._setup_post_process_join_molecule,
            'FreezeRegion': self._setup_post_process_freeze_region,
            'ChangeParticleType': self._setup_change_particle_type,
            'ATRPActivator': self._setup_atrp_activator
        }
        if post_process_type not in pp_type_to_cfg:
            raise RuntimeError('Wrong post_process_type {}'.format(post_process_type))
        return pp_type_to_cfg[post_process_type]

    def _setup_post_process_change_neighbour(self, cfg):
        """Setup PostProcessChangeNeighbourProperties"""
        pp = espressopp.integrator.PostProcessChangeNeighboursProperty(self.tm)
        type_transfers = [
            x.split('->') for x in cfg['type_transfers'].split(',')]
        re_opt = re.compile(r'(?P<type_name>\w+)\(?(?P<options>[a-zA-Z0-9_=,]*)\)?')
        invoke_on = cfg.get('invoke_on')
        print('Setup PostProcessChangeNeighbourProperties, {}'.format(cfg['type_transfers']))
        for old_type, new_type in type_transfers:
            old_type, nb_level = old_type.split(':')
            nb_level = int(nb_level)
            opt_match = re_opt.match(new_type)
            t1_old = self.name2type[old_type]
            if opt_match:
                new_type, options = opt_match.groups()
                t1_new = self.name2type[new_type]
                new_property_def = self.topol.gt.atomtypes[new_type]
                if 'state' not in new_property_def:
                    raise RuntimeError(
                        ('Please define initial atom state in [ atomstate ]'
                         'section of your topology for atom type {}'.format(new_type)))
                new_properties_args = {
                    'type': t1_new,
                    'mass': new_property_def['mass'],
                    'q': new_property_def['charge'],
                    'state': new_property_def['state']
                }
                if options:
                    additional_properties = {}
                    exec(options, {}, additional_properties)
                    new_properties_args.update(additional_properties)
                new_property = espressopp.integrator.TopologyParticleProperties(**new_properties_args)

                self.dynamic_types.add(t1_old)
                self.dynamic_types.add(t1_new)
                print('PostProcessChangeNeighbourProperties: {}->{}: {} at {}'.format(t1_old, t1_new, pp, nb_level))
                pp.add_change_property(t1_old, new_property, nb_level)
                #pp.add_change_property(t1_old, new_property, nb_level+1)

        return output_triplet(pp, invoke_on, EXT_POSTPROCESS)

    def _setup_post_process_remove_neighbour_bonds(self, cfg):
        """Setup PostProcessRemoveNeighbourBonds"""
        pp = espressopp.integrator.PostProcessRemoveNeighbourBond(self.tm)
        bond_types = [
            x.split('->') for x in cfg['bonds_to_remove'].split(',')
        ]
        invoke_on = cfg.get('invoke_on', 'both')
        # bonds_to_remove=opls_220->opls_220:opls_154:1,opls_268->opls_268:opls_270:1
        for anchor_type, pairs_to_remove in bond_types:
            anchor_type_id = self.topol.used_atomsym_atomtype[anchor_type]
            type_name1, type_name2, nb_level = pairs_to_remove.split(':')
            print('Remove bond anchored to {} at distance {} between {}-{} invoke_on={}'.format(
                anchor_type, nb_level, type_name1, type_name2, invoke_on
            ))
            nb_level = int(nb_level)
            type_pid1 = self.topol.used_atomsym_atomtype[type_name1]
            type_pid2 = self.topol.used_atomsym_atomtype[type_name2]
            print((anchor_type_id, nb_level, type_pid1, type_pid2))
            pp.add_bond_to_remove(anchor_type_id, nb_level, type_pid1, type_pid2)
            self.observed_bondtypes.add(tuple(sorted([type_pid1, type_pid2])))
        return output_triplet(pp, invoke_on, EXT_POSTPROCESS)

    def _setup_post_process_freeze_region(self, cfg):
        """Setup freeze region."""
        directions = cfg.get('directions', '-x,x,-y,y,-z,z').split(',')
        target_type = cfg['target_type']
        target_type_id = self.topol.atomsym_atomtype[target_type]
        if cfg.get('stats_file'):
            stats_file = cfg.get('stats_file')
        else:
            stats_file = '{}_{}_freeze_stats.dat'.format(
                self.simulation_args.output_prefix, self.simulation_args.rng_seed)
        final_type_id = max(self.topol.atomsym_atomtype.values()) + 1
        print('Freeze region with particles of type {}, change type to {}'.format(target_type_id, final_type_id))
        self.topol.atomsym_atomtype['FREEZE_{}'.format(final_type_id)] = final_type_id
        boxL = self.system.bc.boxL
        if cfg.get('width_type', 'static') == 'ratio':
            width = float(cfg['width']) * boxL
        else:
            width = espressopp.Real3D(float(cfg['width']))

        remove_particles = eval(cfg.get('remove_particles', 'False'))
        prob = float(cfg.get('prob')) if cfg.get('prob') else None
        p_num = int(cfg.get('p_num')) if cfg.get('p_num') else None
        p_percentage = float(cfg.get('p_percentage')) if cfg.get('p_percentage') else None
        if p_percentage and (p_percentage > 1.0 or p_percentage < 0.0):
            raise RuntimeError('p_percentage not in the range (0.0, 1.0)')

        if prob:
            print('Freeze in region with prob: {}'.format(prob))
        elif p_num:
            print('Freeze in region with p_num: {}'.format(p_num))
        elif p_percentage:
            print('Freeze in region with percentage: {}'.format(p_percentage))

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
                self.system, particle_region, prob=prob, p_num=p_num, p_num_percentage=p_percentage)
            change_in_region.set_particle_properties(
                target_type_id, espressopp.integrator.TopologyParticleProperties(type=final_type_id))
            change_in_region.set_flags(target_type_id, reset_velocity=True, reset_force=True,
                                       remove_particle=remove_particles)
            change_in_region.stats_filename = stats_file

        return output_triplet(change_in_region, None, EXT_INTEGRATOR)

    def _setup_post_process_release_molecule(self, cfg):
        """Setup release molecules."""
        host_type = cfg['host_type']
        target_type = cfg['target_type']
        eq_length = float(cfg['eq_length'])
        alpha = float(cfg['alpha'])
        init_res = float(cfg['init_res'])
        final_type = cfg.get('final_type', target_type)
        cache_file = cfg.get('cache_file')

        replicate = int(cfg.get('replicate', 1))
        release_on = cfg.get('release_on', 'type')  # bond or type
        if release_on not in ['bond', 'type']:
            raise RuntimeError('Wrong keyword release_on {}, only: bond or type'.format(release_on))
        release_count = int(cfg.get('release_count', 1))
        release_host = cfg.get('release_host', 'both')
        if release_host not in ['type_1', 'type_2', 'both']:
            raise RuntimeError('Wrong keyword release_host {}, only left, right, both'.format(release_host))

        # Generate dummy molecules
        max_pid = max(self.topol.atoms)
        dummy_type_id = max(self.topol.atomsym_atomtype.values()) + 1
        self.topol.add_new_atomtype(dummy_type_id, 'DUMMY_{}'.format(dummy_type_id), False)
        host_pids = sorted([x for x, v in self.topol.atoms.items() if v['type'] == host_type])
        target_type_id = self.topol.atomsym_atomtype[target_type]
        target_properties = self.topol.gt.atomtypes[target_type]
        print('Generate {} of dummy particles (type: {}) linked to {}'.format(
            len(host_pids) * replicate, dummy_type_id, host_type))

        # Creates list of dummy particles.
        particle_list = []
        fix_list = []
        props = ['id', 'type', 'pos', 'mass', 'res_id', 'lambda_adr', 'state']
        if cache_file is None or not os.path.exists(cache_file):
            dummy_idx = max_pid + 1
            for idx, host_pid in enumerate(host_pids):
                host_p = self.system.storage.getParticle(host_pid)
                for _ in range(replicate):
                    dummy_pos = host_p.pos + espressopp.Real3D(eq_length, 0.0, 0.0)
                    fix_list.append((host_pid, dummy_idx, eq_length))
                    particle_list.append((
                        dummy_idx,
                        dummy_type_id,
                        dummy_pos,
                        target_properties['mass'],
                        dummy_idx,
                        init_res,
                        target_properties.get('state', 0)))
                    dummy_idx += 1
        if cache_file:
            if os.path.exists(cache_file):
                print('Read dummy particles from {}'.format(cache_file))
                with open(cache_file, 'rb') as in_file:
                    particle_list, fix_list, props = cPickle.load(in_file)
            else:
                print('Save data to {}'.format(cache_file))
                with open(cache_file, 'wb') as in_file:
                    cPickle.dump((particle_list, fix_list, props), in_file)

        self.system.storage.addParticles(particle_list, *props)
        self.system.storage.decompose()

        reaction_post_process = None

        if release_on == 'type':
            print('Release particle on type change {}-{}'.format(host_type, dummy_type_id))
            fix_distance = espressopp.integrator.FixDistances(
                self.system,
                fix_list,
                self.topol.atomsym_atomtype[host_type],
                dummy_type_id)
        else:  # do not remove fix when change of type
            print('Release particle on new bond created, count: {}'.format(release_count))
            fix_distance = espressopp.integrator.FixDistances(self.system, fix_list)
            # Remove by post process in the reaction
            reaction_post_process = espressopp.integrator.PostProcessReleaseParticles(fix_distance, release_count)
        self.fix_distances.append(fix_distance)

        fxd_post_process = espressopp.integrator.PostProcessChangeProperty()
        fxd_post_process.add_change_property(
            dummy_type_id,
            espressopp.integrator.TopologyParticleProperties(
                type=target_type_id,
                mass=target_properties['mass'],
                lambda_adr=0.0))
        fix_distance.add_postprocess(fxd_post_process)
        self.system.integrator.addExtension(fix_distance)

        basic_dynamic_res = espressopp.integrator.BasicDynamicResolution(self.system, {target_type_id: alpha})
        # If the final_type != target_type then we have to change the type of molecules after resolution reaches
        # 1.0
        final_type_id = target_type_id
        if target_type != final_type:
            final_type_id = self.topol.atomsym_atomtype[final_type]
            final_properties = self.topol.gt.atomtypes[final_type]
            final_particle_properties = espressopp.integrator.TopologyParticleProperties(
                type=final_type_id,
                mass=final_properties['mass'],
                q=final_properties['charge'],
                state=final_properties.get('state', 0),
                lambda_adr=1.0)
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

        return output_triplet(reaction_post_process, release_host, EXT_POSTPROCESS)

    def _setup_post_process_join_molecule(self, cfg):
        host_type = cfg['host_type']
        host_type_id = self.topol.atomsym_atomtype[host_type]
        target_type = cfg['target_type']
        target_type_id = self.topol.atomsym_atomtype[target_type]
        target_properties = self.topol.gt.atomtypes[target_type]
        eq_length = float(cfg['eq_length'])
        init_res = float(cfg['init_res'])
        final_type = cfg.get('final_type', target_type)
        final_type_id = self.topol.atomsym_atomtype[final_type]

        dummy_type_id = max(self.topol.atomsym_atomtype.values()) + 1
        self.topol.add_new_atomtype(dummy_type_id, 'DUMMY_{}'.format(dummy_type_id), False)
        print('PostProcessJoinMolecule: create dummy type DUMMY_{}'.format(dummy_type_id))

        fd = espressopp.integrator.FixDistances(self.system, [], host_type_id, dummy_type_id)
        fxd_post_process = espressopp.integrator.PostProcessChangeProperty()
        fxd_post_process.add_change_property(
            dummy_type_id,
            espressopp.integrator.TopologyParticleProperties(
                type=target_type_id,
                mass=target_properties['mass'],
                state=target_properties.get('state', 0),
                lambda_adr=init_res))
        fd.add_postprocess(fxd_post_process)
        self.system.integrator.addExtension(fd)
        self.fix_distances.append(fd)

        pp_join_particles = espressopp.integrator.PostProcessJoinParticles(fd, eq_length)
        dummy_pp = espressopp.integrator.PostProcessChangeProperty()
        dummy_pp.add_change_property(
            final_type_id,
            espressopp.integrator.TopologyParticleProperties(
                type=dummy_type_id, lambda_adr=init_res, state=target_properties.get('state', 0)))

        self.use_thermal_group = True

        return [
            output_triplet(pp_join_particles, 'type_1', EXT_POSTPROCESS),
            output_triplet(dummy_pp, 'type_2', EXT_POSTPROCESS)]

    def _setup_change_particle_type(self, cfg):
        interval = int(cfg['interval'])
        old_type_id = int(cfg['type_id'])
        new_type_id = int(cfg['new_type_id'])
        num_particles = int(cfg['num_particles'])

        change_type = espressopp.integrator.ChangeParticleType(
            self.system,
            interval,
            num_particles,
            old_type_id,
            new_type_id)

        return output_triplet(change_type, None, EXT_INTEGRATOR)

    def _setup_atrp_activator(self, cfg):
        interval = int(cfg['interval'])
        num_particles = int(cfg['num_particles'])
        select_from_all = int(cfg.get('select_from_all', 1))
        ratio_activator = float(cfg['ratio_activator'])
        ratio_deactivator = float(cfg['ratio_deactivator'])
        delta_catalyst = float(cfg['delta_catalyst'])
        k_activate = float(cfg['k_activate'])
        k_deactivate = float(cfg['k_deactivate'])
        stats_file = cfg.get('stats_file', '{}_{}_atrp_stats.dat'.format(
            self.simulation_args.output_prefix, self.simulation_args.rng_seed))

        atrp_activator = espressopp.integrator.ATRPActivator(
            self.system, interval, num_particles, ratio_activator, ratio_deactivator,
            delta_catalyst, k_activate, k_deactivate)
        atrp_activator.stats_filename = stats_file
        atrp_activator.select_from_all = select_from_all
        options = [x.split('->') for x in cfg['options'].split(';')]
        print('Settings ATRP activator extension')
        print('ATRPActivator.interval={} num_part={} select_from_all={}'.format(
            interval, num_particles, select_from_all))
        re_reactant = re.compile(r'(?P<name>\w+)\((?P<state>\d+),\s*(?P<flag>[AD]{1,2})\)')
        re_product = re.compile(r'(?P<new_type>\w+)\((?P<delta>[0-9-]+)\)')
        for to_process, after_process in options:
            reactant = re_reactant.match(to_process).groupdict()
            product = re_product.match(after_process).groupdict()
            reactant_type_id = self.topol.atomsym_atomtype[reactant['name']]
            reactant_state = int(reactant['state'])
            reactant_flag = reactant['flag'] == 'DA'
            delta_state = int(product['delta'])
            product_type_id = self.topol.atomsym_atomtype[product['new_type']]
            product_property = self.topol.gt.atomtypes[product['new_type']]
            if reactant['flag'] not in ['A', 'DA']:
                raise RuntimeError('Flag {} not "A" or "DA"'.format(reactant['flag']))
            atrp_activator.add_reactive_center(
                type_id=reactant_type_id,
                state=reactant_state,
                is_activator=reactant_flag,
                new_property=espressopp.integrator.TopologyParticleProperties(type=product_type_id,
                                                                              mass=product_property['mass'],
                                                                              q=product_property['charge']),
                delta_state=delta_state)
            print('ATRPActivator: added {}->{} state={} is_activator={} delta_state={}'.format(
                to_process, after_process, reactant_state, reactant_flag, delta_state))

        return output_triplet(atrp_activator, None, EXT_INTEGRATOR)
