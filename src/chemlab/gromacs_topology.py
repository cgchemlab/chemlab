#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import collections
import os
import espressopp
import math

class FileBuffer():
    def __init__(self):
        self.linecount = 0
        self.lines = []
        self.pos = 0

    def appendline(self, line):
        self.lines.append(line)

    def readline(self):
        try:
            line = self.lines[self.pos]
        except:
            return ''
        self.pos += 1
        return line

    def readlastline(self):
        try:
            line = self.lines[self.pos-1]
        except:
            return ''
        return line

    def seek(self, p):
        self.pos = p

    def tell(self):
        return self.pos


def FillFileBuffer(fname, filebuffer, cwd=None, defines=None):
    if cwd is None:
        cwd = '.'
    if defines is None:
        defines = {}
    f = open(os.path.join(cwd, fname), 'r')
    for line in f:
        if line.startswith(';'):
            continue
        if "include" in line:
            name = line.split()[1].strip('\"')
            cwd_name = os.path.dirname(name)
            if cwd_name != '':
                cwd = cwd_name
            FillFileBuffer(name, filebuffer, cwd, defines)
        elif 'define' in line:
            t = line.strip().split()
            if len(t) > 2:
                defines[t[1]] = ' '.join(t[2:])
        else:
            l = line.rstrip('\n')
            if l:
                filebuffer.appendline(l)

    f.close()
    return


def PostProcessFileBuffer(filebuffer, defines):
    """Replace all defines with the value from the dictionary."""
    ret_fb = FileBuffer()
    define_keys = set(defines)
    for line in filebuffer.lines:
        line = line.strip()
        if line:
            if not (line.startswith(';') or line.startswith('#define')
                    or line.startswith('#include') or line.startswith('#ifdef')
                    or line.startswith('#ifndef')):
                def_key = set.intersection(set(map(str.strip, line.split())), define_keys)
                if def_key:
                    def_key = def_key.pop()
                    ret_fb.appendline(
                        line.replace(def_key, defines[def_key]))
                else:
                    ret_fb.appendline(line)
            else:
                ret_fb.appendline(line)
    return ret_fb


def convertc6c12(c6, c12, cr):
    if cr == 1:
        if c12 == 0.0:
            return 1.0, 0.0
        sig = pow(c12/c6, 1.0/6.)
        if sig > 0.0:
            eps = 0.25*c6*pow(sig, -6.0)
        else:
            eps = 0.0
        return sig, eps
    else:
        return c6, c12


class GromacsTopology:
    def __init__(self, input_topol):
        self.input_file = input_topol
        self.content = None
        self.data = {}

        self.atomsym_atomtype = {}
        self.atomtype_atomsym = {}

        self.bonds = {}
        self.angles = {}
        self.dihedrals = {}

        self.bondparams = {}
        self.angleparams = {}
        self.dihedralparams = {}

    def read(self):
        fb = FileBuffer()
        defines = {}
        FillFileBuffer(self.input_file, fb, defines=defines)
        f = PostProcessFileBuffer(fb, defines)
        self.gt = espressopp.tools.chemlab.files_io.GROMACSTopologyFile(self.input_file)
        self.topol = self.gt
        self.gt.content = f.lines
        self.gt.read()
        if len(self.gt.molecules) > 1:
            raise RuntimeError('Multiple molecules are not supported')

        self._prepare_data()

    def _prepare_data(self):
        # Generate atom types from atom symbols.
        self.atomsym_atomtype = {
            k: idx for idx, k in enumerate(sorted(self.gt.atomtypes))
        }
        self.atomtype_atomsym = {v: k for k, v in self.atomsym_atomtype.items()}
        self.atomparams = {}
        self.atom_id_params = {}
        self.atom_type_params = {}
        combinationrule = self.topol.defaults['combinationrule']
        for at_id in sorted(self.gt.atoms):
            at_data = self.gt.atoms[at_id]
            at_type = self.gt.atomtypes[at_data.atom_type]
            at_key = '{}-{}'.format(at_data.chain_name, at_data.name)
            self.atomparams[at_key] = {
                'molecule': at_data.chain_name,
                'type': at_data.atom_type,
                'sig': at_type['sigma'],
                'eps': at_type['epsilon'],
                'type_id': self.atomsym_atomtype[at_data.atom_type],
                'state': at_type.get('state', 0)
            }
            if at_data.charge:
                self.atomparams[at_key]['charge'] = at_data.charge
            else:
                self.atomparams[at_key]['charge'] = at_type['charge']
            if at_data.mass:
                self.atomparams[at_key]['mass'] = at_data.mass
            else:
                self.atomparams[at_key]['mass'] = at_type['mass']
            sig, eps = convertc6c12(
                at_type['sigma'], at_type['epsilon'], combinationrule)
            self.atomparams[at_key]['sig'] = sig
            self.atomparams[at_key]['eps'] = eps
            self.atom_id_params[at_id] = self.atomparams[at_key]
            self.atom_type_params[self.atomparams[at_key]['type_id']] = self.atomparams[at_key]

        # Update non_bonded params
        for k, v in self.topol.nonbond_params.items():
            if v['func'] == 1 and self.topol.defaults['combinationrule'] == 1:
                c6 = v['params'][0]
                c12 = v['params'][1]
                sig, eps = convertc6c12(c6, c12, combinationrule)
                v['params'][0] = sig
                v['params'][1] = eps

        self._prepare_bondedparams()
        self._prepare_bondedlists()
        self._prepare_exclusionlists()

    def _prepare_bondedlists(self):
        """Replicate bonded lists."""
        n_atoms = len(self.gt.atoms)
        n_mols = self.gt.molecules.values()[0]
        if len(self.gt.molecules.values()) > 1:
            raise RuntimeError('Single molecule is supported, found: {}'.format(self.gt.molecules))

        self.atoms = {
            k+(mol*n_atoms): v for mol in range(n_mols)
            for k, v in self.atom_id_params.items()
        }

        self.bonds = self._replicate_lists(
            n_mols, n_atoms, self.gt.bonds)
        self.angles = self._replicate_lists(
            n_mols, n_atoms, self.gt.angles)
        self.dihedrals = self._replicate_lists(
            n_mols, n_atoms, self.gt.dihedrals)
        self.pairs = self._replicate_lists(
            n_mols, n_atoms, self.gt.pairs)

    def _prepare_exclusionlists(self):
        self.exclusions = {tuple(sorted(x)) for x in self.bonds.keys()[:]}
        self.exclusions.update(
            {tuple(sorted([x[0], x[2]])) for x in self.angles})
        self.exclusions.update(
            {tuple(sorted([x[0], x[3]])) for x in self.dihedrals})

    def _prepare_bondedparams(self):
        """Prepares bonded params to use with FixedListTypes interaction."""
        bcount = 0
        for i in self.gt.bondtypes:
            for j in self.gt.bondtypes[i]:
                bcount += 1
                params = self.gt.bondtypes[i][j]
                t1 = self.atomsym_atomtype[i]
                t2 = self.atomsym_atomtype[j]
                self.bondparams[(t1, t2)] = params
        assert bcount == len(self.bondparams)
        acount = 0
        for i in self.gt.angletypes:
            for j in self.gt.angletypes[i]:
                for k in self.gt.angletypes[i][j]:
                    acount += 1
                    params = self.gt.angletypes[i][j][k]
                    t1 = self.atomsym_atomtype[i]
                    t2 = self.atomsym_atomtype[j]
                    t3 = self.atomsym_atomtype[k]
                    self.angleparams[(t1, t2, t3)] = params
        assert acount == len(self.angleparams)
        dcount = 0
        for i in self.gt.dihedraltypes:
            for j in self.gt.dihedraltypes[i]:
                for k in self.gt.dihedraltypes[i][j]:
                    for l in self.gt.dihedraltypes[i][j][k]:
                        dcount += 1
                        params = self.gt.dihedraltypes[i][j][k][l]
                        t1 = self.atomsym_atomtype[i]
                        t2 = self.atomsym_atomtype[j]
                        t3 = self.atomsym_atomtype[k]
                        t4 = self.atomsym_atomtype[l]
                        self.dihedralparams[(t1, t2, t3, t4)] = params
        assert dcount == len(self.dihedralparams)

    def _replicate_lists(self, n_mols, n_atoms, input_list, shift=0):
        return {
            tuple(map(lambda x: shift+x+(mol*n_atoms), l)): v
            for mol in range(n_mols) for l, v in input_list.items()
            }

Molecule = collections.namedtuple('Molecule', ['pid', 'pos', 'mass', 'type'])


def combination(sig_1, eps_1, sig_2, eps_2, cr):
    if cr == 2:
        sig = 0.5*(sig_1 + sig_2)
        eps = (eps_1*eps_2)**(1.0/2.0)
    else:
        sig = (sig_1*sig_2)**(1.0/2.0)
        eps = (eps_1*eps_2)**(1.0/2.0)

    return sig, eps


def setNonbondedInteractions(system, gt, vl, lj_cutoff, tab_cutoff=None):  #NOQA
    defaults = gt.gt.defaults
    atomparams = gt.gt.atomtypes
    atomsym_atomtype = gt.atomsym_atomtype

    if tab_cutoff is None:
        tab_cutoff = lj_cutoff

    combinationrule = int(defaults['combinationrule'])
    print('Settings up LJ interactions')
    type_pairs = set()
    for type_1 in atomsym_atomtype:
        for type_2 in atomsym_atomtype:
            type_pairs.add(tuple(sorted([type_1, type_2])))

    lj_interaction = espressopp.interaction.VerletListLennardJones(vl)
    tab_interaction = espressopp.interaction.VerletListTabulated(vl)

    # Special case for MultiTabulated
    cr_multi = collections.defaultdict(list)
    cr_mix_tab = collections.defaultdict(list)
    cr_observs = {}

    print('Number of non-bonded type pairs: {}'.format(len(type_pairs)))
    for type_1, type_2 in type_pairs:
        t1 = atomsym_atomtype[type_1]
        t2 = atomsym_atomtype[type_2]
        param = gt.gt.nonbond_params.get((type_1, type_2))
        table_name = None
        cr_type = None
        cr_min, cr_max = 0, 0
        cr_total = 0
        sig_1, eps_1, sig_2, eps_2 = 0, 0, 0, 0
        sig, eps = -1, -1
        if param:
            print('Using defined non-bonded cross params')
            func = param['func']
            if func == 1:
                sig = float(param['params'][0])
                eps = float(param['params'][1])
            elif func == 8:
                table_name = 'table_{}_{}.xvg'.format(type_1, type_2)
            elif func == 1:
                sig_1, eps_1 = atomparams[type_1]['sigma'], atomparams[type_1]['epsilon']
                sig_2, eps_2 = atomparams[type_2]['sigma'], atomparams[type_2]['epsilon']
                sig, eps = combination(sig_1, eps_1, sig_2, eps_2, combinationrule)
            elif func == 9:
                tab_name = 'table_{}_{}.xvg'.format(param['params'][1], param['params'][0])
                cr_type = atomsym_atomtype[param['params'][2]]
                cr_total = int(param['params'][3])
                cr_min = float(param['params'][4])
                cr_max = float(param['params'][5])
                cr_default = bool(int(param['params'][6])) if len(param['params']) > 6 else False
                if (cr_type, cr_total) not in cr_observs:
                    cr_observs[(cr_type, cr_total)] = espressopp.analysis.ChemicalConversion(
                        system, cr_type, cr_total)
                cr_multi[(t1, t2)].append([
                    cr_observs[(cr_type, cr_total)],
                    tab_name,
                    cr_min,
                    cr_max,
                    cr_default])
            elif func == 10:
                tab1 = param['params'][0]
                tab2 = param['params'][1]
                cr_type = atomsym_atomtype[param['params'][2]]
                cr_total = int(param['params'][3])
                if (cr_type, cr_total) not in cr_observs:
                    cr_observs[(cr_type, cr_total)] = espressopp.analysis.ChemicalConversion(
                        system, cr_type, cr_total
                    )
                cr_mix_tab[(t1, t2)].append([
                    cr_observs[(cr_type, cr_total)],
                    tab1,
                    tab2])

        else:
            sig_1, eps_1 = atomparams[type_1]['sigma'], atomparams[type_1]['epsilon']
            sig_2, eps_2 = atomparams[type_2]['sigma'], atomparams[type_2]['epsilon']
            sig, eps = combination(sig_1, eps_1, sig_2, eps_2, combinationrule)
        # Standard interaction.
        if sig > 0 and eps > 0:
            print('Set lj potential {}-{}, eps={}, sig={}'.format(type_1, type_2, eps, sig))
            ljpot = espressopp.interaction.LennardJones(
                epsilon=eps, sigma=sig, cutoff=lj_cutoff)
            lj_interaction.setPotential(type1=t1, type2=t2, potential=ljpot)
        elif table_name is not None:
            print('Set tab potential {}-{}: {}'.format(type_1, type_2, table_name))
            espp_tab_name = '{}.pot'.format(table_name.replace('.xvg', ''))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(table_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(table_name, espp_tab_name)
            tab_interaction.setPotential(
                type1=t1, type2=t2,
                potential=espressopp.interaction.Tabulated(
                    itype=2, filename=espp_tab_name, cutoff=tab_cutoff))

    system.addInteraction(lj_interaction, 'lj')
    system.addInteraction(tab_interaction, 'lj-tab')

    # Conversion dependent tabulated potentials.
    if cr_multi:
        multi_tab_interaction = espressopp.interaction.VerletListMultiTabulated(vl)
        for (mt1, mt2), data in cr_multi.items():
            mp_tab = espressopp.interaction.MultiTabulated(cutoff=tab_cutoff)
            for cr_obs, tab_name, cr_min, cr_max, cr_default in data:
                espp_tab_name = '{}.pot'.format(tab_name.replace('.xvg', ''))
                if not os.path.exists(espp_tab_name):
                    print('Convert {} to {}'.format(tab_name, espp_tab_name))
                    espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
                mp_tab.register_table(espp_tab_name, 2, cr_obs, cr_min, cr_max, cr_default)
            print('Set multi tabulated potential {}-{}'.format(mt1, mt2))
            multi_tab_interaction.setPotential(
                type1=mt1, type2=mt2, potential=mp_tab)
        system.addInteraction(multi_tab_interaction, 'lj-mtab')

    # Mixed Tabulated potentials.
    if cr_mix_tab:
        mixed_tab_interaction = espressopp.interaction.VerletListMixedTabulated(vl)
        for (mt1, mt2), data in cr_mix_tab.items():
            for cr_obs, tab1, tab2 in data:
                espp_tab1_name = '{}.pot'.format(tab1.replace('.xvg', ''))
                if not os.path.exists(espp_tab1_name):
                    print('Convert {} to {}'.format(tab1, espp_tab1_name))
                    espressopp.tools.convert.gromacs.convertTable((tab1, espp_tab1_name))
                espp_tab2_name = '{}.pot'.format(tab2.replace('.xvg', ''))
                if not os.path.exists(espp_tab2_name):
                    print('Convert {} to {}'.format(tab2, espp_tab2_name))
                    espressopp.tools.convert.gromacs.convertTable((tab2, espp_tab2_name))
                print('Set mixed tabulated potential {}-{}'.format(mt1, mt2))
                mixed_tab_interaction.setPotential(
                    type1=mt1,
                    type2=mt2,
                    potential=espressopp.interaction.MixedTabulated(
                        2, espp_tab1_name, espp_tab2_name, cr_obs, cutoff=tab_cutoff
                    ))
        system.addInteraction(mixed_tab_interaction, 'lj-mix_tab')

    return cr_observs


def setBondInteractions(system, gt, only_interaction=False, name='bonds'):
    """Set bonded interactions."""

    def convert_params(func, raw_data):
        if func == 1:
            return {'K': float(raw_data[1])/2.0, 'r0': float(raw_data[0])}
        elif func == 8:
            espp_tab_name = 'table_b{}.pot'.format(int(raw_data[0]))
            tab_name = 'table_b{}.xvg'.format(int(raw_data[0]))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(tab_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
            return {'itype': 2, 'filename': espp_tab_name}
        else:
            raise RuntimeError('Unknown func type')

    func2interaction = {
        1: (espressopp.interaction.FixedPairListTypesHarmonic, espressopp.interaction.Harmonic),
        8: (espressopp.interaction.FixedPairListTypesTabulated, espressopp.interaction.Tabulated)
    }

    if not only_interaction:
        fpls = []
        bonds_by_func = collections.defaultdict(dict)

        bonds_by_params = collections.defaultdict(list)
        for b, parameters in gt.bonds.items():
            a1, a2 = map(gt.atoms.get, b)
            t1, t2 = a1['type_id'], a2['type_id']
            if parameters:
                func = int(key_tuple[0])
                params = tuple(map(float, parameters))
            else:
                bp = gt.bondparams[(t1, t2)]
                params = tuple([int(bp['func'])] + map(float, bp['params']))
                func = int(bp['func'])
            if params not in bonds_by_func[func]:
                bonds_by_func[func][params] = collections.defaultdict(list)
            bonds_by_func[func][params][(t1, t2)].append(b)



        for types, param in gt.bondparams.items():
            params = tuple([int(param['func'])] + map(float, param['params']))
            func = int(param['func'])
            print params, types
            if params not in bonds_by_func[func]:
                bonds_by_func[func][params] = {types: []}
            if types not in bonds_by_func[func][params]:
                bonds_by_func[func][params][types] = []

        bond_count = 0
        for params, bond_list in bonds_by_params.items():
            func_type = int(params[0])
            if func_type in func2interaction:
                interaction_class, potential_class = func2interaction.get(func_type)
                fpl = espressopp.FixedPairList(system.storage)
                fpls.append(fpl)
                fpl.addBonds(bond_list)
                interaction = interaction_class(system, fpl, potential_class(**convert_params(func_type, params[1:])))
                system.addInteraction(interaction, '{}_{}'.format(name, bond_count))
                print('Set static bond potential "{}_{}" ({}) func_type={} params={}'.format(
                    name, bond_count, len(bond_list), func_type, params[1:]))
                bond_count += 1
            else:
                raise RuntimerError('Bonded interaction of func {} not defined yet!'.format(func_type))
    else:
        fpls = collections.defaultdict(list)
        func2interaction = {
            1: (espressopp.interaction.FixedPairListTypesHarmonic, espressopp.interaction.Harmonic),
            8: (espressopp.interaction.FixedPairListTypesTabulated, espressopp.interaction.Tabulated)
        }
        bondtypes_by_func = collections.defaultdict(list)
        for types, param in gt.bondparams.items():
            bondtypes_by_func[param['func']].append(types)

        bond_count = 0
        for func_type, bondtypes in bondtypes_by_func.items():
            if func_type in func2interaction:
                interaction_class, potential_class = func2interaction.get(func_type)
                fpl = espressopp.FixedPairList(system.storage)
                fpls[func_type].append(fpl)
                interaction = interaction_class(system, fpl)
                for t in bondtypes:
                    param = gt.bondparams[t]['params']
                    interaction.setPotential(
                        type1=t[0], type2=t[1],
                        potential=potential_class(**convert_params(func_type, param))
                    )
                    print('Set dynamic angular potential "angle_{}" type: {}-{} with params: {}'.format(
                        bond_count, t[0], t[1], param))
                system.addInteraction(interaction, 'bond_{}'.format(bond_count))
                bond_count += 1
            else:
                raise RuntimerError('Bonded interaction of func {} not defined yet!'.format(func_type))

    return fpls


def setAngleInteractions(system, gt, only_interaction=False, name='angles'):
    """Set angle interactions."""
    def convert_params(func, raw_data):
        if func == 1:
            return {'K': float(raw_data[1]), 'theta0': float(raw_data[0])*2*math.pi/360}
        elif func == 8:
            espp_tab_name = 'table_a{}.pot'.format(int(raw_data[0]))
            tab_name = 'table_a{}.xvg'.format(int(raw_data[0]))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(tab_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
            return {'itype': 2, 'filename': espp_tab_name}
        else:
            raise RuntimeError('Unknown func type')

    if not only_interaction:
        ftls = []
        func2interaction = {
            1: (espressopp.interaction.FixedTripleListAngularHarmonic, espressopp.interaction.AngularHarmonic),
            8: (espressopp.interaction.FixedTripleListTabulatedAngular, espressopp.interaction.TabulatedAngular)
        }
        angles_by_params = collections.defaultdict(list)
        for b, parameters in gt.angles.items():
            a1, a2, a3 = map(gt.atoms.get, b)
            if parameters:
                key_tuple = tuple(map(float, parameters))
            else:
                params = gt.angleparams[(a1['type_id'], a2['type_id'], a3['type_id'])]
                key_tuple = tuple([params['func']] + map(float, params['params']))
            angles_by_params[key_tuple].append(b)

        angle_count = 0
        for params, angle_list in angles_by_params.items():
            func_type = int(params[0])
            if func_type in func2interaction:
                interaction_class, potential_class = func2interaction.get(func_type)
                ftl = espressopp.FixedTripleList(system.storage)
                ftls.append(ftl)
                ftl.addTriples(angle_list)
                interaction = interaction_class(system, ftl, potential_class(**convert_params(func_type, params[1:])))
                print('Set static angular potential "{}_{}" ({}) func_type={} params={}'.format(
                    name, angle_count, len(angle_list), func_type, params[1:]))
                system.addInteraction(interaction, '{}_{}'.format(name, angle_count))
                angle_count += 1
            else:
                raise RuntimerError('Angular interaction of func {} not defined yet!'.format(func_type))
    else:
        ftls = collections.defaultdict(list)
        func2interaction = {
            1: (espressopp.interaction.FixedTripleListTypesAngularHarmonic, espressopp.interaction.AngularHarmonic),
            8: (espressopp.interaction.FixedTripleListTypesTabulatedAngular, espressopp.interaction.TabulatedAngular)
        }
        angletypes_by_func = collections.defaultdict(list)
        for types, param in gt.angleparams.items():
            angletypes_by_func[param['func']].append(types)

        angle_count = 0
        for func_type, angletypes in angletypes_by_func.items():
            if func_type in func2interaction:
                interaction_class, potential_class = func2interaction.get(func_type)
                ftl = espressopp.FixedTripleList(system.storage)
                ftls[func_type].append(ftl)
                interaction = interaction_class(system, ftl)
                for t in angletypes:
                    param = gt.angleparams[t]['params']
                    interaction.setPotential(
                        type1=t[0], type2=t[1], type3=t[2],
                        potential=potential_class(**convert_params(func_type, param))
                    )
                    print('Set dynamic angular potential "angle_{}" type: {}-{}-{} with params: {}'.format(
                        angle_count, t[0], t[1], t[2], param))
                system.addInteraction(interaction, 'angle_{}'.format(angle_count))
                angle_count += 1
            else:
                raise RuntimerError('Angle interaction of func {} not defined yet!'.format(func_type))
    return ftls


def setDihedralInteractions(system, gt, only_interaction=False, name='dihedrals'):
    """Set dihedral interactions."""
    def convert_params(func, raw_data):
        if func == 1:
            return {'K': float(raw_data[1]), 'phi0': float(raw_data[0])*2*math.pi/360, 'multiplicity': int(raw_data[2])}
        elif func == 8:
            espp_tab_name = 'table_a{}.pot'.format(int(raw_data[0]))
            tab_name = 'table_a{}.xvg'.format(int(raw_data[0]))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(tab_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
            return {'itype': 2, 'filename': espp_tab_name}
        else:
            raise RuntimeError('Unknown func type')

    if not only_interaction:
        fqls = []
        func2interaction = {
            1: (espressopp.interaction.FixedQuadrupleListDihedralHarmonicNCos,
                espressopp.interaction.DihedralHarmonicNCos),
            8: (espressopp.interaction.FixedQuadrupleListTabulatedDihedral, espressopp.interaction.TabulatedDihedral)
        }
        dihedrals_by_params = collections.defaultdict(list)
        for b, parameters in gt.dihedrals.items():
            a1, a2, a3, a4 = map(gt.atoms.get, b)
            if parameters:
                key_tuple = tuple(map(float, parameters))
            else:
                params = gt.dihedralparams[(a1['type_id'], a2['type_id'], a3['type_id'], a4['type_id'])]
                key_tuple = tuple([params['func']] + map(float, params['params']))
            dihedrals_by_params[key_tuple].append(b)

        dihedral_count = 0
        for params, dihedral_list in dihedrals_by_params.items():
            func_type = int(params[0])
            if func_type in func2interaction:
                interaction_class, potential_class = func2interaction.get(func_type)
                fql = espressopp.FixedQuadrupleList(system.storage)
                fqls.append(fql)
                fql.addQuadruples(dihedral_list)
                interaction = interaction_class(system, fql, potential_class(**convert_params(func_type, params[1:])))
                system.addInteraction(interaction, '{}_{}'.format(name, dihedral_count))
                print('Set static dihedral potential "{}_{}" ({}) func_type={} params={}'.format(
                    name, dihedral_count, len(dihedral_list), func_type, params[1:]))
                dihedral_count += 1
            else:
                raise RuntimerError('Dihedral interaction of func {} not defined yet!'.format(func_type))
    else:
        fqls = collections.defaultdict(list)
        func2interaction = {
            1: (espressopp.interaction.FixedQuadrupleListTypesDihedralHarmonicNCos,
                espressopp.interaction.DihedralHarmonicNCos),
            8: (espressopp.interaction.FixedQuadrupleListTypesTabulatedDihedral,
                espressopp.interaction.TabulatedDihedral)
        }
        dihedraltypes_by_func = collections.defaultdict(list)
        for types, param in gt.dihedralparams.items():
            dihedraltypes_by_func[param['func']].append(types)

        dihedral_count = 0
        for func_type, dihedraltypes in dihedraltypes_by_func.items():
            if func_type in func2interaction:
                interaction_class, potential_class = func2interaction.get(func_type)
                fql = espressopp.FixedQuadrupleList(system.storage)
                fqls[func_type].append(fql)
                interaction = interaction_class(system, fql)
                for t in dihedraltypes:
                    param = gt.dihedralparams[t]['params']
                    interaction.setPotential(
                        type1=t[0], type2=t[1], type3=t[2], type4=t[3],
                        potential=potential_class(**convert_params(func_type, param))
                    )
                    print('Set dynamic dihedral potential "dihedral_{}" type: {}-{}-{}-{} with params: {}'.format(
                        dihedral_count, t[0], t[1], t[2], t[3], param))
                system.addInteraction(interaction, 'dihedral_{}'.format(dihedral_count))
                dihedral_count += 1
            else:
                raise RuntimerError('Dihedral interaction of func {} not defined yet!'.format(func_type))
    return fqls


def genParticleList(coordinate, topol):
    """Generates particle list
    Args:
        coordinate: The input coordinate file.
        topol: Topology file.
    Returns:
        List of property names and particle list.
    """
    props = ['id', 'type', 'pos', 'mass', 'q', 'res_id', 'state']
    particle_list = []

    for atom_id in sorted(coordinate.atoms):
        data = coordinate.atoms[atom_id]
        top_data = topol.atomparams['{}-{}'.format(data.chain_name, data.name)]
        particle_list.append(
            [atom_id,
             top_data['type_id'],
             espressopp.Real3D(data.position),
             top_data['mass'],
             top_data['charge'],
             data.chain_idx,
             top_data.get('state', 0)]
        )

    return props, particle_list
