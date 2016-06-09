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

import espressopp
import collections
import os
import math

import files_io

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
        self.gt = files_io.GROMACSTopologyFile(self.input_file)
        self.topol = self.gt
        self.gt.content = f.lines
        self.gt.read()

        self.master_topol = files_io.GROMACSTopologyFile(self.input_file)
        self.master_topol.read()
        if len(self.gt.molecules) > 1:
            raise RuntimeError('Multiple molecules are not supported')

        self._prepare_data()

    def _prepare_data(self):
        # Generate atom types from atom symbols.
        self.atomsym_atomtype = {}  # atom_symbol -> type_id

        self.atomparams = {}
        self.atom_id_params = {}
        self.used_atomtypes = set()
        self.used_atomnr = set()
        self.used_atomsym_atomtype = {}

        self.used_atomnr2atom_type = collections.defaultdict(list)

        combinationrule = self.topol.defaults['combinationrule']
        atype_id = 0
        for at_id in sorted(self.gt.atoms):
            at_data = self.gt.atoms[at_id]
            at_type = self.gt.atomtypes[at_data.atom_type]
            at_key = '{}-{}'.format(at_data.chain_name, at_data.name)
            if at_data.atom_type not in self.atomsym_atomtype:
                self.atomsym_atomtype[at_data.atom_type] = atype_id
                atype_id += 1
            self.atomparams[at_key] = {
                'molecule': at_data.chain_name,
                'type': at_data.atom_type,
                'sig': at_type['sigma'],
                'eps': at_type['epsilon'],
                'type_id': self.atomsym_atomtype[at_data.atom_type],
                'state': at_type.get('state', 0)
            }
            self.used_atomtypes.add(at_data.atom_type)
            self.used_atomnr.add(self.gt.atom_name2atomnr[at_data.atom_type])
            self.used_atomnr2atom_type[self.gt.atom_name2atomnr[at_data.atom_type]].append(at_data.atom_type)
            self.used_atomsym_atomtype[at_data.atom_type] = self.atomsym_atomtype[at_data.atom_type]

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

        # Update non_bonded params
        for k, v in self.topol.nonbond_params.items():
            if v['func'] == 1 and self.topol.defaults['combinationrule'] == 1:
                c6 = float(v['params'][0])
                c12 = float(v['params'][1])
                sig, eps = convertc6c12(c6, c12, combinationrule)
                v['params'][0] = sig
                v['params'][1] = eps

        # Gets the atomtypes from master_topol and set it as the used_atomtypes.
        # This is required to work with chemical reactions although it is not
        # nessecary to work with standard simulation.
        for at_name, at_data in self.master_topol.atomtypes.items():
            self.used_atomtypes.add(at_name)
            self.used_atomnr.add(self.master_topol.atom_name2atomnr[at_name])
            self.used_atomnr2atom_type[self.master_topol.atom_name2atomnr[at_name]].append(at_name)
            if at_name not in self.atomsym_atomtype:
                self.atomsym_atomtype[at_name] = atype_id
                atype_id += 1
            self.used_atomsym_atomtype[at_name] = self.atomsym_atomtype[at_name]

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
            for k, v in self.atom_id_params.items()}

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
                if i in self.used_atomnr and j in self.used_atomnr:
                    bcount += 1
                    params = self.gt.bondtypes[i][j]
                    for ti in self.used_atomnr2atom_type[i]:
                        for tj in self.used_atomnr2atom_type[j]:
                            t1 = self.atomsym_atomtype[ti]
                            t2 = self.atomsym_atomtype[tj]
                            self.bondparams[(t1, t2)] = params
        acount = 0
        for i in self.gt.angletypes:
            for j in self.gt.angletypes[i]:
                for k in self.gt.angletypes[i][j]:
                    if i in self.used_atomnr and j in self.used_atomnr and k in self.used_atomnr:
                        acount += 1
                        params = self.gt.angletypes[i][j][k]
                        for ti in self.used_atomnr2atom_type[i]:
                            for tj in self.used_atomnr2atom_type[j]:
                                for tk in self.used_atomnr2atom_type[k]:
                                    t1 = self.atomsym_atomtype[ti]
                                    t2 = self.atomsym_atomtype[tj]
                                    t3 = self.atomsym_atomtype[tk]
                                    self.angleparams[(t1, t2, t3)] = params
        dcount = 0
        for i in self.gt.dihedraltypes:
            for j in self.gt.dihedraltypes[i]:
                for k in self.gt.dihedraltypes[i][j]:
                    for l in self.gt.dihedraltypes[i][j][k]:
                        if i in self.used_atomnr and j in self.used_atomnr and k in self.used_atomnr and l in self.used_atomnr:
                            dcount += 1
                            params = self.gt.dihedraltypes[i][j][k][l]
                            for ti in self.used_atomnr2atom_type[i]:
                                for tj in self.used_atomnr2atom_type[j]:
                                    for tk in self.used_atomnr2atom_type[k]:
                                        for tl in self.used_atomnr2atom_type[l]:
                                            t1 = self.atomsym_atomtype[ti]
                                            t2 = self.atomsym_atomtype[tj]
                                            t3 = self.atomsym_atomtype[tk]
                                            t4 = self.atomsym_atomtype[tl]
                                            self.dihedralparams[(t1, t2, t3, t4)] = params

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


def setNonbondedInteractions(system, gt, vl, lj_cutoff=None, tab_cutoff=None, tables=None, cr_observs=None):  #NOQA
    """Sets the non-bonded interactions
    Args:
        system: The espressopp.System object.
        gt: The GromacsTopology object.
        vl: The VerletList object.
        lj_cutoff: The cutoff of LJ interaction.
        tab_cutoff: The cutoff of tabulated interactions.
        tables: The list of atom type names that should use tabulated potentials.
        cr_observs: The input list of ConversionObservables dictionary.

    Returns:
        The updated list of Conversion observables.
    """
    defaults = gt.gt.defaults
    atomparams = gt.gt.atomtypes
    atomsym_atomtype = gt.used_atomsym_atomtype

    if tab_cutoff is None:
        tab_cutoff = lj_cutoff

    # The list of tabuleated types.
    if tables is None:
        tables = []

    if cr_observs is None:
        cr_observs = {}

    combinationrule = int(defaults['combinationrule'])
    print('Settings up LJ interactions')
    type_pairs = set()
    for type_1 in atomsym_atomtype:
        for type_2 in atomsym_atomtype:
            type_pairs.add(tuple(sorted([type_1, type_2])))

    has_lj_interaction = False
    lj_interaction = espressopp.interaction.VerletListLennardJones(vl)
    has_tab_interaction = False
    tab_interaction = espressopp.interaction.VerletListTabulated(vl)

    # Special case for MultiTabulated
    cr_multi = collections.defaultdict(list)
    cr_mix_tab = collections.defaultdict(list)

    print('Number of non-bonded type pairs: {}'.format(len(type_pairs)))
    defined_types = set()
    all_type_pairs = set()
    for type_1, type_2 in type_pairs:
        t1 = atomsym_atomtype[type_1]
        t2 = atomsym_atomtype[type_2]
        all_type_pairs.add((t1, t2))
        param = gt.gt.nonbond_params.get((type_1, type_2))
        table_name = None
        sig, eps = -1, -1
        if param:
            func = param['func']
            print('Using defined non-bonded cross params {} {}'.format(func, param['params']))
            if func == 1:
                if param['params']:
                    sig = float(param['params'][0])
                    eps = float(param['params'][1])
                else:
                    sig_1, eps_1 = atomparams[type_1]['sigma'], atomparams[type_1]['epsilon']
                    sig_2, eps_2 = atomparams[type_2]['sigma'], atomparams[type_2]['epsilon']
                    sig, eps = combination(sig_1, eps_1, sig_2, eps_2, combinationrule)
            elif func == 8:
                table_name = 'table_{}_{}.xvg'.format(type_1, type_2)
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
        elif type_1 in tables and type_2 in tables:
            table_name = 'table_{}_{}.xvg'.format(type_1, type_2)
        else:
            sig_1, eps_1 = atomparams[type_1]['sigma'], atomparams[type_1]['epsilon']
            sig_2, eps_2 = atomparams[type_2]['sigma'], atomparams[type_2]['epsilon']
            sig, eps = combination(sig_1, eps_1, sig_2, eps_2, combinationrule)

        # Standard interaction.
        if table_name is not None:
            print('Set tab potential {}-{}: {}'.format(type_1, type_2, table_name))
            espp_tab_name = '{}.pot'.format(table_name.replace('.xvg', ''))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(table_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(table_name, espp_tab_name)
            has_tab_interaction = True
            tab_interaction.setPotential(
                type1=t1, type2=t2,
                potential=espressopp.interaction.Tabulated(
                    itype=2, filename=espp_tab_name, cutoff=tab_cutoff))
            defined_types.add((t1, t2))
        elif eps > 0.0  and sig > 0.0:
            print('Set LJ potential {}-{}, eps={}, sig={}'.format(type_1, type_2, eps, sig))
            ljpot = espressopp.interaction.LennardJones(
                epsilon=eps, sigma=sig, cutoff=lj_cutoff)
            lj_interaction.setPotential(type1=t1, type2=t2, potential=ljpot)
            has_lj_interaction = True
            defined_types.add((t1, t2))

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
            defined_types.add((mt1, mt2))
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
                defined_types.add((mt1, mt2))
        system.addInteraction(mixed_tab_interaction, 'lj-mix_tab')

    if has_lj_interaction:
        print('Adding lj interaction')
        system.addInteraction(lj_interaction, 'lj')
    if has_tab_interaction:
        print('Adding lj-tab interaction')
        system.addInteraction(tab_interaction, 'lj-tab')

    return cr_observs


def set_bonded_interactions(system, gt, name='bonds'):
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
            raise RuntimeError('Unknown func type {}'.format(func))

    func2interaction_dynamic = {
        1: (espressopp.interaction.FixedPairListTypesHarmonic, espressopp.interaction.Harmonic),
        8: (espressopp.interaction.FixedPairListTypesTabulated, espressopp.interaction.Tabulated)
    }

    func2interaction_static = {
        1: (espressopp.interaction.FixedPairListHarmonic, espressopp.interaction.Harmonic),
        8: (espressopp.interaction.FixedPairListTabulated, espressopp.interaction.Tabulated)
    }

    # Sort existing bond lists by functional type.
    bonds_by_func = collections.defaultdict(dict)
    dynamics_bonds_by_func = collections.defaultdict(list)
    for b, parameters in gt.bonds.items():
        ptypes = map(lambda x: gt.atoms[x]['type_id'], b)
        if parameters:
            func = int(parameters[0])
            params = tuple(map(float, parameters[1:]))
            if params not in bonds_by_func[func]:
                bonds_by_func[func][params] = []
            bonds_by_func[func][params].append(b)
        else:
            pparams = gt.bondparams[tuple(ptypes)]
            if not pparams:
                pparams = gt.bonds[tuple(reversed(ptypes))]
            func = int(pparams['func'])
            dynamics_bonds_by_func[func].append(b)

    # Set first static bonds, those one that has explicitly the parameters
    static_fpls = []
    bond_count = 0
    for func in bonds_by_func:
        interaction_class, potential_class = func2interaction_static.get(func)
        for params, b_list in bonds_by_func[func].items():
            fpl = espressopp.FixedPairList(system.storage)
            fpl.addBonds(b_list)
            static_fpls.append(fpl)
            interaction = interaction_class(system, fpl, potential_class(**convert_params(func, params)))
            system.addInteraction(interaction, '{}_{}'.format(name, bond_count))
            #print('Set static bond potential "{}_{}" ({}) func_type={} params={}'.format(
            #    name, bond_count, len(b_list), func, params[1:]))
            bond_count += 1

    dynamics_fpls = collections.defaultdict(dict)
    bondparams_func = collections.defaultdict(list)
    for pt, p in gt.bondparams.items():
        bondparams_func[p['func']].append((pt, p))
        if p['func'] not in dynamics_bonds_by_func:
            dynamics_bonds_by_func[p['func']] = []

    for func, b_list in dynamics_bonds_by_func.items():
        fpl = espressopp.FixedPairList(system.storage)
        fpl.addBonds(b_list)
        interaction_class, potential_class = func2interaction_dynamic.get(func)
        interaction = interaction_class(system, fpl)
        for t, params in bondparams_func[func]:
            interaction.setPotential(
                type1=t[0], type2=t[1],
                potential=potential_class(**convert_params(func, params['params']))
            )
        system.addInteraction(interaction, 'bond_{}'.format(bond_count))
        bond_count += 1
        dynamics_fpls[func] = fpl

    print('Set up bond interactions')
    return dynamics_fpls, static_fpls


def set_angle_interactions(system, gt, name='angles'):
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

    func2interaction_dynamic = {
        1: (espressopp.interaction.FixedTripleListTypesAngularHarmonic, espressopp.interaction.AngularHarmonic),
        8: (espressopp.interaction.FixedTripleListTypesTabulatedAngular, espressopp.interaction.TabulatedAngular)
    }

    func2interaction_static = {
        1: (espressopp.interaction.FixedTripleListAngularHarmonic, espressopp.interaction.AngularHarmonic),
        8: (espressopp.interaction.FixedTripleListTabulatedAngular, espressopp.interaction.TabulatedAngular)
    }

    # Sort existing angle lists by functional type.
    angles_by_func = collections.defaultdict(dict)
    dynamics_angles_by_func = collections.defaultdict(list)
    for b, parameters in gt.angles.items():
        ptypes = map(lambda x: gt.atoms[x]['type_id'], b)
        if parameters:
            func = int(parameters[0])
            params = tuple(map(float, parameters[1:]))
            if params not in angles_by_func[func]:
                angles_by_func[func][params] = []
            angles_by_func[func][params].append(b)
        else:
            pparams = gt.angleparams[tuple(ptypes)]
            if not pparams:
                pparams = gt.angles[tuple(reversed(ptypes))]
            func = int(pparams['func'])
            dynamics_angles_by_func[func].append(b)

    # Set first static angles, those one that has explicitly the parameters
    static_ftls = []
    angle_count = 0
    for func in angles_by_func:
        interaction_class, potential_class = func2interaction_static.get(func)
        for params, b_list in angles_by_func[func].items():
            ftl = espressopp.FixedTripleList(system.storage)
            ftl.addTriples(b_list)
            static_ftls.append(ftl)
            interaction = interaction_class(system, ftl, potential_class(**convert_params(func, params)))
            system.addInteraction(interaction, '{}_{}'.format(name, angle_count))
            angle_count += 1

    dynamics_ftls = collections.defaultdict(dict)
    angleparams_func = collections.defaultdict(list)
    for pt, p in gt.angleparams.items():
        angleparams_func[p['func']].append((pt, p))
        if p['func'] not in dynamics_angles_by_func:
            dynamics_angles_by_func[p['func']] = []

    for func, b_list in dynamics_angles_by_func.items():
        ftl = espressopp.FixedTripleList(system.storage)
        ftl.addTriples(b_list)
        interaction_class, potential_class = func2interaction_dynamic.get(func)
        interaction = interaction_class(system, ftl)
        for t, params in angleparams_func[func]:
            print t, params
            interaction.setPotential(
                type1=t[0], type2=t[1], type3=t[2],
                potential=potential_class(**convert_params(func, params['params']))
            )
        system.addInteraction(interaction, 'angle_{}'.format(angle_count))
        angle_count += 1
        dynamics_ftls[func] = ftl
    print('Set up angle interactions')
    return dynamics_ftls, static_ftls


def set_dihedral_interactions(system, gt, name='dihedrals'):
    """Set dihedral interactions."""
    def convert_params(func, raw_data):
        if func == 1:
            return {'K': float(raw_data[1]), 'phi0': float(raw_data[0])*2*math.pi/360, 'multiplicity': int(raw_data[2])}
        elif func == 3:
            t = map(float, raw_data[1:])
            return {'K{}'.format(i): v for i, v in enumerate(t[1:])}
        elif func == 8:
            espp_tab_name = 'table_d{}.pot'.format(int(raw_data[0]))
            tab_name = 'table_d{}.xvg'.format(int(raw_data[0]))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(tab_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
            return {'itype': 2, 'filename': espp_tab_name}
        else:
            raise RuntimeError('Unknown func type')

    func2interaction_dynamic = {
        1: (espressopp.interaction.FixedQuadrupleListTypesDihedralHarmonicNCos,
            espressopp.interaction.DihedralHarmonicNCos),
        3: (espressopp.interaction.FixedQuadrupleListTypesDihedralRB,
            espressopp.interaction.DihedralRB),
        8: (espressopp.interaction.FixedQuadrupleListTypesTabulatedDihedral,
            espressopp.interaction.TabulatedDihedral)
    }

    func2interaction_static = {
        1: (espressopp.interaction.FixedQuadrupleListDihedralHarmonicNCos,
            espressopp.interaction.DihedralHarmonicNCos),
        3: (espressopp.interaction.FixedQuadrupleListDihedralRB,
            espressopp.interaction.DihedralRB),
        8: (espressopp.interaction.FixedQuadrupleListTabulatedDihedral,
            espressopp.interaction.TabulatedDihedral)
    }

    # Sort existing dihedral lists by functional type.
    dihedrals_by_func = collections.defaultdict(dict)
    dynamics_dihedrals_by_func = collections.defaultdict(list)
    for b, parameters in gt.dihedrals.items():
        ptypes = map(lambda x: gt.atoms[x]['type_id'], b)
        if parameters:
            func = int(parameters[0])
            params = tuple(map(float, parameters[1:]))
            if params not in dihedrals_by_func[func]:
                dihedrals_by_func[func][params] = []
            dihedrals_by_func[func][params].append(b)
        else:
            pparams = gt.dihedralparams[tuple(ptypes)]
            if not pparams:
                pparams = gt.dihedrals[tuple(reversed(ptypes))]
            func = int(pparams['func'])
            dynamics_dihedrals_by_func[func].append(b)

    # Set first static dihedrals, those one that has explicitly the parameters
    static_fqls = []
    dihedral_count = 0
    for func in dihedrals_by_func:
        interaction_class, potential_class = func2interaction_static.get(func)
        for params, b_list in dihedrals_by_func[func].items():
            fql = espressopp.FixedQuadrupleList(system.storage)
            fql.addQuadruples(b_list)
            static_fqls.append(fql)
            interaction = interaction_class(system, fql, potential_class(**convert_params(func, params)))
            system.addInteraction(interaction, '{}_{}'.format(name, dihedral_count))
            dihedral_count += 1

    dynamics_fqls = collections.defaultdict(dict)
    dihedralparams_func = collections.defaultdict(list)
    for pt, p in gt.dihedralparams.items():
        dihedralparams_func[p['func']].append((pt, p))
        if p['func'] not in dynamics_dihedrals_by_func:
            dynamics_dihedrals_by_func[p['func']] = []

    for func, b_list in dynamics_dihedrals_by_func.items():
        fql = espressopp.FixedQuadrupleList(system.storage)
        fql.addQuadruples(b_list)
        interaction_class, potential_class = func2interaction_dynamic.get(func)
        interaction = interaction_class(system, fql)
        for t, params in dihedralparams_func[func]:
            interaction.setPotential(
                type1=t[0], type2=t[1], type3=t[2], type4=t[3],
                potential=potential_class(**convert_params(func, params['params']))
            )
        system.addInteraction(interaction, '{}_{}'.format(name, dihedral_count))
        dihedral_count += 1
        dynamics_fqls[func] = fql

    print('Set up dihedral interactions')
    return dynamics_fqls, static_fqls


def set_pair_interactions(system, gt, args):
    fudgeLJ = gt.gt.defaults.get('fudgeLJ', 1.0)
    fudgeQQ = gt.gt.defaults.get('fudgeQQ', 1.0)

    # Set LJ pair interactions. If parameters are set then add to static list, otherwise to dynamic.

    # Sort existing pair lists by functional type.
    static_pairs = {}
    dynamics_pairs = []
    for b, parameters in gt.pairs.items():
        func = parameters[0]
        params = tuple(map(float, parameters[1:]))
        if params:
            if params not in static_pairs:
                static_pairs[params] = []
            static_pairs[params].append(b)
        else:
            dynamics_pairs.append(b)

    # Set first static pairs, those one that has explicitly the parameters
    combinationrule = int(gt.gt.defaults['combinationrule'])
    static_fpls = []
    pair_count = 0
    for params, b_list in static_pairs.items():
        fpl = espressopp.FixedPairList(system.storage)
        fpl.addBonds(b_list)
        static_fpls.append(fpl)
        sig, eps = combination(float(params[0]), float(params[1]), combinationrule)
        interaction = espressopp.interaction.FixedPairListLennardJones(
            system, fpl,
            espressopp.interaction.LennardJones(
                epsilon=eps,
                sigma=sig,
                cutoff=args.lj_cutoff))
        system.addInteraction(interaction, 'lj14_{}'.format(pair_count))
        static_fpls.append(fpl)
        pair_count += 1

    # Set dynamic interactions.
    dfpl = espressopp.FixedPairList(system.storage)
    dfpl.addBonds(dynamics_pairs)
    interaction = espressopp.interaction.FixedPairListTypesLennardJones(system, dfpl)
    type_pairs = set()
    for type_1 in gt.used_atomsym_atomtype:
        for type_2 in gt.used_atomsym_atomtype:
            type_pairs.add(tuple(sorted([type_1, type_2])))
    atomsym_atomtype = gt.used_atomsym_atomtype
    atomparams = gt.gt.atomtypes
    for type_1, type_2 in type_pairs:
        t1 = atomsym_atomtype[type_1]
        t2 = atomsym_atomtype[type_2]
        sig_1, eps_1 = atomparams[type_1]['sigma'], atomparams[type_1]['epsilon']
        sig_2, eps_2 = atomparams[type_2]['sigma'], atomparams[type_2]['epsilon']
        sig, eps = combination(sig_1, eps_1, sig_2, eps_2, combinationrule)
        interaction.setPotential(
            type1=t1, type2=t2,
            potential=espressopp.interaction.LennardJones(
                sigma=sig, epsilon=fudgeLJ*eps, cutoff=args.lj_cutoff)
        )
    system.addInteraction(interaction, 'lj14_{}'.format(pair_count))

    print('Set up 1-4 pair interactions')
    # Set coulombic pair interaction
    prefQQ = 138.935485 * fudgeQQ
    qq_count = 0
    if args.coulomb_cutoff > 0.0 and prefQQ > 0.0:
        potQQ = espressopp.interaction.CoulombTruncated(
            prefactor=prefQQ, cutoff=args.coulomb_cutoff)
        for fpl in static_fpls:
            interaction = espressopp.interaction.FixedPairListCoulombTruncated(
                system, fpl,
                potential=potQQ)
            #system.addInteraction(interaction, 'coulomb_14_{}'.format(qq_count))
            qq_count += 1
        interaction = espressopp.interaction.FixedPairListTypesCoulombTruncated(
            system, dfpl)
        for type_1, type_2 in type_pairs:
            t1 = atomsym_atomtype[type_1]
            t2 = atomsym_atomtype[type_2]
            interaction.setPotential(type1=t1, type2=t2, potential=potQQ)
        system.addInteraction(interaction, 'coulomb_14_{}'.format(qq_count))

    return dfpl, static_fpls


def set_coulomb_interactions(system, gt, args):
    return None


def gen_particle_list(coordinate, topol):
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
