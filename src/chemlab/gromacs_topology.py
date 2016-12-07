#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
#
#  This file is part of ChemLab.
#
#  Large part of the file is taken from ESPResSo++.
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
    def __init__(self, input_topol, generate_exclusions=True):
        self.input_file = input_topol
        self.content = None
        self.data = {}

        self.atomsym_atomtype = {}
        self.atomtype_atomsym = {}

        self.atoms = {}

        self.bonds = {}
        self.angles = {}
        self.dihedrals = {}
        self.pairs = {}

        self.bondparams = {}
        self.angleparams = {}
        self.dihedralparams = {}

        self.generate_exclusions = generate_exclusions

    def read(self):
        fb = FileBuffer()
        defines = {}
        FillFileBuffer(self.input_file, fb, defines=defines)
        f = PostProcessFileBuffer(fb, defines)
        self.gt = files_io.GROMACSTopologyFile(self.input_file)
        self.topol = self.gt
        self.gt.content = f.lines
        self.gt.read()

        print('Reading master topology file...')
        self.master_topol = files_io.GROMACSTopologyFile(self.input_file)
        self.master_topol.read()
        #if len(self.gt.molecules) > 1:
        #    raise RuntimeError('Multiple molecules are not supported')
        print('Preparing topology data for simulation...')
        self._prepare_data()

    def add_new_atomtype(self, atype_id, atype_name, is_used=False):
        """Adds new atomtype to the data structures.

        Args:
            atype_id: Atom type id.
            atype_name: Atom type name.
            is_used: Is used in simulation?
        """
        self.atomtype_atomsym[atype_id] = atype_name
        self.atomsym_atomtype[atype_name] = atype_id
        if is_used:
            self.used_atomsym_atomtype[atype_name] = atype_id

    def _prepare_data(self):
        # Generate atom types from atom symbols.
        self.atomsym_atomtype = {}  # atom_symbol -> type_id

        self.atomparams = {}
        self.used_atomtypes = set()
        self.used_atomnr = set()
        self.used_atomsym_atomtype = {}

        self.used_atomnr2atom_type = collections.defaultdict(set)

        combinationrule = self.topol.defaults['combinationrule']
        atype_id = 0
        atom_id_offset = 0
        for molecule_name, molecule_numbers in self.gt.molecules:
            print('Building {} with {} molecules'.format(molecule_name, molecule_numbers))
            atom_id_params = {}
            for at_id in sorted(self.gt.molecules_data[molecule_name]['atoms']):
                at_data = self.gt.molecules_data[molecule_name]['atoms'][at_id]
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
                    'state': at_type.get('state', 0),
                    'charge': 0.0,
                    'mass': 0.0,
                    'molecule_name': at_data.molecule_name,
                    'name': at_data.name,
                    'cgnr': at_data.cgnr,
                    'chain_idx': at_data.chain_idx,
                    'chain_name': at_data.chain_name}
                self.used_atomtypes.add(at_data.atom_type)
                self.used_atomnr.add(self.gt.atom_name2atomnr[at_data.atom_type])
                self.used_atomnr2atom_type[self.gt.atom_name2atomnr[at_data.atom_type]].add(at_data.atom_type)
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
                atom_id_params[at_id] = self.atomparams[at_key]
            # Replicate for this molecule, parameters
            n_atoms = len(self.gt.molecules_data[molecule_name]['atoms'])
            self.atoms.update({
                atom_id_offset + k + (mol * n_atoms): v for mol in range(molecule_numbers)
                for k, v in atom_id_params.items()})
            atom_id_offset += molecule_numbers

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
        # necessary to work with standard simulation.
        for at_name, at_data in self.master_topol.atomtypes.items():
            self.used_atomtypes.add(at_name)
            self.used_atomnr.add(self.master_topol.atom_name2atomnr[at_name])
            self.used_atomnr2atom_type[self.master_topol.atom_name2atomnr[at_name]].add(at_name)
            if at_name not in self.atomsym_atomtype:
                self.atomsym_atomtype[at_name] = atype_id
                atype_id += 1
            self.used_atomsym_atomtype[at_name] = self.atomsym_atomtype[at_name]

        self.atomtype_atomsym = {v: k for k, v in self.atomsym_atomtype.items()}

        self._prepare_bondedparams()
        self._prepare_bondedlists()
        if self.generate_exclusions:
            self._prepare_exclusionlists()

    def _prepare_bondedlists(self):
        """Replicate bonded lists."""
        bonded_lists = [
            ('bonds', self.bonds),
            ('angles', self.angles),
            ('dihedrals', self.dihedrals),
            ('pairs', self.pairs)]

        atom_id_offset = 0
        for molecule_name, n_mols in self.gt.molecules:
            print('Replicate bonded lists: {}'.format(molecule_name))
            n_atoms = len(self.gt.molecules_data[molecule_name]['atoms'])
            for list_name, list_ptr in bonded_lists:
                if list_name in self.gt.molecules_data[molecule_name]:
                    list_ptr.update(self._replicate_lists(
                        n_mols,
                        n_atoms,
                        self.gt.molecules_data[molecule_name][list_name],
                        atom_id_offset
                    ))
            atom_id_offset += n_mols

    def _prepare_exclusionlists(self):
        self.exclusions = {tuple(sorted(x)) for x in self.bonds.keys()}
        # Generate remaining exclusions, based on the nrexcl.
        atom_id_offset = 0
        for molecule_name, n_mols in self.gt.molecules:
            print('Building exclusion list for: {}'.format(molecule_name))
            n_atoms = len(self.gt.molecules_data[molecule_name]['atoms'])
            nrexcl = self.gt.moleculetype[molecule_name]  # This many bonds to be removed.
            if 'bonds' in self.gt.molecules_data[molecule_name]:
                mol_excl = self._generate_exclusions(self.gt.molecules_data[molecule_name]['bonds'], nrexcl)
                self.exclusions.update({
                    tuple(sorted(map(lambda x: atom_id_offset + x + (mol*n_atoms), l)))
                    for mol in range(n_mols) for l in mol_excl
                })
            else:
                print('Molecule {} does not have bonds, no exclusions needed'.format(molecule_name))
            atom_id_offset += n_mols  # Still the offset is neighter the exclusions are generated or not

    def _generate_exclusions(self, bond_list, nrexcl):
        """Generates exclusion list for single molecule to be replicated."""
        # Build  neighbour ist, simple who is neighbour of whom.
        exclusion_list = set(bond_list)

        class Node:
            def __init__(self, pid):
                self.pid = pid
                self.neighbours = []

            def add_neighbours(self, nb):
                self.neighbours.append(nb)

        def get_node_by_id(id, nodes):
            node_list = [n for n in nodes if n.pid == id]
            if len(node_list) > 1:
                print "Error: duplicate nodes", id
                exit()
            elif len(node_list) == 0:
                return None
            return node_list[0]

        def get_next_neighbour(root, nr_neighbours, neighbours, visited_nodes):
            if nr_neighbours == 0:
                return neighbours

            # avoid going back the same path
            visited_nodes.append(root)

            # Loop over next neighbours and add them to the neighbours list
            # Recursively call the function with numberNeighbours-1
            for n in root.neighbours:
                if not n in visited_nodes:
                    if n not in neighbours: neighbours.append(n)  # avoid double counting in rings
                    get_next_neighbour(n, nr_neighbours - 1, neighbours, visited_nodes)

        nodes = []
        # make a Node object for each atom involved in bonds
        for bids in bond_list:
            for i in bids:
                if get_node_by_id(i, nodes) is None:
                    n = Node(i)
                    nodes.append(n)

        # find the next neighbours for each node and append them
        for b in bond_list:
            permutations = [(b[0], b[1]), (b[1], b[0])]
            for p in permutations:
                n = get_node_by_id(p[0], nodes)
                nn = get_node_by_id(p[1], nodes)
                n.add_neighbours(nn)

        # seraches for nrexcl next neighbours
        for n in nodes:
            neighbours = []
            get_next_neighbour(n, nrexcl, neighbours, visited_nodes=[])
            for nb in neighbours:
                # check if the permutation is already in the exclusion list
                exclusion_list.add(tuple(sorted((n.pid, nb.pid))))

        return exclusion_list

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
                            self.bondparams[tuple(sorted([t1, t2]))] = params
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
                                    ptypes = (t1, t2, t3)
                                    if t1 > t3:
                                        ptypes = (t3, t2, t1)
                                    self.angleparams[ptypes] = params

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
                                            ptypes = (t1, t2, t3, t4)
                                            if t4 > t1:
                                                ptypes = (t4, t3, t2, t1)
                                            self.dihedralparams[ptypes] = params

    def _replicate_lists(self, n_mols, n_atoms, input_list, shift=0):
        """Replicate input lists by n_mols, where every mol contains n_atoms.

        Args:
            n_mols: Number of mols.
            n_atoms: Number of atoms in single mol.
            input_list: The input list to replicate.
            shift: Shift the atom ids by this value.

        Returns:
             The replicated list.
        """
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


def set_nonbonded_interactions(system, gt, vl, lj_cutoff=None, tab_cutoff=None, tables=None, cr_observs=None):  #NOQA
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
        The list of ParticlePairScaling objects.
    """
    defaults = gt.gt.defaults
    atomparams = gt.gt.atomtypes
    atomsym_atomtype = gt.used_atomsym_atomtype
    atomtype_atomsym = gt.atomsym_atomtype

    if tab_cutoff is None:
        tab_cutoff = lj_cutoff

    # The list of tabuleated types.
    if tables is None:
        tables = []

    if cr_observs is None:
        cr_observs = {}

    particle_pair_scales = []

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
    has_tab_capped_interaction = False
    tab_capped_interaction = espressopp.interaction.VerletListTabulatedCapped(vl)

    # Special case for MultiTabulated
    cr_multi = collections.defaultdict(list)
    cr_mix_tab = collections.defaultdict(list)
    tab_scaled = {} # increment -> (t1, t2)
    dynamic_interactions = collections.defaultdict(dict)

    Func10 = collections.namedtuple('Func10', ['cr_observers', 'tab1', 'tab2'])
    Func12 = collections.namedtuple('Func12', ['mix_value', 'tab1', 'tab2'])

    print('Number of non-bonded type pairs: {}'.format(len(type_pairs)))
    defined_types = set()
    all_type_pairs = set()
    for type_1, type_2 in type_pairs:
        t1 = atomsym_atomtype[type_1]
        t2 = atomsym_atomtype[type_2]
        all_type_pairs.add((t1, t2))
        param = gt.gt.nonbond_params.get((type_1, type_2))
        table_name = None
        table_cap = None
        sig, eps = -1, -1
        if param:
            func = param['func']
            #print('Using defined non-bonded cross params {} {}'.format(func, param['params']))
            if func == 1:
                if param['params']:
                    sig = float(param['params'][0])
                    eps = float(param['params'][1])
                else:
                    sig_1, eps_1 = atomparams[type_1]['sigma'], atomparams[type_1]['epsilon']
                    sig_2, eps_2 = atomparams[type_2]['sigma'], atomparams[type_2]['epsilon']
                    sig, eps = combination(sig_1, eps_1, sig_2, eps_2, combinationrule)
            elif func == 8:
                if param['params']:
                    table_name = param['params'][0]
                else:
                    table_name = 'table_{}_{}.xvg'.format(type_1, type_2)
            elif func == 9:
                tab_name = param['params'][0]
                cr_type = atomsym_atomtype[param['params'][1]]
                cr_total = int(param['params'][2])
                cr_min = float(param['params'][3])
                cr_max = float(param['params'][4])
                cr_default = bool(int(param['params'][5])) if len(param['params']) > 5 else False
                if (cr_type, cr_total) not in cr_observs:
                    cr_observs[(cr_type, cr_total)] = espressopp.analysis.ChemicalConversion(
                        system, cr_type, cr_total)
                cr_multi[(t1, t2)].append([
                    cr_observs[(cr_type, cr_total)],
                    tab_name,
                    cr_min,
                    cr_max,
                    cr_default,
                    cr_type,
                    cr_total])
            elif func == 10:
                tab1 = param['params'][0]
                tab2 = param['params'][1]
                cr_type = atomsym_atomtype[param['params'][2]]
                cr_total = int(param['params'][3])
                if (cr_type, cr_total) not in cr_observs:
                    cr_observs[(cr_type, cr_total)] = espressopp.analysis.ChemicalConversion(
                        system, cr_type, cr_total
                    )
                cr_mix_tab[(t1, t2)].append(Func10(cr_observs[(cr_type, cr_total)], tab1, tab2))
            elif func == 11:  # Special case, dynamic interactions with tables
                max_force = -1
                if param['params']:
                    tn = param['params'][0]
                    if len(param['params']) == 2:
                        max_force = float(param['params'][1])
                else:
                    tn = 'table_{}_{}.xvg'.format(type_1, type_2)
                sig, eps = 0.0, 0.0
                dynamic_interactions[func][(t1, t2)] = (tn, max_force)
            elif func == 12:  # MixedTabulated with static x
                tab1 = param['params'][0]
                tab2 = param['params'][1]
                mix_value = float(param['params'][2])
                cr_mix_tab[(t1, t2)].append(Func12(mix_value, tab1, tab2))
                print t1, t2, param
            elif func == 13:
                table_name = param['params'][0]
                table_cap = float(param['params'][1])
            elif func == 14:  # Special case of tabulated potential that has additional scaling function.
                tab1 = param['params'][0]
                scale_increment = float(param['params'][1])
                if len(param['params']) == 3:
                    max_force = float(param['params'][2])
                else:
                    max_force = -1
                if scale_increment not in tab_scaled:
                    tab_scaled[scale_increment] = {}
                tab_scaled[scale_increment][(t1, t2)] = (tab1, max_force)

        elif type_1 in tables and type_2 in tables:
            table_name = 'table_{}_{}.xvg'.format(type_1, type_2)
        else:
            sig_1, eps_1 = atomparams[type_1]['sigma'], atomparams[type_1]['epsilon']
            sig_2, eps_2 = atomparams[type_2]['sigma'], atomparams[type_2]['epsilon']
            sig, eps = combination(sig_1, eps_1, sig_2, eps_2, combinationrule)

        # Standard interaction.
        if table_name is not None and table_cap is not None:
            print('Set tab potential {}-{}: {} caprad={}'.format(type_1, type_2, table_name, table_cap))
            espp_tab_name = '{}.pot'.format(table_name.replace('.xvg', ''))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(table_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(table_name, espp_tab_name)
            has_tab_capped_interaction = True
            tab_capped_interaction.setPotential(
                type1=t1, type2=t2,
                potential=espressopp.interaction.TabulatedCapped(
                    itype=1, filename=espp_tab_name, cutoff=tab_cutoff, caprad=table_cap))
            defined_types.add((t1, t2))
        elif table_name is not None:
            print('Set tab potential {}-{}: {}'.format(type_1, type_2, table_name))
            espp_tab_name = '{}.pot'.format(table_name.replace('.xvg', ''))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(table_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(table_name, espp_tab_name)
            has_tab_interaction = True
            tab_interaction.setPotential(
                type1=t1, type2=t2,
                potential=espressopp.interaction.Tabulated(
                    itype=1, filename=espp_tab_name, cutoff=tab_cutoff))
            defined_types.add((t1, t2))
        elif sig > 0.0:
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
            for cr_obs, tab_name, cr_min, cr_max, cr_default, cr_type, cr_total in data:
                espp_tab_name = '{}.pot'.format(tab_name.replace('.xvg', ''))
                if not os.path.exists(espp_tab_name):
                    print('Convert {} to {}'.format(tab_name, espp_tab_name))
                    espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
                mp_tab.register_table(espp_tab_name, 1, cr_obs, cr_min, cr_max, cr_default)
                print('Set multi tabulated potential {}-{} range=[{},{}) type obs: {}'.format(
                    mt1, mt2, cr_min, cr_max, cr_type))
            multi_tab_interaction.setPotential(
                type1=mt1, type2=mt2, potential=mp_tab)
            defined_types.add((mt1, mt2))
        system.addInteraction(multi_tab_interaction, 'lj-mtab')

    # Mixed Tabulated potentials.
    if cr_mix_tab:
        mixed_tab_interaction = espressopp.interaction.VerletListMixedTabulated(vl)
        for (mt1, mt2), data in cr_mix_tab.items():
            for mix_func_data in data:
                cr_obs, tab1, tab2 = mix_func_data
                espp_tab1_name = '{}.pot'.format(tab1.replace('.xvg', ''))
                if not os.path.exists(espp_tab1_name):
                    print('Convert {} to {}'.format(tab1, espp_tab1_name))
                    espressopp.tools.convert.gromacs.convertTable((tab1, espp_tab1_name))
                espp_tab2_name = '{}.pot'.format(tab2.replace('.xvg', ''))
                if not os.path.exists(espp_tab2_name):
                    print('Convert {} to {}'.format(tab2, espp_tab2_name))
                    espressopp.tools.convert.gromacs.convertTable((tab2, espp_tab2_name))
                if isinstance(mix_func_data, Func10):
                    print('Set mixed tabulated potential {}-{} with conversion observable (U=x*{} + (1-x)*{})'.format(
                        mt1, mt2, espp_tab1_name, espp_tab2_name))
                    mixed_tab_interaction.setPotential(
                        type1=mt1,
                        type2=mt2,
                        potential=espressopp.interaction.MixedTabulated(
                            1, espp_tab1_name, espp_tab2_name, cr_obs, cutoff=tab_cutoff
                        ))
                elif isinstance(mix_func_data, Func12):
                    print('Set mixed tabulated potential {}-{} with static scaling x={x} (U={x}*{table1}+(1-{x})*{table2})'.format(
                        mt1, mt2, x=cr_obs, table1=espp_tab1_name, table2=espp_tab2_name))
                    mixed_tab_interaction.setPotential(
                        type1=mt1,
                        type2=mt2,
                        potential=espressopp.interaction.MixedTabulated(
                            1, table1=espp_tab1_name, table2=espp_tab2_name, mix_value=cr_obs))
                else:
                    raise RuntimeError('Wrong type of data: {}'.format(type(data)))
                defined_types.add((mt1, mt2))
        system.addInteraction(mixed_tab_interaction, 'lj-mix_tab')

    if tab_scaled:
        for scale_increment, data_list in tab_scaled.items():
            particle_pair_scale_map = espressopp.esutil.ParticlePairScaling(
                0.0, scale_increment, vl, system.integrator)
            particle_pair_scales.append(particle_pair_scale_map)
            max_forces_group = collections.defaultdict(dict)
            for (t1, t2), (tab_name, max_force) in data_list.items():
                max_forces_group[max_force][(t1, t2)] = tab_name
            bn = 0
            for max_force, data in max_forces_group.items():
                tab_scaled_interaction = espressopp.interaction.VerletListScaleTabulated(vl, particle_pair_scale_map)
                for (t1, t2), tab_name in data.items():
                    espp_tab_name = '{}.pot'.format(tab_name.replace('.xvg', ''))
                    if not os.path.exists(espp_tab_name):
                        print('Convert {} to {}'.format(tab_name, espp_tab_name))
                        espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
                    print('Set scaled potential {}-{} (max force: {} scale_increment: {})'.format(
                        t1, t2, max_force, scale_increment))
                    tab_scaled_interaction.setPotential(
                        type1=t1,
                        type2=t2,
                        potential=espressopp.interaction.Tabulated(1, espp_tab_name, cutoff=tab_cutoff))
                if max_force != -1:
                    tab_scaled_interaction.setMaxForce(max_force)
                system.addInteraction(tab_scaled_interaction, 'tab-scaled_{}'.format(bn))
                bn += 1

    if dynamic_interactions:
        print('Set up DynamicResolution non-bonded interactions')
        for func, data_list in dynamic_interactions.items():
            if func == 11:
                max_forces_group = collections.defaultdict(dict)
                for (t1, t2), (tab_name, max_force) in data_list.items():
                    max_forces_group[max_force][(t1, t2)] = tab_name

                bn = 0
                for max_force, data in max_forces_group.items():
                    interDynamicTab = espressopp.interaction.VerletListDynamicResolutionTabulated(vl, False)
                    for (t1, t2), tab_name in data.items():
                        espp_tab_name = '{}.pot'.format(tab_name.replace('.xvg', ''))
                        if not os.path.exists(espp_tab_name):
                            print('Convert {} to {}'.format(tab_name, espp_tab_name))
                            espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
                        print('Set dynamic resolution potential {}-{} (max force: {})'.format(
                            t1, t2, max_force))
                        interDynamicTab.setPotential(
                            type1=t1,
                            type2=t2,
                            potential=espressopp.interaction.Tabulated(1, espp_tab_name, cutoff=tab_cutoff))
                    if max_force != -1:
                        interDynamicTab.setMaxForce(max_force)
                    system.addInteraction(interDynamicTab, 'tab-dynamic_{}'.format(bn))
                    bn += 1
            else:
                raise RuntimeError('Currently {} not supported'.format(func))

    if has_lj_interaction:
        print('Adding lj interaction')
        system.addInteraction(lj_interaction, 'lj')
    if has_tab_interaction:
        print('Adding lj-tab interaction')
        system.addInteraction(tab_interaction, 'lj-tab')

    if has_tab_capped_interaction:
        print('Adding lj-tab_cap interaction')
        system.addInteraction(tab_capped_interaction, 'lj-tab_cap')

    return cr_observs, particle_pair_scales


def set_bonded_interactions(system, gt, dynamic_type_ids, change_bond_types=set(), name='bonds'):
    """Set bonded interactions.

        Args:
            system: The espressopp.System object.
            gt: The GromacsTopology object.
            dynamic_type_ids: The set of particle types that can change during the simulation.
            change_bond_types: The set of bond types that can be updated during the simulation.
            name: The prefix for the interaction. (default: 'bonds')

        Return:
            tuple with dictionary of dynamic bond fixed pair list and list of static fixed pair list.
    """
    def convert_params(func, raw_data):
        if func == 1:
            return {'K': float(raw_data[1])/2.0, 'r0': float(raw_data[0])}
        elif func == 8:
            espp_tab_name = 'table_b{}.pot'.format(int(raw_data[0]))
            tab_name = 'table_b{}.xvg'.format(int(raw_data[0]))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(tab_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
            return {'itype': 1, 'filename': espp_tab_name}
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

    dynamics_bonds_by_func = collections.defaultdict(list)

    bondparams_func = collections.defaultdict(list)
    dynamic_ptypes = {}
    for pt, p in gt.bondparams.items():
        # Bond types are not in the list of dynamic bond types or in set of change
        if set(pt) - dynamic_type_ids == set(pt) and tuple(pt) not in change_bond_types:
            continue
        params = p['params']
        if tuple(params) not in set(dynamic_ptypes.values()):
            dynamic_ptypes[tuple(sorted(pt))] = tuple(params)
            bondparams_func[p['func']].append((pt, p))
            if p['func'] not in dynamics_bonds_by_func:
                dynamics_bonds_by_func[p['func']] = []
        else:
            dynamic_ptypes = {k: v for k, v in dynamic_ptypes.items() if v != tuple(params)}
            bondparams_func[p['func']] = [x for x in bondparams_func[p['func']][:] if x != (pt, p)]

    # Sort existing bond lists by functional type and select if it is dynamic bond or static.
    bonds_by_func = collections.defaultdict(dict)
    for b, parameters in gt.bonds.items():
        ptypes = tuple(sorted(map(lambda x: gt.atoms[x]['type_id'], b)))
        is_dynamic_bond = ptypes in dynamic_ptypes
        if parameters:  # Has parameters on the list.
            func = int(parameters[0])
            params = tuple(map(float, parameters[1:]))
        else:  # Without parameters, take from bondtypes
            params = gt.bondparams[ptypes]
            if not params:
                params = gt.bonds[tuple(reversed(ptypes))]
            func = int(params['func'])
            params = tuple(map(float, params['params']))
        if is_dynamic_bond:
            dynamics_bonds_by_func[func].append(b)
        else:
            if params not in bonds_by_func[func]:
                bonds_by_func[func][params] = []
            bonds_by_func[func][params].append(b)

    bond_count = 0
    # Process the dynamic bonds, those are stored in the FixedPairTypesList where
    # potential is choose depends on the type of particles.
    dynamics_fpls = collections.defaultdict(dict)
    static_fpls = []
    for func, b_list in dynamics_bonds_by_func.items():
        fpl = espressopp.FixedPairList(system.storage)
        fpl_params = collections.defaultdict(dict)
        fpl.params = fpl_params
        fpl.addBonds(b_list)
        interaction_class, potential_class = func2interaction_dynamic.get(func)
        interaction = interaction_class(system, fpl)
        observe_list = False
        for t, params in bondparams_func[func]:
            observe_list = observe_list or (t in change_bond_types)
            interaction.setPotential(
                type1=t[0], type2=t[1],
                potential=potential_class(**convert_params(func, params['params'])))
            fpl_params[t[0]][t[1]] = params
            fpl_params[t[1]][t[0]] = params
        system.addInteraction(interaction, 'dyn_{}_{}'.format(name, bond_count))
        bond_count += 1
        dynamics_fpls[collections.namedtuple('dfpls', ['func', 'is_observe_list'])(func, observe_list)] = fpl

    # Set first static bonds, those one that has explicitly the parameters
    for func in bonds_by_func:
        interaction_class, potential_class = func2interaction_static.get(func)
        for params, b_list in bonds_by_func[func].items():
            if b_list:
                fpl = espressopp.FixedPairList(system.storage)
                fpl.addBonds(b_list)
                fpl.params = (func, params)
                static_fpls.append(fpl)
                interaction = interaction_class(system, fpl, potential_class(**convert_params(func, params)))
                system.addInteraction(interaction, '{}_{}'.format(name, bond_count))
                bond_count += 1

    print('Set up bond interactions')
    return dynamics_fpls, static_fpls


def set_angle_interactions(system, gt, dynamic_type_ids, change_angle_types=set(), name='angles'):
    """Set angle interactions."""
    def convert_params(func, raw_data):
        if func == 1:
            return {'K': float(raw_data[1])/2.0, 'theta0': float(raw_data[0])*2*math.pi/360}
        elif func == 8:
            espp_tab_name = 'table_a{}.pot'.format(int(raw_data[0]))
            tab_name = 'table_a{}.xvg'.format(int(raw_data[0]))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(tab_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
            return {'itype': 1, 'filename': espp_tab_name}
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

    dynamics_angles_by_func = collections.defaultdict(list)

    angleparams_func = collections.defaultdict(list)
    dynamic_ptypes = {}
    for pt, p in gt.angleparams.items():
        # Bond types are not in the list of dynamic angle types or in set of change
        if (set(pt) - dynamic_type_ids == set(pt) and ((pt[0], pt[1]) not in change_angle_types
            or (pt[1], pt[2]) not in change_angle_types
            or (pt[1], pt[0]) not in change_angle_types
            or (pt[2], pt[1]) not in change_angle_types)):
            continue
        params = p['params']
        if tuple(params) not in set(dynamic_ptypes.values()):
            dynamic_ptypes[tuple((pt))] = tuple(params)
            angleparams_func[p['func']].append((pt, p))
            if p['func'] not in dynamics_angles_by_func:
                dynamics_angles_by_func[p['func']] = []
        else:
            dynamic_ptypes = {k: v for k, v in dynamic_ptypes.items() if v != tuple(params)}
            angleparams_func[p['func']] = [x for x in angleparams_func[p['func']][:] if x != (pt, p)]

    # Sort existing angle lists by functional type and select if it is dynamic angle or static.
    angles_by_func = collections.defaultdict(dict)
    for b, parameters in gt.angles.items():
        ptypes = tuple(map(lambda x: gt.atoms[x]['type_id'], b))
        is_dynamic_angle = ptypes in dynamic_ptypes
        if parameters:  # Has parameters on the list.
            func = int(parameters[0])
            params = tuple(map(float, parameters[1:]))
        else:  # Without parameters, take from angletypes
            params = gt.angleparams[ptypes]
            if not params:
                params = gt.angles[tuple(reversed(ptypes))]
            func = int(params['func'])
            params = tuple(map(float, params['params']))
        if is_dynamic_angle:
            dynamics_angles_by_func[func].append(b)
        else:
            if params not in angles_by_func[func]:
                angles_by_func[func][params] = []
            angles_by_func[func][params].append(b)

    angle_count = 0
    # Process the dynamic angles, those are stored in the FixedTripleTypesList where
    # potential is choose depends on the type of particles.
    dynamics_ftls = collections.defaultdict(dict)
    static_ftls = []
    for func, b_list in dynamics_angles_by_func.items():
        ftl = espressopp.FixedTripleList(system.storage)
        ftl_params = lambda: collections.defaultdict(ftl_params)
        ftl.params = ftl_params()
        ftl.addTriples(b_list)
        interaction_class, potential_class = func2interaction_dynamic.get(func)
        interaction = interaction_class(system, ftl)
        observe_list = False
        for t, params in angleparams_func[func]:
            observe_list = observe_list or (
                (t[0], t[1]) in change_angle_types or (t[1], t[2]) in change_angle_types
                or (t[1], t[0]) in change_angle_types or (t[2], t[1]) in change_angle_types)
            interaction.setPotential(
                type1=t[0], type2=t[1], type3=t[2],
                potential=potential_class(**convert_params(func, params['params'])))
            ftl.params[t[0]][t[1]][t[2]] = params
            ftl.params[t[2]][t[1]][t[0]] = params
        system.addInteraction(interaction, 'dyn_{}_{}'.format(name, angle_count))
        angle_count += 1
        dynamics_ftls[collections.namedtuple('dftls', ['func', 'is_observe_list'])(func, observe_list)] = ftl

    # Set first static angles, those one that has explicitly the parameters
    for func in angles_by_func:
        interaction_class, potential_class = func2interaction_static.get(func)
        for params, b_list in angles_by_func[func].items():
            if b_list:
                ftl = espressopp.FixedTripleList(system.storage)
                ftl.addTriples(b_list)
                ftl.params = (func, params)
                static_ftls.append(ftl)
                interaction = interaction_class(system, ftl, potential_class(**convert_params(func, params)))
                system.addInteraction(interaction, '{}_{}'.format(name, angle_count))
                angle_count += 1

    print('Set up angle interactions')
    return dynamics_ftls, static_ftls

def set_dihedral_interactions(system, gt, dynamic_type_ids, change_dihedral_types=set(), name='dihedrals'):
    """Set dihedral interactions."""
    def convert_params(func, raw_data):
        if func == 1:
            return {'K': float(raw_data[1]),
                    'phi0': float(raw_data[0])*2*math.pi/360,
                    'multiplicity': int(raw_data[2])}
        elif func == 3:
            t = map(float, raw_data[1:])
            return {'K{}'.format(i): v for i, v in enumerate(t[1:])}
        elif func == 8:
            espp_tab_name = 'table_d{}.pot'.format(int(raw_data[0]))
            tab_name = 'table_d{}.xvg'.format(int(raw_data[0]))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(tab_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
            return {'itype': 1, 'filename': espp_tab_name}
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

    dynamics_dihedrals_by_func = collections.defaultdict(list)

    dihedralparams_func = collections.defaultdict(list)
    dynamic_ptypes = {}
    for pt, p in gt.dihedralparams.items():
        # Bond types are not in the list of dynamic dihedral types or in set of change
        if (set(pt) - dynamic_type_ids == set(pt) and ((pt[0], pt[1]) not in change_dihedral_types
            or (pt[1], pt[2]) not in change_dihedral_types
            or (pt[2], pt[3]) not in change_dihedral_types
            or (pt[1], pt[0]) not in change_dihedral_types
            or (pt[2], pt[1]) not in change_dihedral_types
            or (pt[3], pt[2]) not in change_dihedral_types)):
            continue
        params = p['params']
        if tuple(params) not in set(dynamic_ptypes.values()):
            dynamic_ptypes[tuple((pt))] = tuple(params)
            dihedralparams_func[p['func']].append((pt, p))
            if p['func'] not in dynamics_dihedrals_by_func:
                dynamics_dihedrals_by_func[p['func']] = []
        else:
            dynamic_ptypes = {k: v for k, v in dynamic_ptypes.items() if v != tuple(params)}
            dihedralparams_func[p['func']] = [x for x in dihedralparams_func[p['func']][:] if x != (pt, p)]

    # Sort existing dihedral lists by functional type and select if it is dynamic dihedral or static.
    dihedrals_by_func = collections.defaultdict(dict)
    for b, parameters in gt.dihedrals.items():
        ptypes = tuple(map(lambda x: gt.atoms[x]['type_id'], b))
        is_dynamic_dihedral = ptypes in dynamic_ptypes
        if parameters:  # Has parameters on the list.
            func = int(parameters[0])
            params = tuple(map(float, parameters[1:]))
        else:  # Without parameters, take from dihedraltypes
            params = gt.dihedralparams[ptypes]
            if not params:
                params = gt.dihedrals[tuple(reversed(ptypes))]
            func = int(params['func'])
            params = tuple(map(float, params['params']))
        if is_dynamic_dihedral:
            dynamics_dihedrals_by_func[func].append(b)
        else:
            if params not in dihedrals_by_func[func]:
                dihedrals_by_func[func][params] = []
            dihedrals_by_func[func][params].append(b)

    dihedral_count = 0
    # Process the dynamic dihedrals, those are stored in the FixedTripleTypesList where
    # potential is choose depends on the type of particles.
    dynamics_fqls = collections.defaultdict(dict)
    static_fqls = []
    for func, b_list in dynamics_dihedrals_by_func.items():
        fql = espressopp.FixedTripleList(system.storage)
        fql_params = collections.defaultdict(dict)
        fql.params = fql_params
        fql.addBonds(b_list)
        interaction_class, potential_class = func2interaction_dynamic.get(func)
        interaction = interaction_class(system, fql)
        observe_list = False
        for t, params in dihedralparams_func[func]:
            observe_list = observe_list or (
                (t[0], t[1]) in change_dihedral_types
                or (t[1], t[2]) in change_dihedral_types
                or (t[2], t[3]) in change_dihedral_types
                or (t[1], t[0]) in change_dihedral_types
                or (t[2], t[1]) in change_dihedral_types
                or (t[3], t[2]) in change_dihedral_types)
            interaction.setPotential(
                type1=t[0], type2=t[1], type3=t[2],
                potential=potential_class(**convert_params(func, params['params'])))
            fql.params[t[0]][t[1]][t[2]][t[3]] = params
            fql.params[t[3]][t[2]][t[1]][t[0]] = params
        system.addInteraction(interaction, 'dyn_{}_{}'.format(name, dihedral_count))
        dihedral_count += 1
        dynamics_fqls[collections.namedtuple('dfqls', ['func', 'is_observe_list'])(func, observe_list)] = fql

    # Set first static dihedrals, those one that has explicitly the parameters
    for func in dihedrals_by_func:
        interaction_class, potential_class = func2interaction_static.get(func)
        for params, b_list in dihedrals_by_func[func].items():
            if b_list:
                fql = espressopp.FixedTripleList(system.storage)
                fql.addTriples(b_list)
                fql.params = (func, params)
                static_fqls.append(fql)
                interaction = interaction_class(system, fql, potential_class(**convert_params(func, params)))
                system.addInteraction(interaction, '{}_{}'.format(name, dihedral_count))
                dihedral_count += 1

    print('Set up dihedral interactions')
    return dynamics_fqls, static_fqls


def set_pair_interactions(system, gt, args, dynamic_type_ids):
    fudgeLJ = gt.gt.defaults.get('fudgeLJ', 1.0)
    fudgeQQ = gt.gt.defaults.get('fudgeQQ', 1.0)

    # Sort existing pair lists by functional type.
    static_pairs = {}
    dynamics_pairs = []
    atomparams = gt.gt.atomtypes
    combinationrule = int(gt.gt.defaults['combinationrule'])
    for b, parameters in gt.pairs.items():
        ptypes = map(lambda x: gt.atoms[x]['type_id'], b)
        s_ptypes = set(ptypes)
        is_dynamic_bond = s_ptypes - dynamic_type_ids != s_ptypes
        if is_dynamic_bond:
            dynamics_pairs.append(b)
        else:
            params = tuple(map(float, parameters[1:]))
            if params:
                if params not in static_pairs:
                    static_pairs[params] = []
            else:
                params = gt.pairparams.get(tuple(ptypes))
                if not params:
                    params = gt.pairparams.get(tuple(reversed(ptypes)))
                    if not params and gt.gt.defaults['gen-pairs']:
                        sig_1, eps_1 = atomparams[ptypes[0]]['sigma'], atomparams[ptypes[0]]['epsilon']
                        sig_2, eps_2 = atomparams[ptypes[0]]['sigma'], atomparams[ptypes[0]]['epsilon']
                        sig, eps = combination(sig_1, eps_1, sig_2, eps_2, combinationrule)
                        params = [sig, fudgeLJ*eps]
            static_pairs[params].append(b)

    # Set first static pairs, those one that has explicitly the parameters
    static_fpls = []
    pair_count = 0
    for params, b_list in static_pairs.items():
        if b_list:
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
            if type_1 in dynamic_type_ids or type_2 in dynamic_type_ids:
                type_pairs.add(tuple(sorted([type_1, type_2])))
    atomsym_atomtype = gt.used_atomsym_atomtype
    if type_pairs:
        print('Set up 1-4 pair interactions')
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
        system.addInteraction(interaction, 'dyn_lj14_{}'.format(pair_count))

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
    props = ['id', 'type', 'pos', 'mass', 'q', 'res_id', 'state', 'lambda_adr']
    particle_list = []

    for atom_id in sorted(coordinate.atoms):
        data = coordinate.atoms[atom_id]
        top_data = topol.atoms[atom_id]
        particle_list.append(
            [atom_id,
             top_data['type_id'],
             espressopp.Real3D(data.position),
             top_data['mass'],
             top_data['charge'],
             data.chain_idx,
             top_data.get('state', 0),
             1.0]
        )

    return props, particle_list
