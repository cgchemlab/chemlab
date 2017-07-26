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
import logging
import os
import re

import numpy

__doc__ = "Set of I/O classes and functions."""

logger = logging.getLogger(__name__)

Atom = collections.namedtuple(
    'Atom', [
        'atom_id',
        'name',
        'chain_name',
        'chain_idx',
        'position',
        'velocity'
    ])


class TopoAtom(object):
    """Atom object used in TopologyFiles."""
    atom_id = None
    atom_type = None
    chain_idx = None
    chain_name = None
    name = None
    cgnr = None
    charge = None
    mass = None
    active_site = None

    def __init__(self, atom_id=None, atom_type=None, chain_idx=None, chain_name=None, name=None,
                 cgnr=None, charge=None, mass=None, active_site=None):
        self.atom_id = atom_id
        self.atom_type = atom_type
        self.chain_idx = chain_idx
        self.chain_name = chain_name
        self.name = name
        self.cgnr = cgnr
        self.charge = charge
        self.mass = mass
        self.active_site = active_site

    def __repr__(self):
        return '{} ({}): {} ({}) q={}, m={}, as={}'.format(
            self.atom_id, self.chain_idx, self.name, self.chain_name, self.charge,
            self.mass, self.active_site)


def prepare_path(file_path):
    """Prepare the file to open.

    Args:
      file_path: The file path.

    Returns:
      The path to the file.
    """

    if os.path.exists(file_path):
        file_name = os.path.basename(file_path)
        dir_name = os.path.dirname(file_path)
        if not dir_name:
            dir_name = '.'
        existing_copies = [x for x in os.listdir(dir_name) if x.startswith('_%s' % file_name)]
        if existing_copies:
            max_copy_id = max([int(x.strip('_').split('.')[-1]) for x in existing_copies])
        else:
            max_copy_id = 0
        new_file_name = '_%s.%d_' % (file_name, max_copy_id+1)
        new_file_path = os.path.join(dir_name, new_file_name)
        logger.warning('Found: %s, backup on: %s', file_path, new_file_path)
        os.rename(file_path, new_file_path)

    return file_path


class CoordinateFile(object):
    """Coordinate file object."""

    def __init__(self, file_name):
        self.file_name = file_name
        self.title = None
        self.atoms_updated = False
        self.atoms = {}
        self.fragments = collections.defaultdict(dict)
        self.content = None
        self.box = None
        self.data = None
        self.scale_factor = 1.0
        self.file = None
        self.atoms = {}
        self.chains = {}

    def init(self):
        self.__init__(self.file_name)
        logger.info('Init of coordinate file')
        self.content = None
        self.box = None
        self.file = None
        self.atoms_updated = False


class TopologyFile(object):
    """Reader for GROMACS .top files.

    Args:
        file_name: The input topology file.
    """

    def __init__(self, file_name):
        self.file_name = file_name
        self.title = None
        self.atoms_updated = False
        self.new_data = {
            'bonds': {},
            'angles': {},
            'dihedrals': {},
            'improper_dihedrals': {},
            'pairs': {}
            }
        self.parsers = {}
        self.writers = {}

        self.content = None
        self.file = None

    def init(self):
        self.__init__(self.file_name)
        logger.info('Init of topology file.')
        self.chains = {}
        self.content = None
        self.file = None
        self.atoms_updated = False


class GROFile(CoordinateFile):
    scale_factor = 1.0

    def read(self):
        """Reads the .gro file and return the atom list.

        Returns:
          The dict with atoms (key: atom_id, value: atom object).
        """

        self.file = open(self.file_name, 'r')
        if not self.content:
            self.content = self.file.readlines()

        self.title = self.content[0].replace('\r\n', '').replace('\n', '')
        number_of_atoms = int(self.content[1])

        logger.info('Reading GRO file %s', self.file_name)

        for line in self.content[2:number_of_atoms + 2]:
            chain_idx = int(line[0:5].strip())
            chain_name = line[5:10].strip()
            at_name = line[10:15].strip()
            at_id = int(line[15:20].strip())
            # Nedd to rescale.
            pos_x = float(line[20:28].strip()) * self.scale_factor
            pos_y = float(line[28:36].strip()) * self.scale_factor
            pos_z = float(line[36:44].strip()) * self.scale_factor

            vel = (None, None, None)
            # Velocity in the file
            if len(line) > 45:
                v_x = float(line[44:52].strip()) * self.scale_factor
                v_y = float(line[52:60].strip()) * self.scale_factor
                v_z = float(line[60:68].strip()) * self.scale_factor
                vel = numpy.array([v_x, v_y, v_z])

            self.atoms[at_id] = (
                Atom(
                    atom_id=at_id,
                    name=at_name,
                    chain_name=chain_name,
                    chain_idx=chain_idx,
                    position=numpy.array([pos_x, pos_y, pos_z]),
                    velocity=vel
                ))
            self.fragments[chain_name][at_name] = self.atoms[at_id]
            if chain_name not in self.chains:
                self.chains[chain_name] = {}
            if chain_idx not in self.chains[chain_name]:
                self.chains[chain_name][chain_idx] = {}
            self.chains[chain_name][chain_idx][at_name] = self.atoms[at_id]

        # Reads the box size, the last line.
        self.box = numpy.array(
            map(float, filter(None, self.content[number_of_atoms + 2].split(' ')))
            ) * self.scale_factor

    def write(self, file_name=None, force=False):
        """Writes the content to the output file.

        Args:
          file_name: The new file name, otherwise the old one will be used.
          force: Force to save even if any atoms were not updated.
        """

        if self.atoms_updated or force:
            output = []
            if self.title:
                output.append(self.title)
            else:
                output.append('XXX of molecules')
            # Puts the number of atoms
            output.append('%d' % len(self.atoms))
            # Puts the definition of the atoms, fixed format.
            fmt = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
            fmt_vel = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f"
            for at_id in sorted(self.atoms):
                at = self.atoms[at_id]
                if at.velocity[0] is None:
                    output.append(fmt % (
                        at.chain_idx,
                        at.chain_name,
                        at.name,
                        at.atom_id,
                        at.position[0], at.position[1], at.position[2]))
                else:
                    output.append(fmt_vel % (
                        at.chain_idx,
                        at.chain_name,
                        at.name,
                        at.atom_id,
                        at.position[0], at.position[1], at.position[2],
                        at.velocity[0], at.velocity[1], at.velocity[2]))

            output.append('%f %f %f\n' % tuple(self.box))
            write_file_path = prepare_path(file_name if file_name else self.file_name)
            logger.info('Writing GRO file %s', write_file_path)
            output_file = open(write_file_path, 'w')
            output_file.writelines('\n'.join(output))
            output_file.close()
            self.atoms_updated = False

    def update_position(self, system, unfolded=False):
        """Updates position based on curent state

        Args:
            system: The espressopp.System object.
            unfolded: If set to True then write in unfoded state.
        """
        if unfolded:
            boxL = system.bc.boxL
            for pid in self.atoms:
                p = system.storage.getParticle(pid)
                self.atoms[pid] = self.atoms[pid]._replace(
                    position=[p.pos[x] + p.imageBox[x]*boxL[x] for x in range(3)],
                    velocity=[p.v[x] for x in range(3)]
                )
        else:
            for pid in self.atoms:
                p = system.storage.getParticle(pid)
                self.atoms[pid] = self.atoms[pid]._replace(position=p.pos)

    @classmethod
    def load_data(cls, system, file_name, name_seq, chain_name_seq, particle_ids):
        """Load position of particles and creates new GROFile object

        Args:
            system: The espressopp.System object.
            file_name: The name of file.
            name_seq: The sequence of atom names.
            chain_name_seq: The sequence of chain names.
            particle_ids: The list of particle ids to include in the file

        Returns:
            GROFile object.
        """
        f = cls(file_name)
        idx = 0
        for pid in particle_ids:
            p = system.storage.getParticle(pid)
            f.atoms[pid] = Atom(
                atom_id=pid,
                name=name_seq[idx % len(name_seq)],
                chain_name=chain_name_seq[idx % len(chain_name_seq)],
                chain_idx=1 + idx / len(chain_name_seq),
                position=p.pos
            )
            idx += 1
        f.box = system.bc.boxL
        return f


class PDBFile(CoordinateFile):
    scale_factor = 0.1  # PDB is expressed in Angstrome and the program use nm

    def read(self):
        """Reads the file and return atom list."""

        self.file = open(self.file_name, 'r')

        if not self.content:
            self.content = self.file.readlines()

        logger.info('Reading PDB file %s', self.file_name)

        for line in self.content:
            if line.startswith('CRYST1'):
                # Box size
                self.box = numpy.array(
                    map(float, filter(None, line.split(' '))[1:4])
                    ) * self.scale_factor
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                atom_id = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                chain_name = line[17:20].strip()  # Residue name
                chain_idx = line[22:26].strip()
                pos_x = float(line[30:38]) * self.scale_factor
                pos_y = float(line[38:46]) * self.scale_factor
                pos_z = float(line[46:54]) * self.scale_factor
                self.atoms[atom_id] = (
                    Atom(
                        atom_id=atom_id,
                        name=atom_name,
                        chain_name=chain_name,
                        chain_idx=chain_idx,
                        position=numpy.array([pos_x, pos_y, pos_z])
                        ))
                self.fragments[chain_name][atom_name] = self.atoms[atom_id]

        if len([x for x in self.box if x == self.box[0]]) != 3:
            raise ValueError('The box size in all direction should be the same')

    def write(self, file_name=None, force=False):
        """Write the file again."""
        if self.atoms_updated or force:
            output = []
            output.append('REMARK generate by YAPT')
            output.append('MODEL 1')
            # Writing the box coordinates
            # Following http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM
            # Boxes are orthorhombic for now
            output.append('%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n' % (
                'CRYST1',
                self.box[0] / self.scale_factor,
                self.box[1] / self.scale_factor,
                self.box[2] / self.scale_factor,
                90.00,
                90.00,
                90,
                'P 1',
                1
                ))

            # Puts the number of atoms
            output.append('%d' % len(self.atoms))
            # Puts the definition of the atoms, fixed format.
            fmt = '%-6s%5d %4s %-3s  %4d    %8.3f%8.3f%8.3f                      %2s'
            for at_id in sorted(self.atoms):
                at = self.atoms[at_id]
                output.append(fmt % (
                    'ATOM  ',
                    int(at.atom_id) % 100000,
                    at.name,
                    at.chain_name,
                    int(at.chain_idx) % 10000,
                    at.position[0] / self.scale_factor,
                    at.position[1] / self.scale_factor,
                    at.position[2] / self.scale_factor,
                    at.name
                    ))

            output.append('TER')
            output.append('ENDMDL')

            write_file_path = prepare_path(file_name if file_name else self.file_name)
            logger.info('Writing PDB file %s', write_file_path)
            output_file = open(write_file_path, 'w')
            output_file.writelines('\n'.join(output))
            output_file.close()
            self.atoms_updated = False


class GROMACSTopologyFile(TopologyFile):
    """Very basic representation of topology file."""

    def __init__(self, file_name):
        super(GROMACSTopologyFile, self).__init__(file_name)

        self.parsers = {
            'defaults': self._parse_defaults,
            'atomtypes': self._parse_atomtypes,
            'atomstate': self._parse_atomstate,
            'nonbond_params': self._parse_nonbond_params,
            'bondtypes': self._parse_bondtypes,
            'angletypes': self._parse_angletypes,
            'dihedraltypes': self._parse_dihedraltypes,
            'atoms': self._parse_atoms,
            'bonds': self._parse_bonds,
            'dihedrals': self._parse_dihedrals,
            'improper_dihedrals': self._parse_improper_dihedrals,
            'pairs': self._parse_pairs,
            'angles': self._parse_angles,
            'moleculetype': self._parse_moleculetype,
            'system': self._parse_system,
            'molecules': self._parse_molecules
            }

        self.writers = {
            'defaults': self._write_defaults,
            'moleculetype': self._write_moleculetype,
            'system': self._write_system,
            'molecules': self._write_molecules,
            'atomtypes': self._write_atomtypes,
            'atomstate': self._write_atomstate,
            'bondtypes': self._write_bondtypes,
            'angletypes': self._write_angletypes,
            'dihedraltypes': self._write_dihedraltypes,
            'nonbond_params': self._write_nonbondparams,
            'atoms': self._write_atoms,
            'bonds': self._write_bonds,
            'angles': self._write_angles,
            'dihedrals': self._write_dihedrals,
            'pairs': self._write_pairs,
            'cross_bonds': lambda: self._write_default(
                [self.new_data.get('cross_bonds'), self.cross_bonds]),
            'cross_angles': lambda: self._write_default(
                [self.new_data.get('cross_angles'), self.cross_angles]),
            'cross_dihedrals': lambda: self._write_default(
                [self.new_data.get('cross_dihedrals'), self.cross_dihedrals]),
            'cross_pairs': lambda: self._write_default(
                [self.new_data.get('cross_pairs'), self.cross_pairs])
            }
        self.current_charges = {}
        self.atomtypes = {}
        self.atomnr2atom_name = collections.defaultdict(list)
        self.atom_name2atomnr = {}
        self.atomstate = {}
        self.nonbond_params = {}
        self.bondtypes = {}
        self.angletypes = {}
        self.dihedraltypes = {}
        self.header_section = []
        self.defaults = None
        self.moleculetype = {}
        self.molecules = []
        self.system_name = None

        self.atoms = {}
        # Store bond, angle, dihedral, pairs list by the moleculetype.
        self.current_molecule = None
        self.molecules_data = collections.defaultdict(dict)

    def init(self):
        """Reset the class properties without creating the object again."""
        logger.info('Init of topology file.')
        self.new_data = {
            'bonds': {},
            'angles': {},
            'dihedrals': {},
            'improper_dihedrals': {},
            'pairs': {}
            }
        if '__state' in self.__dict__:
            del self.__dict__['__state']
        self.current_charges = {}
        self.atoms_updated = False

    def update_position(self, pdbfile):
        """Reads the position data from the coordinate file and update the atoms.

        Args:
          pdbfile: The pdb file.
        """
        raise NotImplemented('method not implemented')

        # logger.info('Update position from file %s', pdbfile.file_name)
        # for k, v in pdbfile.atoms.iteritems():
        #     self.atoms[k].position = v.position

    def read(self):
        """Reads the topology file."""

        if not self.content:
            self.file = open(self.file_name, 'r')
            self.content = self.file.readlines()

        logger.info('Reading top file %s', self.file_name)

        # New version
        current_parser = None
        visited_sections = set()
        section_name = None
        previous_section = None
        for raw_line in self.content:
            line = re.sub(';.*$', '', raw_line.strip())
            if line.startswith(';') or line.startswith('#') or len(line) == 0:
                continue
            elif line.startswith('['):  # Section
                previous_section = section_name
                section_name = line.replace('[', '').replace(']', '').strip()
                # Hack for GROMACS improper_dihedrals
                if previous_section == 'dihedrals' and section_name == 'dihedrals':
                    section_name = 'improper_dihedrals'
                current_parser = self.parsers.get(section_name)
                if current_parser is not None:
                    print('{}: Reading section {}'.format(self.file_name, section_name))
                else:
                    print('Parser for section {} not defined'.format(section_name))
                visited_sections.add(previous_section)
            else:
                if current_parser is not None: # and section_name not in visited_sections:
                    raw_data = filter(None, line.split())
                    if raw_data:
                        current_parser(raw_data)  # pylint:disable=E1102

    def write(self, filename=None):  #NOQA
        """Updates the topology file.

        Args:
          filename: The optional output filename.
        """
        if filename is None:
            filename = self.file_name
        output_file = open(prepare_path(filename), 'w')

        new_data = []
        current_section = None
        skip_lines = False
        section_writer = None

        if self.content is None:
            sections = []
            if self.defaults:
                sections.append('defaults')
            if self.atomtypes:
                sections.append('atomtypes')
            if self.bondtypes:
                sections.append('bondtypes')
            if self.angletypes:
                sections.append('angletypes')
            if self.dihedraltypes:
                sections.append('dihedraltypes')
            if self.nonbond_params:
                sections.append('nonbond_params')
            if self.atomstate:
                sections.append('atomstate')

            sections.extend([
                'moleculetype',
                'atoms',
                'bonds',
                'angles',
                'dihedrals',
                'pairs',
                'system',
                'molecules'])
            self.content = []
            for s in sections:
                self.content.append('[ %s ]\n' % s)
                self.content.append('\n')

        new_data.extend(self.header_section)

        for line in self.content:
            tmp_line = line.strip()
            if tmp_line.startswith('['):  # section part
                new_data.append('\n')  # space before section name
                new_data.append(line)
                previous_section = current_section
                current_section = tmp_line.replace('[', '').replace(']', '').strip()
                if previous_section == 'dihedrals' and current_section == 'dihedrals':
                    current_section = 'improper_dihedrals'
                section_writer = self.writers.get(current_section)
                skip_lines = False
            elif tmp_line.startswith(';') or tmp_line.startswith('#'):
                new_data.append(line)
            else:
                if section_writer is None:  # there is no special writer, simply copy the line
                    new_data.append(line)
                elif not skip_lines:
                    output_writer = section_writer()
                    if output_writer:
                        new_data.extend(['%s\n' % x for x in output_writer])
                    new_data.extend(['\n'])
                    skip_lines = True
        # Add missing new lines if not present
        new_data = ['{}\n'.format(x) if not x.endswith('\n') else x for x in new_data]

        logger.debug('Writing topology file %s...', filename)
        output_file.writelines(new_data)
        output_file.close()
        self.atoms_updated = False

    # Parsers for the data.
    def _parse_defaults(self, raw_data):
        self.defaults = {
            'func': int(raw_data[0]),
            'combinationrule': int(raw_data[1])
            }
        if len(raw_data) > 2:
            self.defaults['gen-pairs'] = raw_data[2] == 'yes'
            self.defaults['fudgeLJ'] = float(raw_data[3])
            self.defaults['fudgeQQ'] = float(raw_data[4])
        else:
            self.defaults['gen-pairs'] = False
            self.defaults['fudgeLJ'] = 1.0
            self.defaults['fudgeQQ'] = 1.0
        self.defaults['nbfunc'] = 1

    def _parse_atomtypes(self, raw_data):
        has_data = True
        if len(raw_data) == 7:
            atom_name = raw_data[0]
            atom_nr = raw_data[0]
            atom_mass = float(raw_data[2])
            atom_q = float(raw_data[3])
            atom_type = raw_data[4]
            sigma = float(raw_data[5])
            epsilon = float(raw_data[6])
        elif len(raw_data) == 6:
            atom_name = raw_data[0]
            atom_nr = raw_data[0]
            atom_mass = float(raw_data[1])
            atom_q = float(raw_data[2])
            atom_type = raw_data[3]
            sigma = float(raw_data[4])
            epsilon = float(raw_data[5])
        elif len(raw_data) == 8 and raw_data[0].startswith('opls'):
            atom_name = raw_data[0]
            atom_nr = raw_data[1]
            atom_mass = float(raw_data[3])
            atom_q = float(raw_data[4])
            atom_type = raw_data[5]
            # Data directly in sigma, epsilon, no need to convert.
            sigma = float(raw_data[6])
            epsilon = float(raw_data[7])
        else:
            print('Skip atom type {}'.format(raw_data[0]))
            has_data = False

        if has_data:
            self.atom_name2atomnr[atom_name] = atom_nr
            self.atomnr2atom_name[atom_nr].append(atom_name)
            self.atomtypes[atom_name] = {
                'name': atom_name,
                'mass': atom_mass,
                'charge': atom_q,
                'type': atom_type,
                'sigma': sigma,
                'epsilon': epsilon
            }
            #print self.atomtypes[atom_name]
            if atom_name in self.atomstate:
                self.atomtypes[atom_name]['state'] = self.atomstate[atom_name]


    def _parse_nonbond_params(self, raw_data):
        k = tuple(sorted(raw_data[:2]))
        if k in self.nonbond_params:
            raise RuntimeError('{} already exists, wrong topology'.format(k))
        self.nonbond_params[k] = {
            'func': int(raw_data[2]),
            'params': raw_data[3:]}

    def _parse_atomstate(self, raw_data):
        atom_type = raw_data[0]
        atom_state = int(raw_data[1])
        if atom_type in self.atomtypes:
            self.atomtypes[atom_type]['state'] = atom_state
        self.atomstate[atom_type] = atom_state

    def _parse_bondtypes(self, raw_data):
        i, j = raw_data[:2]
        if i not in self.bondtypes:
            self.bondtypes[i] = {}
        if j not in self.bondtypes:
            self.bondtypes[j] = {}

        self.bondtypes[i][j] = {
            'func': int(raw_data[2]),
            'params': raw_data[3:]
        }
        self.bondtypes[j][i] = self.bondtypes[i][j]

    def _parse_angletypes(self, raw_data):
        i, j, k = raw_data[:3]
        if i not in self.angletypes:
            self.angletypes[i] = {}
        if j not in self.angletypes[i]:
            self.angletypes[i][j] = {}
        if k not in self.angletypes:
            self.angletypes[k] = {}
        if j not in self.angletypes[k]:
            self.angletypes[k][j] = {}

        self.angletypes[i][j][k] = {
            'func': int(raw_data[3]),
            'params': raw_data[4:]
        }
        self.angletypes[k][j][i] = self.angletypes[i][j][k]

    def _parse_dihedraltypes(self, raw_data):
        i, j, k, l = raw_data[:4]
        if i not in self.dihedraltypes:
            self.dihedraltypes[i] = {}
        if j not in self.dihedraltypes[i]:
            self.dihedraltypes[i][j] = {}
        if k not in self.dihedraltypes[i][j]:
            self.dihedraltypes[i][j][k] = {}
        if l not in self.dihedraltypes:
            self.dihedraltypes[l] = {}
        if k not in self.dihedraltypes[l]:
            self.dihedraltypes[l][k] = {}
        if j not in self.dihedraltypes[l][k]:
            self.dihedraltypes[l][k][j] = {}

        try:
            self.dihedraltypes[i][j][k][l] = {
                'func': int(raw_data[4]),
                'params': raw_data[5:]
            }
            self.dihedraltypes[l][k][j][i] = self.dihedraltypes[i][j][k][l]
        except ValueError:
            print('Skip {}'.format(raw_data))

    def _parse_atoms(self, raw_data):
        if self.current_molecule is None:
            raise RuntimeError('Wrong order, before bonds there should be a moleculetype section')

        at = TopoAtom()
        at.atom_id = int(raw_data[0])
        at.atom_type = raw_data[1]
        at.chain_idx = int(raw_data[2])
        at.chain_name = raw_data[3]
        at.name = raw_data[4]
        at.cgnr = int(raw_data[5])
        at.molecule_name = self.current_molecule

        if len(raw_data) > 6:
            at.charge = float(raw_data[6])
        if len(raw_data) > 7:
            at.mass = float(raw_data[7])

        if 'atoms' not in self.molecules_data[self.current_molecule]:
            self.molecules_data[self.current_molecule]['atoms'] = {}
        self.molecules_data[self.current_molecule]['atoms'][at.atom_id] = at

    def _parse_bonds(self, raw_data):
        atom_tuple = tuple(map(int, raw_data[0:2]))
        if self.current_molecule is None:
            raise RuntimeError('Wrong order, before bonds there should be a moleculetype section')
        if 'bonds' not in self.molecules_data[self.current_molecule]:
            self.molecules_data[self.current_molecule]['bonds'] = {}
        self.molecules_data[self.current_molecule]['bonds'][atom_tuple] = raw_data[2:]

    def _parse_angles(self, raw_data):
        atom_tuple = tuple(map(int, raw_data[0:3]))

        if self.current_molecule is None:
            raise RuntimeError('Wrong order, before bonds there should be a moleculetype section')
        if 'angles' not in self.molecules_data[self.current_molecule]:
            self.molecules_data[self.current_molecule]['angles'] = {}
        self.molecules_data[self.current_molecule]['angles'][atom_tuple] = raw_data[3:]

    def _parse_dihedrals(self, raw_data):
        atom_tuple = tuple(map(int, raw_data[0:4]))

        if self.current_molecule is None:
            raise RuntimeError('Wrong order, before bonds there should be a moleculetype section')
        if 'dihedrals' not in self.molecules_data[self.current_molecule]:
            self.molecules_data[self.current_molecule]['dihedrals'] = {}
        self.molecules_data[self.current_molecule]['dihedrals'][atom_tuple] = raw_data[4:]

    def _parse_improper_dihedrals(self, raw_data):
        atom_tuple = tuple(map(int, raw_data[0:4]))

        if self.current_molecule is None:
            raise RuntimeError('Wrong order, before bonds there should be a moleculetype section')
        if 'improper_dihedrals' not in self.molecules_data[self.current_molecule]:
            self.molecules_data[self.current_molecule]['improper_dihedrals'] = {}
        self.molecules_data[self.current_molecule]['improper_dihedrals'][atom_tuple] = raw_data[3:]

    def _parse_pairs(self, raw_data):
        atom_tuple = tuple(map(int, raw_data[0:2]))

        if self.current_molecule is None:
            raise RuntimeError('Wrong order, before bonds there should be a moleculetype section')
        if 'pairs' not in self.molecules_data[self.current_molecule]:
            self.molecules_data[self.current_molecule]['pairs'] = {}
        self.molecules_data[self.current_molecule]['pairs'][atom_tuple] = raw_data[2:]


    def _parse_moleculetype(self, raw_data):
        self.current_molecule = raw_data[0]
        self.moleculetype[raw_data[0]] = int(raw_data[1])
        print('Reading molecule {}, nrexcl={}'.format(raw_data[0], raw_data[1]))

    def _parse_molecules(self, raw_data):
        # Orders matter.
        print('{} - {} molecules'.format(raw_data[0], int(raw_data[1])))
        self.molecules.append((raw_data[0], int(raw_data[1])))

    def _parse_system(self, raw_data):
        self.system_name = raw_data[0]

    # Writers
    def _write_atoms(self):
        return_data = []
        atoms = {}
        for mol_data in self.molecules_data.values():
            atoms.update(mol_data['atoms'])

        for atom_id in sorted(atoms):
            x = atoms[atom_id]
            return_data.append('%s %s %s %s %s %s %s %s' % (
                               x.atom_id,
                               x.atom_type,
                               x.chain_idx,
                               x.chain_name,
                               x.name,
                               x.cgnr,
                               x.charge if x.charge is not None else '0.0',
                               x.mass if x.mass is not None else ''
                               ))
        return return_data

    def _write_atomtypes(self):
        return_data = []
        for atom_type, values in self.atomtypes.items():
            return_data.append('{name} {mass} {charge} {type} {sigma} {epsilon}'.format(
                **values))
        return return_data

    def _write_nonbondparams(self):
        return_data = []
        for types, params in self.nonbond_params.items():
            return_data.append('{} {} {} {}'.format(
                types[0], types[1], params['func'], ' '.join(map(str, params['params']))))

        return return_data

    def _write_atomstate(self):
        return_data = []
        for val in self.atomstate.items():
            return_data.append('{} {}'.format(*val))
        return return_data

    def _write_bondtypes(self):
        return_data = []
        for i in self.bondtypes:
            for j, params in self.bondtypes[i].items():
                return_data.append('{} {} {} {}'.format(i, j, params['func'], ' '.join(params['params'])))
        return return_data

    def _write_angletypes(self):
        return_data = []
        for i in self.angletypes:
            for j in self.angletypes[i]:
                for k, params in self.angletypes[i][j].items():
                    return_data.append('{} {} {} {} {}'.format(
                        i, j, k, params['func'], ' '.join(params['params'])))
        return return_data

    def _write_dihedraltypes(self):
        return_data = []
        for i in self.dihedraltypes:
            for j in self.dihedraltypes[i]:
                for k in self.dihedraltypes[i][j]:
                    for l, params in self.dihedraltypes[i][j][k].items():
                        return_data.append('{} {} {} {} {} {}'.format(
                            i, j, k, l, params['func'], ' '.join(params['params'])))
        return return_data

    def _write_bonds(self):  # pylint:disable=R0201
        return_data = []
        total_bonds = {}
        for mol_data in self.molecules_data.values():
            return_data.extend(self._write_default(mol_data.get('bonds', {})))
            total_bonds.update(mol_data.get('bonds', {}))
        return_data.extend(self._write_default(self.new_data['bonds'], total_bonds))
        return return_data

    def _write_pairs(self):  # pylint:disable=R0201
        return_data = []
        total_pairs = {}
        for mol_data in self.molecules_data.values():
            return_data.extend(self._write_default(mol_data.get('pairs', {})))
            total_pairs.update(mol_data.get('pairs', {}))
        return_data.extend(self._write_default(self.new_data['pairs'], total_pairs))
        return return_data

    def _write_angles(self):
        return_data = []
        total_angles = {}
        for mol_data in self.molecules_data.values():
            return_data.extend(self._write_default(mol_data.get('angles', {})))
            total_angles.update(mol_data.get('angles', {}))
        return_data.extend(self._write_default(self.new_data['angles'], total_angles))
        return return_data

    def _write_dihedrals(self):
        return_data = []
        total_dihedrals = {}
        for mol_data in self.molecules_data.values():
            return_data.extend(self._write_default(mol_data.get('dihedrals', {})))
            total_dihedrals.update(mol_data.get('dihedrals', {}))
        return_data.extend(self._write_default(self.new_data['dihedrals'], total_dihedrals))
        return return_data

    def _write_improper_dihedrals(self):
        return_data = []
        total_impropers = {}
        for mol_data in self.molecules_data.values():
            return_data.extend(self._write_default(mol_data.get('improper_dihedrals', {})))
            total_impropers.update(mol_data.get('improper_dihedrals', {}))

        return_data.extend(self._write_default(self.new_data['improper_dihedrals'], total_impropers))
        return return_data

    def _write_defaults(self):
        if self.defaults:
            defaults = self.defaults.copy()
            defaults['gen-pairs'] = 'yes' if defaults['gen-pairs'] else 'no'
            return ['{nbfunc} {combinationrule} {gen-pairs} {fudgeLJ} {fudgeQQ}'.format(**defaults)]
        return []

    def _write_moleculetype(self):
        return_data = []
        for mol_name, nrexcl in self.moleculetype.items():
            return_data.append('{} {}'.format(mol_name, nrexcl))

        return return_data

    def _write_system(self):
        return [self.system_name]

    def _write_molecules(self):
        return ['{} {}'.format(*x) for x in self.molecules]

    def _write_default(self, datas=None, check_in=None):  # pylint:disable=R0201
        if check_in is None:
            check_in = []

        if datas is None:
            return False

        if not isinstance(datas, list):
            datas = [datas]

        if None in datas:
            return False

        flat_data = []
        for data in datas:
            for key, values in data.iteritems():
                rev_key = tuple(reversed(key))
                if tuple(key) not in check_in or rev_key not in check_in or rev_key not in data:
                    flat_data.append(list(key) + list(values))
        flat_data.sort()
        return ['%s' % ' '.join(map(str, x)) for x in flat_data]


class LammpsReader(object):
    """Very simple LAMMPS data file and input parser."""
    def __init__(self):
        self.data_parsers = {
            'Atoms': self._read_atom,
            'Velocities': self._read_velocity,
            'Masses': self._read_mass,
            'Bonds': self._read_bond,
            'Angles': self._read_angle,
            'Dihedrals': self._read_dihedral,
            'Impropers': self._read_improper,
            'header': self._read_header,
            'coeffs': self._read_coeff,
        }
        self.force_field = collections.defaultdict(dict)
        self.init()
        # We shift box to the origin and we need to do the same with atom position.
        self._box_translate = {}

    def init(self):
        self.current_section = 'header'
        self.previous_section = None

        # Data structures
        self._item_counters = {}
        self._type_counters = {}
        self._mass_type = {}
        self._section_line = None

        self.box = {}
        self.atoms = collections.defaultdict(dict)
        self.topology = {
            'bonds': collections.defaultdict(list),
            'angles': collections.defaultdict(list),
            'dihedrals': collections.defaultdict(list),
            'impropers': collections.defaultdict(list),
            }
        self.distance_scale_factor = 0.1

    def read_data(self, file_name, scale_factor=None, update=False):
        """Reads data file written with write_data command.

        Arsgs:
            file_name: The name of data file to read.
            scale_factor: The factor by which every distance quantity will be multiply.
        """
        if update:
            self.init()

        if scale_factor is not None:
            self.distance_scale_factor = scale_factor

        with open(file_name, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                section_line = line.split('#')[0].strip()
                if section_line in self.data_parsers:
                    print('{}: Reading section {}'.format(file_name, line))
                    self.previous_section = self.current_section
                    self.current_section = section_line
                elif 'Coeff' in section_line:
                    print('{}: Reading coefficient section {}'.format(file_name, line))
                    self.previous_section = self.current_section
                    self.current_section = 'coeffs'
                    self._section_line = section_line
                elif self.current_section is not None:
                    self.data_parsers[self.current_section](line)
                else:
                    self.previous_section = self.current_section
                    self.current_section = None

    def read_input(self, file_name):
        """Reads LAMMPS input script. Only take cares on *_style and pair_coeff

        Args:
            file_name: Name of input file.
        """
        with open(file_name, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:  # skip empty lines
                    continue
                if '_style' in line:
                    sp_line = line.split()
                    self.force_field[sp_line[0]] = sp_line[1:]
                elif 'bond_coeff' in line or 'angle_coeff' in line or 'dihedral_coeff' in line:
                    sp_line = line.split()
                    stype = sp_line[0].replace('_coeff', '')
                    btype = sp_line[1].strip()
                    self.force_field[stype][btype] = sp_line[2:]
                elif 'pair_coeff' in line:
                    sp_line = line.split()
                    at_1, at_2 = sp_line[1:3]
                    if '*' not in sp_line[1]:
                        at_1 = int(at_1)
                    if '*' not in sp_line[2]:
                        at_2 = int(at_2)
                    self.force_field['pair_coeff'][tuple(sorted((at_1, at_2)))] = sp_line[3:]
                elif 'units' in line:
                    self.units = line.split()[1]
                    if self.units == 'real':
                        self.distance_scale_factor = 10**-1

    def update_atoms(self, file_name):
        """Reads LAMMPS data file and updates only atom section.

        Args:
            file_name: Input data file with Atoms section.
        """
        with open(file_name, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                section_line = line.split('#')[0].strip()
                if section_line in self.data_parsers:
                    self.previous_section = self.current_section
                    self.current_section = section_line
                    if section_line == 'Atoms':
                        print('Found "Atoms" section, updating atoms')
                elif self.current_section is not None and self.current_section == 'Atoms':
                    self._read_atom(line, update=True)
                else:
                    self.current_section = None

    # Parsers section
    def _read_header(self, input_line):
        """Parses header of data file."""
        sp_line = input_line.split()
        if 'types' in sp_line:
            self._type_counters[sp_line[1]] = int(sp_line[0])
        elif 'xhi' in sp_line or 'yhi' in sp_line or 'zhi' in sp_line:
            lo, hi = map(float, sp_line[:2])
            lo *= self.distance_scale_factor
            hi *= self.distance_scale_factor
            tag = sp_line[-1].replace('hi', '')
            self._box_translate[tag] = lo
            self.box[tag] = hi - lo
        elif ('atoms' in sp_line or 'bonds' in sp_line or 'angles' in sp_line
              or 'dihedrals' in sp_line):
            self._item_counters[sp_line[1]] = int(sp_line[0])

    def _read_coeff(self, input_line):
        coeff_type, _ = self._section_line.split()
        coeff_type = coeff_type.lower()

        sp_line = input_line.split()
        ff_type = int(sp_line[0])
        self.force_field[coeff_type][ff_type] = sp_line[1:]

    def _read_atom(self, input_line, update=False):
        sp_line = input_line.split()
        # Set type
        sp_line[:3] = map(int, sp_line[:3])
        sp_line[3:7] = map(float, sp_line[3:7])
        if len(sp_line) == 10:
            sp_line[7:10] = map(int, sp_line[7:10])
            at_id, at_tag, at_type, q, x, y, z, nx, ny, nz = sp_line
        else:
            at_id, at_tag, at_type, q, x, y, z = sp_line
            nx, ny, nz = None, None, None

        if at_id > self._item_counters['atoms']:
            raise RuntimeError(
                ('Number of atoms in "header" section does not '
                 'correspond to number of atoms in "Atoms" section.'))

        if at_type > self._type_counters['atom']:
            raise RuntimeError(('Atom type {} not found.'.format(at_type)))

        if update:  # Update
            if at_id not in self.atoms:
                raise RuntimeError('Cannot update atom with id {}. Not found.'.format(at_id))
            update_dict = {
                'position': (
                    x*self.distance_scale_factor - self._box_translate['x'],
                    y*self.distance_scale_factor - self._box_translate['y'],
                    z*self.distance_scale_factor - self._box_translate['z']),
                'atom_type': at_type,
                'res_id': at_tag,
                'charge': q
                }
            if nx is not None:
                update_dict['image'] = (nx, ny, nz)
            self.atoms[at_id].update(update_dict)
        else:  # New entry
            if at_id in self.atoms:
                raise RuntimeError('Cannot overwrite atom with id {} if update=False'.format(at_id))
            self.atoms[at_id] = {
                'atom_type': at_type,
                'res_id': at_tag,
                'position': (
                    x*self.distance_scale_factor - self._box_translate['x'],
                    y*self.distance_scale_factor - self._box_translate['y'],
                    z*self.distance_scale_factor - self._box_translate['z']),
                'image': (nx, ny, nz),
                'charge': q,
                'vel': (0.0, 0.0, 0.0),
                'mass': self._mass_type.get(at_type, 0.0)
            }

    def _read_velocity(self, input_line):
        sp_line = input_line.split()
        sp_line[0] = int(sp_line[0])
        sp_line[1:] = map(float, sp_line[1:])
        at_id, vx, vy, vz = sp_line
        self.atoms[at_id]['vel'] = (
            vx*self.distance_scale_factor,
            vy*self.distance_scale_factor,
            vz*self.distance_scale_factor)

    def _read_bond(self, input_line):
        idd, bond_type, at_1, at_2 = map(int, input_line.split())
        if idd > self._item_counters['bonds']:
            raise RuntimeError('Number of bond is wrong.')

        if at_1 not in self.atoms or at_2 not in self.atoms:
            raise RuntimeError('{} or {} not found in list of atoms.'.format(at_1, at_2))
        self.topology['bonds'][bond_type].append(tuple(sorted((at_1, at_2))))

    def _read_angle(self, input_line):
        idd, angle_type, at_1, at_2, at_3 = map(int, input_line.split())
        if idd > self._item_counters['angles']:
            raise RuntimeError('Number of angle is wrong.')

        if at_1 not in self.atoms or at_2 not in self.atoms or at_3 not in self.atoms:
            raise RuntimeError('{}, {} or {} not found in list of atoms.'.format(at_1, at_2, at_3))
        self.topology['angles'][angle_type].append((at_1, at_2, at_3))

    def _read_dihedral(self, input_line):
        idd, dihedral_type, at_1, at_2, at_3, at_4 = map(int, input_line.split())
        if idd > self._item_counters['dihedrals']:
            raise RuntimeError('Number of dihedrals is wrong.')

        if (at_1 not in self.atoms or at_2 not in self.atoms or at_3 not in self.atoms or
                at_4 not in self.atoms):
            raise RuntimeError('{}, {}, {} or {} not found in list of atoms.'.format(
                at_1, at_2, at_3, at_4))
        self.topology['dihedrals'][dihedral_type].append((at_1, at_2, at_3, at_4))

    def _read_improper(self, input_line):
        idd, dihedral_type, at_1, at_2, at_3, at_4 = map(int, input_line.split())
        if idd > self._item_counters['improperss']:
            raise RuntimeError('Number of impropers is wrong.')

        if (at_1 not in self.atoms or at_2 not in self.atoms or at_3 not in self.atoms or
                at_4 not in self.atoms):
            raise RuntimeError('{}, {}, {} or {} not found in list of atoms.'.format(
                at_1, at_2, at_3, at_4))
        self.topology['impropers'][dihedral_type].append((at_1, at_2, at_3, at_4))

    def _read_mass(self, input_line):
        sp_line = input_line.split()
        at_id, mass = sp_line
        self._mass_type[int(at_id)] = float(mass)


def read_coordinates(file_name):
    file_suffix_class = {
        'pdb': PDBFile,
        'gro': GROFile
        }
    file_suffix = file_name.split('.')[-1]
    return file_suffix_class[file_suffix](file_name)


def read_topology(file_name, settings):
    file_suffix_class = {
        'top': GROMACSTopologyFile,
        }
    file_suffix = file_name.split('.')[-1]
    return file_suffix_class[file_suffix](file_name, settings)
