#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
#
#  This file is part of ChemLab
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

import os
import unittest
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))
print sys.path[0]

import chemlab

class TestTopologyReader(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.topol_file = 'topol.top'
        cls.gt = chemlab.gromacs_topology.GromacsTopology(cls.topol_file, generate_exclusions=True)
        cls.gt.read()

    def test_replicated_molecules(self):
        """Test the molecule replication"""
        total_nr_atoms = len(self.gt.atoms)
        expected_nr_atoms = 0
        for mol_name, nmols in self.gt.gt.molecules:
            mol_atoms = len(self.gt.gt.molecules_data[mol_name]['atoms'])
            expected_nr_atoms += nmols * mol_atoms
        self.assertEqual(total_nr_atoms, expected_nr_atoms)
        
        total_nr_bonds = len(self.gt.bonds)
        expected_nr_bonds = 0
        for mol_name, nmols in self.gt.gt.molecules:
            mol_bonds = len(self.gt.gt.molecules_data[mol_name].get('bonds', []))
            expected_nr_bonds += nmols * mol_bonds
        self.assertEqual(total_nr_bonds, expected_nr_bonds)

        total_nr_angles = len(self.gt.angles)
        expected_nr_angles = 0
        for mol_name, nmols in self.gt.gt.molecules:
            mol_angles = len(self.gt.gt.molecules_data[mol_name].get('angles', []))
            expected_nr_angles += nmols * mol_angles
        self.assertEqual(total_nr_angles, expected_nr_angles)
        
        total_nr_dihedrals = len(self.gt.dihedrals)
        expected_nr_dihedrals = 0
        for mol_name, nmols in self.gt.gt.molecules:
            mol_dihedrals = len(self.gt.gt.molecules_data[mol_name].get('dihedrals',[]))
            expected_nr_dihedrals += nmols * mol_dihedrals
        self.assertEqual(total_nr_dihedrals, expected_nr_dihedrals)

        total_nr_pairs = len(self.gt.pairs)
        expected_nr_pairs = 0
        for mol_name, nmols in self.gt.gt.molecules:
            mol_pairs = len(self.gt.gt.molecules_data[mol_name].get('pairs', []))
            expected_nr_pairs += nmols * mol_pairs
        self.assertEqual(total_nr_pairs, expected_nr_pairs)


if __name__ == '__main__':
    unittest.main()
