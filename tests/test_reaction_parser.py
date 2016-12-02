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

import espressopp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))
print sys.path[0]

import chemlab

class TestTopologyReader(unittest.TestCase):
    def test_parse_exchange_reaction(self):
        input_string = 'C(0,1):E(0,1) + W(0,1) -> A(1):Z(1) + E(1)'
        reactants, r_type = chemlab.reaction_parser.parse_exchange_equation(input_string)
        self.assertEqual(r_type, chemlab.reaction_parser.REACTION_EXCHANGE)

        # Check reactants structure
        self.assertEqual(reactants['type_1']['name'], 'C')
        self.assertEqual(reactants['type_1']['new_type'], 'A')
        self.assertEqual(reactants['type_1']['min'], '0')
        self.assertEqual(reactants['type_1']['max'], '1')
        self.assertEqual(reactants['type_1']['delta'], '1')

        self.assertEqual(reactants['type_2']['name'], 'E')
        self.assertEqual(reactants['type_2']['new_type'], 'E')
        self.assertEqual(reactants['type_2']['min'], '0')
        self.assertEqual(reactants['type_2']['max'], '1')
        self.assertEqual(reactants['type_2']['delta'], '1')

        self.assertEqual(reactants['type_3']['name'], 'W')
        self.assertEqual(reactants['type_3']['new_type'], 'Z')
        self.assertEqual(reactants['type_3']['min'], '0')
        self.assertEqual(reactants['type_3']['max'], '1')
        self.assertEqual(reactants['type_3']['delta'], '1')

if __name__ == '__main__':
    unittest.main()