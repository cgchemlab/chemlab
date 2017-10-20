"""
Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>

This file is part of ChemLab.

ChemLab is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import setuptools.command.test
from setuptools import setup, find_packages

class TestCommand(setuptools.command.test.test):
    """ Setuptools test command explicitly using test discovery. """

    def _test_args(self):
        yield 'discover'
        for arg in super(TestCommand, self)._test_args():
            yield arg

setup(
    name='chemlab',
    author='Jakub Krajniak',
    author_email='jkrajniak@gmail.com',
    version='1.5.0',
    license='GPLv3+',
    test_suite='src.tests')
