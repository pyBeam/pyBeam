#
# pyBeam, an open-source Beam Solver
#
# Copyright (C) 2019 by the authors
#
# File Developers: Ruben Sanchez (SciComp, TU Kaiserslautern)
#                  Tim Albring (SciComp, TU Kaiserslautern)
#
# This file is part of pyBeam.
#
# pyBeam is free software: you can redistribute it and/or
# modify it under the terms of the GNU Affero General Public License
# as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# pyBeam is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU Affero General Public License for more details.
# You should have received a copy of the GNU Affero
# General Public License along with pyBeam.
# If not, see <http://www.gnu.org/licenses/>.
#

project('pyBeam', 'cpp', 
        version: '0.1.0', 
        license: 'GNU Affero 3.0',
        default_options: ['buildtype=release', 
                          'warning_level=0', 
                          'cpp_std=c++11'])

pymod = import('python')
python = pymod.find_installation()
python_dep = dependency('python3-embed')

# check for swig
swig = find_program('swig')
swig_gen = generator(
    swig,
    output: ['@BASENAME@.cxx'],
    arguments: [ '-c++', '-python', '-outdir', meson.current_build_dir(), '-o', './@OUTPUT@', '@INPUT@'],
    depfile: '@BASENAME@.d',
)

eigen_dep = [declare_dependency(include_directories: './pyBeam/externals/Eigen')]

pyBeam_cpp_args = []
pyBeam_deps     = []

# compile common library
subdir('pyMLS/externals/ann')
subdir('pyMLS')
subdir('pyMLS/pylib')
subdir('pyBeam')
subdir('pyBeam/pylib')

