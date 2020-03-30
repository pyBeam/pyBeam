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

beam_src =files(['./src/ad.cpp',
                './src/input.cpp',
                './src/element.cpp',
                './src/rigid_element.cpp',
                './src/structure.cpp',
                './src/rotations.cpp',
                './src/beam.cpp'])
                
codi_dep = [declare_dependency(include_directories: 'externals/CoDiPack/include')]
codi_rev_args = '-DCODI_REVERSE_TYPE'

pyBeam = static_library('pyBeam', 
                        beam_src,
                        install : false, 
                        dependencies : eigen_dep, 
                        cpp_args : ['-Wno-ignored-attributes', '-Wno-deprecated-declarations'])
                        
#pyBeamAD = static_library('pyBeamAD', 
#                           beam_src,
#                           install : false, 
#                           dependencies : [eigen_dep, codi_dep],
#                           cpp_args : ['-Wno-ignored-attributes', '-Wno-deprecated-declarations', codi_rev_args])                        

pyBeam_dep = [declare_dependency(link_with: pyBeam)]
#pyBeamAD_dep = [declare_dependency(link_with: pyBeamAD)]

cpp_source = swig_gen.process('./swig/pyBeam.i')
swig_lib = shared_library(
    '_pyBeam',
    cpp_source,
    dependencies: [pyBeam_dep, python_dep, eigen_dep],
    link_args: [],
    install: true,
    install_dir: ['bin'],
    include_directories : ['./include/'],
    cpp_args : ['-Wno-ignored-attributes', '-Wno-deprecated-declarations', '-DSWIG_TYPE_TABLE=pyBeam'],
    name_prefix : '',
)
if meson.is_subproject()
	meson.add_install_script(meson.source_root()+'/subprojects/pyBeam/install.sh', meson.build_root()+'/subprojects/pyBeam', 'pyBeam.py')
else
	meson.add_install_script(meson.source_root()+'/install.sh', meson.build_root(), 'pyBeam.py')
endif

#cpp_source = swig_gen.process('./swig/pyBeamAD.i')
#swig_lib = shared_library(
#    '_pyBeamAD',
#    cpp_source,
#    dependencies: [pyBeamAD_dep, python_dep, codi_dep, eigen_dep],
#    link_args: [],
#    install: true,
#    install_dir: ['bin'],
#    include_directories : ['./include/'],
#    cpp_args : ['-Wno-ignored-attributes', 
#                '-Wno-deprecated-declarations', 
#                '-DSWIG_TYPE_TABLE=pyBeamAD', 
#                codi_rev_args],
#    name_prefix : '',
#)
#if meson.is_subproject()
#	meson.add_install_script(meson.source_root()+'/subprojects/pyBeam/install.sh', meson.build_root()+'/subprojects/pyBeam', 'pyBeamAD.py')
#else
#	meson.add_install_script(meson.source_root()+'/install.sh', meson.build_root(), 'pyBeamAD.py')
#endif
