#
# pyMLS, an open-source Moving Least Squares library
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

dep_dirs = include_directories(['./externals/libigl'])
mls_src =files(['./src/interface.cpp'])

mls_dep = declare_dependency(include_directories : dep_dirs)


pyMLS = static_library('pyMLS', 
                        mls_src,
                        install : false, 
                        dependencies : [mls_dep, eigen_dep, ann_dep],
                        cpp_args : ['-Wno-ignored-attributes', '-Wno-deprecated-declarations'])

pyMLS_dep = [declare_dependency(link_with: pyMLS)]

cpp_source = swig_gen.process('./swig/pyMLS.i')
swig_lib = shared_library(
    '_pyMLS',
    cpp_source,
    dependencies: [pyMLS_dep, python_dep, eigen_dep, mls_dep],
    install: true,
    include_directories : ['./include/'],
    install_dir: ['lib'],
    cpp_args : ['-Wno-ignored-attributes', '-Wno-deprecated-declarations'],
    name_prefix : '',
)
meson.add_install_script(meson.source_root()+'/install.sh', meson.build_root(), 'pyMLS.py')
