# pyBeam, an open-source Beam Solver
#
# Copyright (C) 2019 by the authors
#
# File Developer: Rocco Bombardieri (Carlos III University Madrid)
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
# Installation tips:

NOTE: this procedure has been verified for Ubuntu and Fedora

- Make sure you have meson installed:
a) Using pip (from terminal):
pip3 install --upgrade meson

b) Regular (from terminal):

i) For Ubuntu
sudo apt-get install meson
ii) For Fedora
sudo dnf meson

- Make sure yopu have Swig installed:
sudo apt-get update
sudo apt-get install swig

- Make sure to initialize the external submodules:
git submodule init
git submodule update

- Compile:
meson build --prefix=$PWD

- Install (into prefix folder):
ninja -C build install

IMPORTANT: If compilation fails check your meson version. Compiling only works for version 0.52.0 (verified) and above. 
If compiling returns error it may be due to an old versino of meson. 

-Upload meson (via pip):
pip3 install --upgrade meson

- Check meson version
meson -v

If unistall is necessary:
sudo apt-get remove meson

or (via pip):
pip uninstall somepackage

and istall it again.

