# pyBeam: an open-source Beam Solver
### Beta version: 0.1

pyBeam is a nonlinear beam finite element solver developed with aeronautical design applications in mind.

It incorporates an AD-based adjoint solver for gradient computation, which relies on the AD library [CoDiPack](https://www.scicomp.uni-kl.de/codi/). Its goal aims towards a fully-functional adjoint-based infrastructure for performing gradient-based optimization of aircraft wing configurations, via coupling with the open-source CFD Suite [SU2](https://su2code.github.io). For an overview of the technical details in pyBeam, please see the following paper presented at EUROGEN 2019:

Bombardieri, R.,_et al_(2019) ["Towards an open-source framework for aero-structural design and optimization within the SU2 suite"](https://www.researchgate.net/publication/335972259_Towards_an_open-source_framework_for_aero-structural_design_and_optimization_within_the_SU2_suite) EUROGEN 2019 Proceedings, Guimarães, Portugal, Sept 12-14

Please also kindly cite our paper above if you are using pyBeam in your research!

## pyBeam Developers

We follow the popular "GitFlow" branching model for scalable development. In this repository, the master branch represents the latest stable major or minor release, and it should only be modified during version releases.

Work that is staged for release is merged into the develop branch via pull requests from various "feature" branches. At release time, the develop branch is pushed into the master branch and tagged as a release.

pyBeam is open for development for anyone wishing to contribute. A list of current contributors can be found in the AUTHORS.md file.

Continuous integration of test cases is provided by Travis CI.

[![Build Status](https://travis-ci.com/pyBeam/pyBeam.svg?branch=develop)](https://travis-ci.com/pyBeam/pyBeam)

## Installation tips:

NOTE: this procedure has been verified for Ubuntu 18.04

- Make sure you have meson installed

 Using pip (from terminal):
```
pip3 install --upgrade meson
```

From package manager:
```
sudo apt-get install meson
```

- Make sure you have Swig installed
```
sudo apt-get update
sudo apt-get install swig
```

- Make sure to initialize the external submodules
```
git submodule init
git submodule update
```

- Compile:
```
meson build --prefix=$PWD
```

- Compile and install (into prefix folder):
```
ninja -C build install
```

IMPORTANT: If compilation fails check your meson version. Compiling only works for version 0.52.0 (verified) and above. 
If compiling returns error it may be due to an old version of meson.

- Upload meson (via pip)
```
pip3 install --upgrade meson
```

- Check meson version
```
meson -v
```

If uninstall is necessary
```
sudo apt-get remove meson
```
or, via pip,
```
pip uninstall meson
```
and install it again.

