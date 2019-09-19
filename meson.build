project('pyBeam', 'cpp', 
        version: '0.1.0', 
        license: 'GNU Affero 3.0',
        default_options: ['buildtype=release', 
                          'warning_level=0', 
                          'cpp_std=c++11'])

pymod = import('python')
python = pymod.find_installation()

# check for swig
swig = find_program('swig')

pyBeam_cpp_args = []
pyBeam_deps     = []

# compile common library
subdir('pyMLS')
subdir('pyBeam')

