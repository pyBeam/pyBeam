project('pyBeam', 'cpp', 
        version: '0.1.0', 
        license: 'GNU Affero 3.0',
        default_options: ['buildtype=release', 
                          'warning_level=0', 
                          'cpp_std=c++11'])

pymod = import('python')
python = pymod.find_installation()
python_dep = dependency('python3')

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

