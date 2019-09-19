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
                        
pyBeamAD = static_library('pyBeamAD', 
                           beam_src,
                           install : false, 
                           dependencies : [eigen_dep, codi_dep],
                           cpp_args : ['-Wno-ignored-attributes', '-Wno-deprecated-declarations', codi_rev_args])                        

pyBeam_dep = [declare_dependency(link_with: pyBeam)]
pyBeamAD_dep = [declare_dependency(link_with: pyBeamAD)]

cpp_source = swig_gen.process('./swig/pyBeam.i')
swig_lib = shared_library(
    '_pyBeam',
    cpp_source,
    dependencies: [pyBeam_dep, python_dep, eigen_dep],
    link_args: [],
    install: true,
    include_directories : ['./include/'],
    cpp_args : ['-Wno-ignored-attributes', '-Wno-deprecated-declarations', '-DSWIG_TYPE_TABLE=pyBeam'],
    name_prefix : '',
)
meson.add_install_script(meson.source_root()+'/install.sh', meson.current_build_dir(), 'pyBeam.py')

cpp_source = swig_gen.process('./swig/pyBeamAD.i')
swig_lib = shared_library(
    '_pyBeamAD',
    cpp_source,
    dependencies: [pyBeamAD_dep, python_dep, codi_dep, eigen_dep],
    link_args: [],
    install: true,
    include_directories : ['./include/'],
    cpp_args : ['-Wno-ignored-attributes', 
                '-Wno-deprecated-declarations', 
                '-DSWIG_TYPE_TABLE=pyBeamAD', 
                codi_rev_args],
    name_prefix : '',
)
meson.add_install_script(meson.source_root()+'/install.sh', meson.current_build_dir(), 'pyBeamAD.py')
