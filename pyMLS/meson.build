dep_dirs = include_directories(['./externals/ann',
                                './externals/libigl',
                                './externals/ann/lib',
                                '../pyBeam/externals/Eigen'])
lib_dirs = [meson.current_source_dir()+'/externals/ann/lib/libANN.a']
mls_src =files(['./src/interface.cpp'])

dep_declar = declare_dependency(include_directories : dep_dirs)
python_dep = dependency('python3')
pyMLS = static_library('pyMLS', 
                        mls_src,
                        install : false, 
                        dependencies : dep_declar, 
                        cpp_args : ['-Wno-ignored-attributes', '-Wno-deprecated-declarations'])

pyMLS_dep = [declare_dependency(link_with: pyMLS)]

swig_gen = generator(
    swig,
    output: ['@BASENAME@.cxx'],
    arguments: [ '-c++', '-python', '-outdir', meson.current_build_dir(), '-o', './@OUTPUT@', '@INPUT@'],
    depfile: '@BASENAME@.d',
)


cpp_source = swig_gen.process('./swig/pyMLS.i')
swig_lib = shared_library(
    '_pyMLS',
    cpp_source,
    dependencies: [pyMLS_dep, python_dep],
    link_args: [lib_dirs],
    install: true,
    include_directories : ['./include/', dep_dirs],
    cpp_args : ['-Wno-ignored-attributes', '-Wno-deprecated-declarations'],
    name_prefix : '',
)
meson.add_install_script(meson.source_root()+'/install.sh', meson.current_build_dir(), 'pyMLS.py')
