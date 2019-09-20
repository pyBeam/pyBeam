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
