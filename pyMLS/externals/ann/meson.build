ann_src =files(['./src/ANN.cpp',
                './src/brute.cpp',
                './src/kd_tree.cpp',
                './src/kd_util.cpp',
                './src/kd_split.cpp',
                './src/kd_search.cpp',
                './src/kd_pr_search.cpp',
                './src/kd_fix_rad_search.cpp',
                './src/kd_dump.cpp',
                './src/bd_tree.cpp',
                './src/bd_search.cpp',
                './src/bd_pr_search.cpp',
                './src/bd_fix_rad_search.cpp',
                './src/perf.cpp'])

libANN = static_library('ANN', 
                        ann_src,
                        install : false, 
                        include_directories : ['./include'],
                        cpp_args : ['-Wno-ignored-attributes', '-Wno-deprecated-declarations', '-O3'])

ann_dep = declare_dependency(include_directories : ['./include'],
                            link_with: libANN)

