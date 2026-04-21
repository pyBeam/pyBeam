/* file : pyMLS.i */
 
%module pyMLS
%{
    #include "../include/interface.hpp"
%}

%include "std_vector.i"

namespace std {
    %template(DoubleVector) vector<double>;
    %template(IntVector)    vector<int>;
};        

%include "../include/interface.hpp";
