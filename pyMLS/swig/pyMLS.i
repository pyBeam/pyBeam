/* file : pyMLS.i */
 
%module pyMLS
%{
    #include "../include/interface.hpp"
%}

%include "std_vector.i"

namespace std {
    %template(DoubleVector)  vector<double>;
};        

%include "../include/interface.hpp";
