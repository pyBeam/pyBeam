/* file : beam.i */

%import "../include/types.h"
 
/* name of module to use*/
%module pyBeam
%{
    #include "../include/beam.h"
%}

%include "../include/beam.h";

