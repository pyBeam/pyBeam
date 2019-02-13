/* file : beam.i */

%import "../include/types.h"
 
/* name of module to use*/
%module pyBeamFM
%{
    #include "../include/beam.h"
%}

%include "../include/beam.h";

