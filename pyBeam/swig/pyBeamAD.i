/* file : beam.i */

%import "../include/types.h"

/* name of module to use*/
%module pyBeamAD
%{
    #include "../include/beam.h"
    #include "../include/input.h"
    #include "../include/property.h"
    #include "../include/geometry.h"
    #include "../include/element.h"
    #include "../include/rigid_element.h"
    #include "../include/CDV.h"
%}

%include "std_string.i"
%include "../include/beam.h";
%include "../include/input.h";
%include "../include/property.h";
%include "../include/geometry.h";
%include "../include/element.h";
%include "../include/rigid_element.h";
%include "../include/CDV.h";
