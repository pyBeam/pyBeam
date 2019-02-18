/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Rocco Bombardieri, Tim Albring, Ruben Sanchez, Rauno Cavallaro
 * 
 * Developers: Tim Albring, Ruben Sanchez (SciComp, TU Kaiserslautern)
 *             Rocco Bombardieri, Rauno Cavallaro (Carlos III University Madrid)
 *
 * This file is part of pyBeam.
 *
 * pyBeam is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * pyBeam is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Affero General Public License for more details.
 * You should have received a copy of the GNU Affero
 * General Public License along with pyBeam.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once
#include <math.h>
#include <map>
#include <stdio.h>      /* printf, fopen */
/*We can think of using #include "./mpi_structure.hpp" and the command SU2_MPI::Error
 * once fully integrated with SU2 core
 */ 
#include <stdlib.h>     /* exit, EXIT_FAILURE */ 

#include "../include/types.h"





class CProperty
{
    
private:
    
protected:
    
    // Units Sys: SI
    
    
    addouble A;				// cross section area
    addouble Iyy, Izz;
    addouble Jt; 				//Torsional Moment of Inertia
    addouble J0; 				//Polar Moment of Inertia
    unsigned long PropertyID  ;        
    
    
public:
    
    
    
    CProperty(int ID) ;
    
    ~CProperty(void);
    
    //void SetWebThickness(passivedouble thickness) {t = thickness; }
    
    //void SetWebHeight(passivedouble height) {h = height; }  
    
    //void SetFlangeWidth(passivedouble FlangeWidth) {b = FlangeWidth; } 
    
    void SetSectionProperties(passivedouble A_in, passivedouble Iyy_in, passivedouble Izz_in, passivedouble Jt_in) ; 
    
    //void SetAuxVect(passivedouble x, passivedouble y, passivedouble z ) {AuxVector(0) = x; AuxVector(1) = y; AuxVector(2) = z; };
    
    //addouble GetWebThickness(void) { return t;}
    
    //addouble GetWebHeight(void) { return h;}
    
    //addouble GetFlangeWidth(void) { return b;}  
    
    addouble GetIyy(void) { return Iyy; }     
    
    addouble GetIzz(void) { return Izz; } 
    
    addouble GetA(void) { return A; }   
    
    addouble GetJt(void) { return Jt; }
    
    addouble GetJ0(void) { return J0; }  
    
    //Vector3dDiff GetAuxVect { return AuxVector;}
    
};