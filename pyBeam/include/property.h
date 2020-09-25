/*
 * pyBeam, an open-source Beam Solver
 *
 * Copyright (C) 2019 by the authors
 *
 * File developers: Rocco Bombardieri (Carlos III University Madrid)
 *                  Ruben Sanchez (SciComp, TU Kaiserslautern)
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
#include <stdio.h>
#include <stdlib.h>


#include "../include/types.h"


class CProperty{

//private:

protected:

    // Units Sys: SI
    addouble A;                 /// cross section area
    addouble Iyy, Izz;          /// Bending Moment of Inertia
    addouble Jt;                /// Torsional Moment of Inertia
    addouble J0;                /// Polar Moment of Inertia
    
    unsigned long PropertyID;   // ID of the property

    addouble t_sk;           // skin  thickness 
    addouble t_sp;           // spar thickness
    addouble A_stiff;        // stiffener Area
    addouble A_fl;           // flanges Area 
    addouble h;              //  box tot height 
    addouble C_wb;           //  box tot length 
    int n_stiff;             // number of stiffener (not taking into account the Spar's flanges)
    
    
public:
    
//    addouble A_2;               ///< cross section area
//    addouble Iyy_2, Izz_2;      ///< Bending Moment of Inertia
//    addouble Jt_2;              ///< Torsional Moment of Inertia
//    addouble J0_2;              ///< Polar Moment of Inertia
//    addouble Iyy_b;
//    addouble Izz_b;
    
    
     
    CProperty(int ID) {PropertyID = ID; }

    ~CProperty(void){}
     
    
    inline void SetSectionProperties(passivedouble A_in, passivedouble Iyy_in,
                                     passivedouble Izz_in, passivedouble Jt_in) {
        A = A_in;
        Iyy = Iyy_in;
        Izz = Izz_in;
        J0 = Iyy + Izz;
        Jt = Jt_in;               
    }
    
    /*passivedouble tsk_in,passivedouble tsp_in,
                                        passivedouble Astiff_in,passivedouble h_in,
                                        passivedouble Cwb_in,passivedouble nstiff_in,passivedouble b_in,VectorXdDiff ys_in
     */
    void SetSectionProperties2();
 

    //inline addouble Get_t_sk(void){return t_sk;}
    
    //inline addouble Get_t_sp(void){return t_sp;}
    
   // inline addouble GetA_stiff(void){return A_stiff;}
    
    //inline addouble Get_h(void){return h;}
    
    //inline addouble GetC_wb(void){return C_wb;}
    
    //inline addouble Get_n_stiff(void){return n_stiff;}
    
    //inline addouble Get_b(void){return b;}
    
    //inline VectorXdDiff Get_ys(void){return ys;}
     
    
    inline addouble GetIyy(void) { return Iyy; }

    inline addouble GetIzz(void) { return Izz; }

    inline addouble GetA(void) { return A; }

    inline addouble GetJt(void) { return Jt; }

    inline addouble GetJ0(void) { return J0; }
    
    
 
};
