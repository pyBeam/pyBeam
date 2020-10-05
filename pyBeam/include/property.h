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
    addouble A_b;
    addouble Iyy_b, Izz_b; 
    
    unsigned long PropertyID;   // ID of the property

    addouble C_wb;           //  box tot length 
    addouble h;              //  box tot height 
    addouble t_sk;           // skin  thickness 
    addouble t_sp;           // spar thickness
    addouble A_fl;           // flanges Area     
    int n_stiff;             // number of stiffener (not taking into account the Spar's flanges)
    addouble A_stiff;        // stiffener Area    
    
public:
    
        
    CProperty(int ID) {PropertyID = ID; }

    ~CProperty(void){}
     
    
    inline void SetSectionProperties(passivedouble A_in, passivedouble Iyy_in,
                                     passivedouble Izz_in, passivedouble Jt_in, passivedouble A_b_in, passivedouble Iyy_b_in, passivedouble Izz_b_in) {
        A = A_in;
        Iyy = Iyy_in;
        Izz = Izz_in;
        J0 = Iyy + Izz;
        Jt = Jt_in;
        A_b=A_b_in;
        Iyy_b= Iyy_b_in;
        Izz_b= Izz_b_in;
    }
    
    /*passivedouble tsk_in,passivedouble tsp_in,
                                        passivedouble Astiff_in,passivedouble h_in,
                                        passivedouble Cwb_in,passivedouble nstiff_in,passivedouble b_in,VectorXdDiff ys_in
     */
    void SetSectionProperties(passivedouble C_wb_, passivedouble h_, passivedouble t_sk_,  
                              passivedouble t_sp_, passivedouble A_fl_, 
                              int n_stiff_, passivedouble A_stiff_);
 
    void SetSectionProperties2();
    
    //void VectorXdDiff GetDesignVariable(void);
         
    inline addouble GetIyy(void) { return Iyy; }

    inline addouble GetIzz(void) { return Izz; }

    inline addouble GetA(void) { return A; }

    inline addouble GetJt(void) { return Jt; }

    inline addouble GetJ0(void) { return J0; }
     
     inline addouble GetA_b(void) { return A_b; }
      
     inline addouble GetIyy_b(void) { return Iyy_b; }
       
     inline addouble GetIzz_b(void) { return Izz_b; }
     
     inline addouble GetC_wb(void) { return C_wb; }
     
     inline addouble Geth(void) { return h; }
     
     inline addouble Gett_sk(void) { return t_sk; }
     
     inline addouble Gett_sp(void) { return t_sp; }
     
     inline addouble GetA_fl(void) { return A_fl; }
     
     inline addouble Getn_stiff(void) { return n_stiff; }
     
     inline addouble GetA_stiff(void) { return A_stiff; }
};
