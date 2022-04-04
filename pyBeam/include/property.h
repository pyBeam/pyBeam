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

private:

protected:

    // Units Sys: SI
    addouble A;                 // cross section area
    addouble Iyy, Izz;          // Bending Moment of Inertia
    addouble Jt;                // Torsional Moment of Inertia
    addouble J0;                // Polar Moment of Inertia
    unsigned long PropertyID;   // ID of the property

public:

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

    inline addouble GetIyy(void) { return Iyy; }

    inline addouble GetIzz(void) { return Izz; }

    inline addouble GetA(void) { return A; }

    inline addouble GetJt(void) { return Jt; }

    inline addouble GetJ0(void) { return J0; }
    
    
    inline void SetIyy(addouble Iyy_) { Iyy = Iyy_ ; }

    inline void SetIzz(addouble Izz_) { Izz = Izz_ ; }
    
    inline void SetA(addouble A_) { A = A_ ; }

    inline void SetJt(addouble Jt_) { Jt = Jt_ ; }    
  

};
