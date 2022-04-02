/*
 * pyBeam, an open-source Beam Solver
 *
 * Copyright (C) 2019 by the authors
 *
 * File developers: Rocco Bombardieri (Carlos III University Madrid)
 *                  Rauno Cavallaro (Carlos III University Madrid)
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
#include <math.h>       /* exp */

#include "../include/types.h"

//#include "../include/element.h"
//#include "../include/rigid_element.h"
//#include "../include/structure.h"
//#include "../include/geometry.h"
//#include "../include/input.h"

class CDV{

private:

//    int ID;                    ///< Global ID number of the DV
   
    std::string TAG;            ///< KEY in the input file 
    int idx;                  ///< position in the container of the entities defined by TAG
    std::string sTAG;         ///< sub KEY within the field     
    
    addouble lB;              ///< lower bound
    addouble uB;              ///< upper bound


protected:

public:

    CDV(void){};

    virtual ~CDV(void){};

    inline void SetDV(std::string TAGi, int idxi, std::string sTAGi,  
                       passivedouble lBi, passivedouble uBi) {
        TAG = TAGi;
        idx = idxi;
        sTAG =  sTAGi;
        lB   = lBi;
        uB = uBi;
    }

    inline std::string GetTAG(void) { return TAG; }
    inline int         Getidx(void) { return idx; }
    inline std::string GetsTAG(void) { return sTAG; }    
    inline addouble GetlB(void) { return lB; }
    inline addouble GetuB(void) { return uB; }   
    
};
