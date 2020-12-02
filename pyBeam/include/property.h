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
    addouble A;                 ///< cross section area
    addouble Iyy, Izz;          ///< Bending Moment of Inertia
    addouble Jt;                ///< Torsional Moment of Inertia
    addouble J0;                ///< Polar Moment of Inertia   
    unsigned long PropertyID;   // ID of the property
    
    addouble A_b;               ///< Cross section area of the booms
    addouble Iyy_b, Izz_b;      ///<    
    addouble C_wb;           ///<  box tot length 
    addouble h;              ///<  box tot height 
    addouble t_sk;           ///< skin  thickness 
    addouble t_sp;           ///< spar thickness
    addouble A_fl;           ///< flanges Area     
    int n_stiff;             ///< number of stiffener (not taking into account the Spar's flanges)
    addouble A_stiff;        ///< stiffener Area    
    
    int isWBDV=0;     ///< flag which determines if the initial inputs were given as WB sizes
    
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
        // If only inertial properties are given, then we assume there are no stringers 
        A_b = 0.0;
        Iyy_b = 0.0;
        Izz_b = 0.0;
        isWBDV = 0;
        // Set to zero the WB sizes 
        C_wb = 0;          
        h = 0;
        t_sk = 0;
        t_sp = 0;
        A_fl = 0;
        n_stiff = 0 ;
        A_stiff = 0;        
    
    }
    
    void SetSectionProperties(passivedouble C_wb_, passivedouble h_, passivedouble t_sk_,  
                              passivedouble t_sp_, passivedouble A_fl_, 
                              int n_stiff_, passivedouble A_stiff_);    
    
    void FromWBtoInertias();
    
    
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
     
    inline int Getn_stiff(void) { return n_stiff; }
     
    inline addouble GetA_stiff(void) { return A_stiff; }
    
    inline int GetisWBDV(void){ return  isWBDV;} 
     
    void RegisterInput_WB(void);      ///<  Registers properties as inputs for sensitivity evaluation

    void GetGradient_WB(void);      ///<  Registers properties as inputs for sensitivity evaluation
    
    void RegisterInput_A(void){
//        if (isWBDV == 1){        
//           AD::RegisterInput(C_wb);}
        if (isWBDV == 0){
           AD::RegisterInput(A);}}        
    
    
    passivedouble  GetGradient_A(void){
//        if (isWBDV == 1){
//            return AD::GetValue(AD::GetDerivative(C_wb));}
        if (isWBDV == 0){
            return AD::GetValue(AD::GetDerivative(A));} }       
    
    void InitializePropDVsVec(addouble*,int);   ///< Initializing prop DVs Vecor from Property
    
    void SetDependencyfromDVVec(addouble*,int); ///< Register the dependency from the prop DVs Vector
    
    void SetA(addouble A_in){A=A_in;} ;
};