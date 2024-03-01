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
#include <math.h>
#include <map>
#include <stdio.h>
#include <stdlib.h>

#include "../include/types.h"


class CInput
{

private:

protected:

    //##################     Numerical Inputs     ###########################

    int nConstr;                    // Total number of constraints
    MatrixXdDiff  Constr_matrix;    // COnstraint matrix [ NODE_ID DOF_ID ]
    
    unsigned long nNodes;           // Number of overall nodes along the wing (no collapsed)
    unsigned long nFEM;             // Number of finite elements
    unsigned long nRBE2;            // Number of RBE2 rigid elements

    unsigned short nDOF;            // Number of degrees of freedom

    addouble penalty = pow(10,8);  // penalty coefficient

    addouble load;                  // [N];              // To be removed
    int follower_flag;              // (0) Nonfollower (1) follower (2) approx follower
    unsigned long loadSteps;        // Number of load steps
    unsigned long nIter;            // Number of iterations

    addouble end_time;              // [sec] for SS calculation                    // To be removed
    addouble dt;                    // [sec] time increment for SS calculation     // To be removed

    addouble tol_LinSol;            // Tolerance of the linear solver
    unsigned short kind_LinSol;     // Tolerance of the linear solver
    unsigned short WriteRestart;    // Write restart option
    unsigned short Restart;         // Restart option

    bool discrete_adjoint;

    //##############    Material inputs (only ONE homogeneous material is allowed by now)  ###########################
    // Units Sys: SI
    addouble E;                     // Elastic modulus [GPa]
    addouble E_dimensional;         // Dimensional Elastic modulus [GPa]
    addouble Poiss;                 // Poisson Ratio
    addouble ro;                    // Beam Density [kg/m^3]
    addouble G;                     // Shear modulus

    //################     Convergence Parameters    ###########################
    addouble convCriteria;


public:

    
    // Double constructor to make it compatible for old cases in which RBE2 were not defined
    CInput(int py_nPoint, int py_nElem);
    
    CInput(int py_nPoint, int py_nElem, int py_nRBE2);

    virtual ~CInput(void);

    void SetParameters();

    void SetDiscreteAdjoint(void) {discrete_adjoint = true; }

    bool GetDiscreteAdjoint(void) {return discrete_adjoint; }

    void SetYoungModulus(passivedouble YoungModulus) {E_dimensional = YoungModulus;}

    void SetPoisson(passivedouble Poisson) {Poiss = Poisson;}

    void RegisterInput_E(void) {AD::RegisterInput(E_dimensional);}

    void RegisterInput_Nu(void) {AD::RegisterInput(Poiss);}

    passivedouble GetGradient_E(void) {return AD::GetValue(AD::GetDerivative(E_dimensional));}

    passivedouble GetGradient_Nu(void) {return AD::GetValue(AD::GetDerivative(Poiss));}

    void SetDensity(passivedouble Density) {ro = Density; }

    void SetFollowerFlag(int FollowerFlag) {follower_flag = 0; /*Forced to be nonfollower*/ }

    void SetLoadSteps(unsigned long LoadSteps) { loadSteps = LoadSteps; }

    void SetNStructIter(unsigned long NStructIter) {nIter = NStructIter; }

    void SetTolerance_LinSol(passivedouble tolerance) {tol_LinSol = tolerance; }

    void SetKind_LinSol(unsigned short kind_solver) {kind_LinSol = kind_solver; }

    void Set_WriteRestartFlag(unsigned short WriteRestartFlag) {WriteRestart = WriteRestartFlag; }

    void Set_RestartFlag(unsigned short RestartFlag) {Restart = RestartFlag; }

    void SetConvCriterium(passivedouble ConvCriterium) {convCriteria = ConvCriterium; }

    void SetnConstr(int nconstr) {nConstr = nconstr; Constr_matrix = MatrixXdDiff::Zero(nconstr,2);}

    void SetSingleConstr(int iConstr, int node_id, int DOF_id) {Constr_matrix(iConstr,1 -1) = node_id;
                                                                Constr_matrix(iConstr,2 -1) = DOF_id;}

    MatrixXdDiff  GetConstrMatrix() {return Constr_matrix;};

    addouble GetYoungModulus() {return E; }

    addouble GetYoungModulus_dimensional() {return E_dimensional; }

    addouble GetPoisson() {return Poiss; }

    void SetShear(addouble val_shear) {G = val_shear;}

    addouble GetShear() {return G; }

    addouble GetDensity() {return ro; }

    addouble GetTolerance_LinSol(void) {return tol_LinSol; }

    unsigned short Get_WriteRestartFlag(void) {return WriteRestart; }

    unsigned short Get_RestartFlag(void) {return Restart; }

    unsigned short GetKind_LinSol(void) {return kind_LinSol; }

    unsigned long Get_nNodes(void) { return nNodes; }

    unsigned long Get_nFEM(void) { return nFEM; }

    unsigned long Get_nRBE2(void) { return nRBE2; }

    unsigned short Get_nDOF(void) { return nDOF; }

    unsigned short Get_FollowerFlag(void) { return follower_flag; }

    unsigned long Get_LoadSteps(void) { return loadSteps;}

    unsigned long Get_nIter(void) { return nIter;}

    addouble Get_ConvCriteria(void) { return convCriteria; }

    addouble GetPenalty(void) { return penalty; }


};


