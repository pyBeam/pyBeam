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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include "../include/types.h"
#include "../include/element.h"
#include "../include/rigid_element.h"
#include "../include/rotations.h"
#include "../include/geometry.h"
#include "../include/input.h"

#include <iostream>

class CStructure
{
    
private:

    bool verbose = true;
    addouble tol_LinSol;
    unsigned short kind_linSol;

    enum KIND_LINEARSOLVER {
        FullPivHouseholderQr = 0,
        PartialPivLu = 1,
        FullPivLu = 2,
        HouseholderQr = 3,
        ColPivHouseholderQr = 4,
        LLT = 5,
        LDLT = 6
    };
    
public:
    
    
    int nNode;                  // Number of structural nodes
    int nfem;                   // number of finite elements   // has to be assigned in the constructor
    int nRBE2;                  // number of finite elements
    
    int DOF;                    // In space, 6
    
    CRBE2 **RBE2;               // Pointer to the first RBE2 element
    CElement **element;         // Pointer to the first finite element
    CNode **node;               // Pointer to the first finite element

    VectorXdDiff cross_term;    // Store the displacement vector
    
    MatrixXdDiff Ksys;
    MatrixXdDiff Ksys_red;      // [relative to masters in case of RBE2]
    MatrixXdDiff K_penal;       // penalty matrix for rigid elements
    
    // In case Rigid with Lagrangian multiplier method
    MatrixXdDiff Ksys_lam;  // Aumented matrix [6*(nNode+nRBE)x6*(nNode+nRBE)]
    VectorXdDiff U_lam;     // Aumented Displacement array [6*(nNode+nRBE)]
    VectorXdDiff Residual_lam; //  // Aumented Residual array [6*(nNode+nRBE)]
    VectorXdDiff dU_lam;           // Aumented Displacement array (increment)[6*(nNode+nRBE)]  
    VectorXdDiff U_lam_adj;     // Aumented Displacement array [6*(nNode+nRBE)]
    
    MatrixXdDiff  Constr_matrix;// COnstraint matrix [ NODE_ID DOF_ID ]
    
    VectorXdDiff U;             // Displacement array
    VectorXdDiff dU;            // Displacement array (increment)
    VectorXdDiff dU_red;        // Displacement array (increment) [relative to masters in case of RBE2]
    VectorXdDiff X;             // Position of the fem nodes in global coordinate system
    VectorXdDiff X0;            // Position of the fem nodes in global coordinate system

    VectorXdDiff U_adj;         // Adjoint of the displacement array (cumulative)
    
    VectorXdDiff Fint;          // Array of internal forces
    VectorXdDiff Fext;          // Array of External Forces
    VectorXdDiff Residual;      // Array of Unbalanced Forces
    
    VectorXdDiff Fnom;          // Array of nominal forces
    VectorXdDiff Fnom_old;          // Array of nominal forces  of previous FSI run  
    
    addouble YoungModulus;
    addouble penalty;
    
    CStructure(CInput *input, CElement **container_element, CNode **container_node);
    
    ~CStructure();
    
    /*##############################################################
     *
     *         External Forces, Residual, Internal Forces
     *
     *###############################################################*/
    
    inline void ReadForces(int nTotalDOF, addouble *loadVector) {
        for (int iLoad = 0; iLoad < nTotalDOF; iLoad++){Fnom(iLoad) = loadVector[iLoad];}
    }
    
    inline void ResetForces() {
        Fnom_old = Fnom;
        Fnom = VectorXdDiff::Zero(nNode*6);
    }    

    inline void SetDimensionalYoungModulus(addouble val_E){ YoungModulus = val_E; }

    // External forces are normalized by the Young Modulus
    //inline void UpdateExtForces(addouble lambda){ Fext = lambda* Fnom / YoungModulus; }
    inline void UpdateExtForces(addouble lambda){       
        Fext = ( lambda* (Fnom -Fnom_old) +Fnom_old )  / YoungModulus;}
    void EvalResidual();

    //===================================================
    //      Assembly RBE2 rigid constraint matrix
    //===================================================

    inline void AddRBE2(CInput *input, CRBE2** container_RBE2) {nRBE2 = input->Get_nRBE2(); RBE2 = container_RBE2;}
    
    // Penalty strategy
    void SetPenalty() {  penalty = 100*Ksys.diagonal().maxCoeff(); };
    
    void RigidResidual();

    void AssemblyRigidPenalty();
    
    // Lagrange Multiplier strategy
    inline void SetRigidLagrangeDimensions() {U_lam = VectorXdDiff::Zero((nNode+nRBE2)*6); dU_lam = VectorXdDiff::Zero((nNode+nRBE2)*6); Residual_lam = VectorXdDiff::Zero((nNode+nRBE2)*6); Ksys_lam =MatrixXdDiff::Zero((nNode+nRBE2)*6,(nNode+nRBE2)*6); U_lam_adj = VectorXdDiff::Zero((nNode+nRBE2)*6); };
    
    void RigidResidualLagrange();
    
    void AssemblyRigidLagrange();
    
    void ImposeBC_RigidLagrangian();
    
    void SolveLinearStaticSystem_RigidLagrangian(int iIter, std::ofstream &history , int print); 
    //===================================================
    //      Assembly System Stiffness Matrix
    //===================================================

    void AssemblyTang(int iIter);

    void EvalSensRot();    // Evaluate the sensitivity of Rotation Matrix - need for Jacobian

    void ImposeBC();  // Impose the boundary condition on the stiffness matrix and Residual

    //===================================================
    //      Solve linear static system
    //===================================================
    // Assembles LHS and RHS and solves the linear static problem

    void SolveLinearStaticSystem(int iIter, std::ofstream &history , int print);

    //===================================================
    //      Update Coordinates
    //===================================================
    /* This member function upates the coordinates (expressed in global reference system) of
     * the finite element nodes.
     * This is necessary for "booking" the position, as the compatiblity and the equations
     * are based on the displacements.
     */

    void UpdateCoord(int nRBE2, int iRigid);

    void RestartCoord();

    void InitialCoord();

    void SetCoord0();

    void UpdateLength();

    void UpdateRotationMatrix();

    void UpdateRotationMatrix_FP();

    //===================================================
    //      INTERNAL FORCES
    //===================================================

    void UpdateInternalForces();

    void UpdateInternalForces_FP();

    void InitializeInternalForces();

    addouble GetDisplacement(int iNode, int iDim) {
        return U(6 * iNode + iDim);
    }

    inline void RegisterSolutionInput(int iRigid) {
        if (nRBE2 != 0 and iRigid == 1) {
            for (unsigned long i = 0; i < (nNode + nRBE2)*6; i++){
                AD::RegisterInput(U_lam(i)); }
        } 
        else {
            for (unsigned long i = 0; i < nNode * 6; i++){
                AD::RegisterInput(U(i));}
        }
    }

    inline void RegisterSolutionOutput(int iRigid) {
        if (nRBE2 != 0 and iRigid == 1) {
            for (unsigned long i = 0; i < (nNode + nRBE2)*6; i++){
                AD::RegisterOutput(U_lam(i)); }
        } 
        else {
            for (unsigned long i = 0; i < nNode * 6; i++){
                AD::RegisterOutput(U(i));}
        }
    }

    inline void ExtractSolutionAdjoint(int iRigid) {
        if (nRBE2 != 0 and iRigid == 1) {
           for (unsigned long i = 0; i < (nNode + nRBE2)*6; i++){
               U_lam_adj(i) = AD::GetDerivative(U_lam(i));}
        }
        else{
        for (unsigned long i = 0; i < nNode * 6; i++){
            U_adj(i) = AD::GetDerivative(U(i));}
        }
    }

    inline void SetSolutionAdjoint(int iRigid) {
        if (nRBE2 != 0 and iRigid == 1) {
        for (unsigned long i = 0; i < nNode * 6; i++){
            AD::SetDerivative(U_lam(i), AD::GetValue(U_lam_adj(i) + cross_term(i))); 
        }    
        for (unsigned long i = nNode * 6; i < (nNode + nRBE2)*6; i++){
            AD::SetDerivative(U_lam(i), AD::GetValue( U_lam_adj(i) ));
        }
        }
        else{
        for (unsigned long i = 0; i < nNode * 6; i++){
            AD::SetDerivative(U(i), AD::GetValue(U_adj(i) + cross_term(i)));}
        }
    }

    inline void StoreDisplacementAdjoint(int iNode, int iDim, passivedouble val_adj) {
        cross_term(6*iNode+iDim) = val_adj;
    }

    inline addouble GetCoordinates(int pos, int index) {return X(3*pos+index);}

    inline addouble GetInitialCoordinates(int pos, int index) {return X0(3*pos+index);}

    inline void SetLowVerbosity(void) { verbose = false; }
    inline void SetHighVerbosity(void) { verbose = true; }
};
