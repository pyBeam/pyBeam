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
    MatrixXdDiff
    int DOF;                    // In space, 6
    
    CRBE2 **RBE2;               // Pointer to the first RBE2 element
    CElement **element;         // Pointer to the first finite element
    CNode **node;               // Pointer to the first finite element

    VectorXdDiff cross_term;    // Store the displacement vector
    
    MatrixXdDiff Ksys;
    MatrixXdDiff Ksys_red;      // [relative to masters in case of RBE2]
    MatrixXdDiff K_penal;       // penalty matrix for rigid elements
    VectorXdDiff V_penal;       // penalty vector for rigid elements
    
    MatrixXdDiff KRBE;          // Kinematic constraint matrix due to the RBE2 elements   [totalDOFs, BossDOFs]
    MatrixXdDiff KRBE_ext;      // Kinematic constraint matrix due to the RBE2 elements   [totalDOFs, BossDOFs]
    
    MatrixXdDiff  Constr_matrix;// COnstraint matrix [ NODE_ID DOF_ID ]
    
    VectorXdDiff U;             // Displacement array
    VectorXdDiff dU;            // Displacement array (increment)
    VectorXdDiff dU_red;        // Displacement array (increment) [relative to masters in case of RBE2]
    VectorXdDiff X;             // Position of the fem nodes in global coordinate system
    VectorXdDiff X0;            // Position of the fem nodes in global coordinate system

    VectorXdDiff U_adj;         // Adjoint of the displacement array (cumulative)
    
    VectorXdDiff Fpenal;        // Array of internal forces
    VectorXdDiff Fint;          // Array of internal forces
    VectorXdDiff Fext;          // Array of External Forces
    VectorXdDiff Residual;      // Array of Unbalanced Forces
    VectorXdDiff Residual_red;  // Array of Unbalanced Forces   [relative to masters in case of RBE2]
    
    VectorXdDiff Fnom;          // Array of nominal forces
    
    addouble YoungModulus;
    
    CStructure(CInput *input, CElement **container_element, CNode **container_node);
    
    ~CStructure();
    
    /*##############################################################
     *
     *         External Forces, Residual, Intenral Forces
     *
     *###############################################################*/
    
    inline void ReadForces(int nTotalDOF, addouble *loadVector) {
        for (int iLoad = 0; iLoad < nTotalDOF; iLoad++){Fnom(iLoad) = loadVector[iLoad];}
    }

    inline void SetDimensionalYoungModulus(addouble val_E){ YoungModulus = val_E; }

    // External forces are normalized by the Young Modulus
    inline void UpdateExtForces(addouble lambda){ Fext = lambda* Fnom / YoungModulus; }

    void EvalResidual(unsigned short rigid);

    void EvalPenaltyForces(addouble penalty);

    //===================================================
    //      Assembly RBE2 rigid constraint matrix
    //===================================================

    inline void AddRBE2(CInput *input, CRBE2** container_RBE2) {nRBE2 = input->Get_nRBE2(); RBE2 = container_RBE2;}

    void AssemblyRigidConstr();

    void AssemblyRigidPenalty(addouble penalty);

    void UpdateRigidConstr(int iIter);

    //===================================================
    //      Assembly System Stiffness Matrix
    //===================================================

    void AssemblyTang(int iIter);

    void EvalSensRot();    // Evaluate the sensitivity of Rotation Matrix - need for Jacobian

    //void EvalSensRotFiniteDifferences();    // Evaluate the sensitivity of Rotation Matrix - need for Jacobian

    //===================================================
    //      Solve linear static system
    //===================================================
    // Assembles LHS and RHS and solves the linear static problem

    void SolveLinearStaticSystem(int iIter);

    void SolveLinearStaticSystem_RBE2(int iIter);

    void SolveLinearStaticSystem_RBE2_penalty(int iIter);

    //===================================================
    //      Update Coordinates
    //===================================================
    /* This member function upates the coordinates (expressed in global reference system) of
     * the finite element nodes.
     * This is necessary for "booking" the position, as the compatiblity and the equations
     * are based on the displacements.
     */

    void UpdateCoord();

    void RestartCoord();

    void UpdateAxvector_RBE2();

    void UpdateCoord_RBE2(int iIter);

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
        return U(6*iNode+iDim);
    }

    inline void RegisterSolutionInput(void) {
        for (unsigned long i = 0; i < nNode * 6; i++)
            AD::RegisterInput(U(i));
    }

    inline void RegisterSolutionOutput(void) {
        for (unsigned long i = 0; i < nNode * 6; i++)
            AD::RegisterOutput(U(i));
    }

    inline void ExtractSolutionAdjoint(void) {
        for (unsigned long i = 0; i < nNode * 6; i++){
            U_adj(i) = AD::GetDerivative(U(i));
        }
    }

    inline void SetSolutionAdjoint(void) {
        for (unsigned long i = 0; i < nNode * 6; i++)
            AD::SetDerivative(U(i), AD::GetValue(U_adj(i) + cross_term(i)));
    }

    inline void StoreDisplacementAdjoint(int iNode, int iDim, passivedouble val_adj) {
        cross_term(6*iNode+iDim) = val_adj;
    }

    inline addouble GetCoordinates(int pos, int index) {return X(3*pos+index);}

    inline addouble GetInitialCoordinates(int pos, int index) {return X0(3*pos+index);}

    inline void SetLowVerbosity(void) { verbose = false; }
    inline void SetHighVerbosity(void) { verbose = true; }
};
