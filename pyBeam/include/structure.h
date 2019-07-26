/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rocco Bombardieri, Rauno Cavallaro 
 * 
 * File developers: Rocco Bombardieri (Carlos III University Madrid)
 *                  Rauno Cavallaro (Carlos III University Madrid)
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
    
    
    int nNode;               // Number of structural nodes  
    int nfem;                // number of finite elements   // has to be assigned in the constructor
    int nRBE2;                // number of finite elements
    
    int DOF;                  // In space, 6
    
    CRBE2 **RBE2;        // Pointer to the first RBE2 element
    CElement **element;  // Pointer to the first finite element
    CNode **node;        // Pointer to the first finite element

    addouble **disp;     // Store the displacement vector
    passivedouble **disp_adj;     // Store the displacement vector
    
    MatrixXdDiff M;      // Recall in Eigen X stays for dynamic, d for addouble:  (nfem+1)*6  X   (nfem+1)*6
    MatrixXdDiff Ksys;
    MatrixXdDiff Ksys_red; // [relative to masters in case of RBE2]       
    MatrixXdDiff K_penal;  // penalty matrix for rigid elements
    VectorXdDiff V_penal;  // penalty vector for rigid elements
    
    MatrixXdDiff KRBE;  // Kinematic constraint matrix due to the RBE2 elements   [totalDOFs, BossDOFs]   
    MatrixXdDiff KRBE_ext;  // Kinematic constraint matrix due to the RBE2 elements   [totalDOFs, BossDOFs]
    
    MatrixXdDiff  Constr_matrix;    // COnstraint matrix [ NODE_ID DOF_ID ]
    
    VectorXdDiff U;             // Displacement array (cumulative)
    VectorXdDiff dU;           // Displacement array (iterative)
    VectorXdDiff dU_red;           // Displacement array (iterative) [relative to masters in case of RBE2]        
    VectorXdDiff X;            // Position of the fem nodes in global coordinate system
    VectorXdDiff X0;            // Position of the fem nodes in global coordinate system

    VectorXdDiff U_adj;         // Adjoint of the displacement array (cumulative)
    
    VectorXdDiff Fpenal;        // Array of internal forces
    VectorXdDiff Fint;        // Array of internal forces
    VectorXdDiff Fext;        // Array of External Forces
    VectorXdDiff Residual;    // Array of Unbalanced Forces
    VectorXdDiff Residual_red;    // Array of Unbalanced Forces   [relative to masters in case of RBE2]      
    
    VectorXdDiff Fnom;        // Array of nominal forces
    
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
    
    inline void UpdateExtForces(addouble lambda){ Fext = lambda* Fnom; }

    void EvalResidual(unsigned short rigid);

    void EvalPenaltyForces(addouble penalty); 
    
    //===================================================
    //      Assembly RBE2 rigid constraint matrix
    //===================================================        
    
    void AddRBE2(CInput *input, CRBE2** container_RBE2) {nRBE2 = input->Get_nRBE2(); RBE2 = container_RBE2;};
    
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
    /* This member function upates the coordinates (expressed in global reference system) of the finite element nodes.
     * This is necessary for "booking" the position, as the compatiblity and the equations are based on the dispalcements.
     */
    
    void UpdateCoord();
    
    void RestartCoord();
    
    void UpdateAxvector_RBE2(); 

    void UpdateCoord_RBE2(int iIter);        
    
    void InitialCoord();
    
    void UpdateLength();
    
    void UpdateRotationMatrix();
    
    void UpdateRotationMatrix_FP();    
    
    //===================================================
    //      INTERNAL FORCES
    //===================================================
        
    void UpdateInternalForces();

    void UpdateInternalForces_FP();    
    
    void InitializeInternalForces();
    
    addouble GetDisplacement(int pos, int index) {
        return disp[pos][index];
    };

    inline void SetDisplacement(int iNode, int iDim) {
        disp[iNode][iDim] = X(3*iNode+iDim) - X0(3*iNode+iDim);
    };

    inline void RegisterDisplacement(int iNode, int iDim) {
        AD::RegisterOutput(disp[iNode][iDim]);
    };

    inline void SetSolutionDependencies(void) {
        for (unsigned long i = 0; i < nNode; i++){
           X(i*3) = U(i*6) ;
           X(i*3+1) = U(i*6+1);
           X(i*3+2) = U(i*6+2);
        }
    };

    inline void RegisterSolutionInput(void) {
        for (unsigned long i = 0; i < nNode * 6; i++)
          AD::RegisterInput(U(i));
    };

    inline void RegisterSolutionOutput(void) {
        for (unsigned long i = 0; i < nNode * 6; i++)
          AD::RegisterOutput(U(i));
    };

    inline void ExtractSolutionAdjoint(void) {
        for (unsigned long i = 0; i < nNode * 6; i++){
          U_adj(i) = AD::GetDerivative(U(i));
        }
    };

    inline void SetSolutionAdjoint(void) {
        for (unsigned long i = 0; i < nNode * 6; i++)
          AD::SetDerivative(U(i), AD::GetValue(U_adj(i)));
    };

    inline void StoreDisplacementAdjoint(int iNode, int iDim, passivedouble val_adj) {
        disp_adj[iNode][iDim] = val_adj;
    };

    inline void SetDisplacementAdjoint(int iNode, int iDim) {
        AD::SetDerivative(disp[iNode][iDim], AD::GetValue(disp_adj[iNode][iDim]));
    };
    
    addouble GetCoordinates(int pos, int index) {return X(3*pos+index);};
    
    addouble GetInitialCoordinates(int pos, int index) {return X0(3*pos+index);};
    
};
