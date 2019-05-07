/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rocco Bombardieri, Rauno Cavallaro 
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

#ifdef DEBG
#include <fstream>
#endif

class CStructure
{
    
private:
    
public:
    
    
    int nNode;               // Number of structural nodes  
    int nfem;                // number of finite elements   // has to be assigned in the constructor
    int nRBE2;                // number of finite elements
    
    int DOF;                  // In space, 6
    int FollFlag;             // Flag for Follower forces (1)
    
    CRBE2 **RBE2;      // Pointer to the first RBE2 element        
    CElement **fem;      // Pointer to the first finite element
    CNode **node;        // Pointer to the first finite element
    
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
    
    VectorXdDiff Fpenal;        // Array of internal forces
    VectorXdDiff Fint;        // Array of internal forces
    VectorXdDiff Fext;        // Array of External Forces
    VectorXdDiff Residual;    // Array of Unbalanced Forces
    VectorXdDiff Residual_red;    // Array of Unbalanced Forces   [relative to masters in case of RBE2]      
    
    Vector3dDiff Ftip;      // (vector read from the input file) - needed as basis to update the Fext
    VectorXdDiff Fnom;        // Array of nominal forces
    
    CStructure(CInput *input, CElement **element, CNode **container_node);
    
    ~CStructure();
    
    /*##############################################################
     *
     *         External Forces, Residual, Intenral Forces
     *
     *###############################################################*/
    
    void ReadForces(int nTotalDOF, addouble *loadVector);
    
    void UpdateExtForces(addouble  );
    
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
    
    void EvalSensRot();  // Evaluate the sensitivity of Rotation Matrix - need for Jacobian
    
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
    
    void UpdateAxvector_RBE2(); 

    void UpdateCoord_RBE2(int iIter);        
    
    void InitialCoord();
    
    void UpdateLength();
    
    void UpdateRotationMatrix();
    
    //===================================================
    //     TOOLS: WRITING COORDINATES
    //===================================================
    
    void EchoCoord();
    
    void EchoDisp();
    
    void EchoRes();
    
    void EchoFext();
    
    void EchoMatrixK();
    
    //===================================================
    //      INTERNAL FORCES
    //===================================================
        
    void UpdateInternalForces();
    
    addouble GetDisplacement(int pos, int index) {
        addouble disp;
        disp = X(3*pos+index) - X0(3*pos+index);
        return disp;
    };
    
    addouble GetCoordinates(int pos, int index) {return X(3*pos+index);};
    
    addouble GetInitialCoordinates(int pos, int index) {return X0(3*pos+index);};
    
};
