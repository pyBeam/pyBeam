/*
 * pyBeam, an open-source Beam Solver
 *
 * Copyright (C) 2019 by the authors
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
#include "../include/types.h"

#include "../include/input.h"
#include "../include/geometry.h"

class CRBE2
{
private:
public:
    
    int iRBE2;
    CNode* node_master;
    CNode* node_slave;
    
    int RBE2dofs = 12;   // RBE2 DOFs

    addouble l_rigid;                    // Rigid Length  
    VectorXdDiff axis_vector;            // axis vector master --> slave (non unitary)
    VectorXdDiff axis_vector0;           // axis vector master --> slave (non unitary)
    VectorXdDiff axis_vector_old;        // axis vector master --> slave (non unitary)

    VectorXi MasterDOFs  =  VectorXi::Zero(6); // Master DOFs
    VectorXi SlaveDOFs  =  VectorXi::Zero(6); // Slave DOFs  

    MatrixXdDiff Kinem_matrix;
    MatrixXdDiff Kinem_matrix0;
    MatrixXdDiff Kinem_matrix_old;    

    MatrixXdDiff MStrans;
    MatrixXdDiff MStrans0;
    MatrixXdDiff MStrans_old;

    VectorXdDiff g;       // constraint equations on LHS
    MatrixXdDiff J;       // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_0;     // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_1;     // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_2;     // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_3;     // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_4;     // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_5;     // Jacobian of the constraint equations on LHS

private:

public:

    // In the constructor we assign to the RBE2 its nodes 
    CRBE2(int RBE2_ID) ;

    ~CRBE2(void);

    void  Initializer(CNode* Node_mast, CNode* Node_slv);

    inline void SetNode_1( CNode* Node_mast) { node_master = Node_mast;}

    inline void SetNode_2( CNode* Node_slv) { node_slave = Node_slv;}

    void setGlobalDOFs();

    void setLength();

    void InitializeAxisVector();

    void InitializeKinemMatrix();

    void UpdateKinemMatirx();

    void EvalConstraintEquation( VectorXdDiff Um,VectorXdDiff Us);

    void EvalJacobian( VectorXdDiff Um,VectorXdDiff Us);

    void EvalHessian( VectorXdDiff Um,VectorXdDiff Us);

};
