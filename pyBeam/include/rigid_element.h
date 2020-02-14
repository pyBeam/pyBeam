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
#include "../include/rotations.h"
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
    VectorXdDiff ax;           // axis vector master --> slave (non unitary)

    VectorXi MasterDOFs  =  VectorXi::Zero(6); // Master DOFs
    VectorXi SlaveDOFs  =  VectorXi::Zero(6); // Slave DOFs  

    MatrixXdDiff MStrans;

    VectorXdDiff g = VectorXdDiff::Zero(6);       // constraint equations on LHS
    MatrixXdDiff G = MatrixXdDiff::Zero(6,12);       // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_0 = MatrixXdDiff::Zero(12,12);     // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_1 = MatrixXdDiff::Zero(12,12);     // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_2 = MatrixXdDiff::Zero(12,12);     // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_3 = MatrixXdDiff::Zero(12,12);     // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_4 = MatrixXdDiff::Zero(12,12);     // Jacobian of the constraint equations on LHS
    MatrixXdDiff H_5 = MatrixXdDiff::Zero(12,12);     // Jacobian of the constraint equations on LHS
    
private:

public:

    // In the constructor we assign the id
    CRBE2(int RBE2_ID) ;

    ~CRBE2(void);

    void  Initializer(CNode* Node_mast, CNode* Node_slv);

    inline void SetNodeMaster( CNode* Node_mast) { node_master = Node_mast;};

    inline void SetNodeSlave( CNode* Node_slv) { node_slave = Node_slv;};

    void setGlobalDOFs();

    void setLength();

    void InitializeAxisVector();
    
    void EvalConstraintEquation( VectorXdDiff Um,VectorXdDiff Us);
      
    void InitializeJacobian( );
    
    void EvalJacobian(VectorXdDiff Um );

    void InitializeHessian();    
    
    };
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
