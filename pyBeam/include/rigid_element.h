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
    VectorXdDiff axis_vector;             // axis vector master --> slave (non unitary)
    
    VectorXdDiff GlobalDOFs  =  VectorXdDiff::Zero(12); // Global DOFs 
    MatrixXdDiff Kinem_matrix;
    
    // Menber functions
private:
    
    
public:
    
    
    // In the constructor we assign to the RBE2 its nodes 
    CRBE2(int RBE2_ID) ; 
    
    ~CRBE2(void);

    void  Initializer(CNode* Node_mast, CNode* Node_slv);     
    
    inline void SetNode_1( CNode* Node_mast) { node_master = Node_mast; };
    
    inline void SetNode_2( CNode* Node_slv) { node_slave = Node_slv;};  
    
    void setGlobalDOFs();    
    
    void setLength();
    
};    