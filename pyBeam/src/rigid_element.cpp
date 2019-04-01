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


#include "../include/rigid_element.h"
#include <iostream>


CRBE2::CRBE2(int RBE2_ID) { iRBE2 = RBE2_ID; };

void CRBE2::setGlobalDOFs(){
    
    int nodeIndexM =  node_master-> GeID();
    int nodeIndexS =  node_slave-> GeID();
    int i = 0;
    for (int i=1; i<=6; i++) {
        // This coefficients will be used to asseble full_to_red and red_to_full vectors which have 0 expressing no relation so it is decided that DOFs will be from 1 to 6 
        MasterDOFs(i-1) = (nodeIndexM-1)*6 + i;// -1;      
        SlaveDOFs(i-1) = (nodeIndexS-1)*6 + i;// -1;       
    };
    
}; 

void CRBE2::setLength() {
    addouble a = node_slave->GetCoordinate0(0) - node_master->GetCoordinate0(0);
    addouble b = node_slave->GetCoordinate0(1) - node_master->GetCoordinate0(1);
    addouble c = node_slave->GetCoordinate0(2) - node_master->GetCoordinate0(2);
    addouble intermediate = pow(a ,2) + pow(b,2) + pow( c ,2) ;
    l_rigid =  sqrt(intermediate );
    
};   

void CRBE2::InitializeAxisVector() {
    
    addouble a = node_slave->GetCoordinate0(0) - node_master->GetCoordinate0(0);
    addouble b = node_slave->GetCoordinate0(1) - node_master->GetCoordinate0(1);
    addouble c = node_slave->GetCoordinate0(2) - node_master->GetCoordinate0(2);
    
    axis_vector = VectorXdDiff::Zero(3);
    axis_vector(0) = a; axis_vector(1) = b; axis_vector(2) = c;
    
    axis_vector0 = VectorXdDiff::Zero(3);
    axis_vector_old = axis_vector0;
    axis_vector = axis_vector0;   
}
    
    void CRBE2::InitializeKinemMatrix() { 
        
        Kinem_matrix = MatrixXdDiff::Identity(6,6);  
        MatrixXdDiff MStrans = MatrixXdDiff::Zero(3,3);  
        
        /*        Kprim = [        0         axis_vector(3)   -axis_vector(2);
         -axis_vector(3)        0           axis_vector(1);
         axis_vector(2) -axis_vector(1)        0            ];*/
        
        MStrans <<   0,                  axis_vector0(3 -1),   -axis_vector0(2 -1),
                -axis_vector0(3 -1),        0,              axis_vector0(1 -1),
                axis_vector0(2 -1),    -axis_vector0(1 -1),        0      ; 
        
        Kinem_matrix.block(1 -1, 4 -1, 3, 3) = MStrans;
        Kinem_matrix0 = Kinem_matrix;
        Kinem_matrix_old = Kinem_matrix;
        
    }   
    
    void CRBE2::UpdateKinemMatirx() {
        
        Kinem_matrix_old = Kinem_matrix;        
        
        Kinem_matrix = MatrixXdDiff::Identity(6,6);        
        MatrixXdDiff MStrans = MatrixXdDiff::Zero(3,3); 
        
        MStrans <<   0,                  axis_vector(3 -1),   -axis_vector(2 -1),
                -axis_vector(3 -1),        0,              axis_vector(1 -1),
                axis_vector(2 -1),    -axis_vector(1 -1),        0      ; 
        
        Kinem_matrix.block(1 -1, 4 -1, 3, 3) = MStrans;        
        
    } 
    
    void CRBE2::Initializer(CNode* Node_mast, CNode* Node_slv){
        
        // Associate the nodes object   
        SetNode_1( Node_mast) ;
        SetNode_2( Node_slv);   
        // Calculate element DOFs
        setGlobalDOFs();    
        // set rigid element length
        setLength();
        // Initialize dimensional axis vector
        InitializeAxisVector();
        
        
        InitializeKinemMatrix();
        
        
    }
    
    CRBE2::~CRBE2(void) {};