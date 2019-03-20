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
    
    int nodeIndexA =  nodeA-> GeID();
    int nodeIndexB =  nodeB-> GeID();
    int i = 0;
    for (int i=1; i<=6; i++) {
        GlobalDOFs(i-1) = (nodeIndexA-1)*6 + i -1;      
        GlobalDOFs(6+i-1) = (nodeIndexB-1)*6 + i -1;       
    };
    
}; 

void CElement::setLength() {
    addouble a = node_slave->GetCoordinate0(0) - node_master->GetCoordinate0(0);
    addouble b = node_slave->GetCoordinate0(1) - node_master->GetCoordinate0(1);
    addouble c = node_slave->GetCoordinate0(2) - node_master->GetCoordinate0(2);
    addouble intermediate = pow(a ,2) + pow(b,2) + pow( c ,2) ;
    l_rigid =  sqrt(intermediate );
    
    axisvector = VectorXdDiff::Zero(6);
    axisvector(0) = a; axisvector(1) = b; axisvector(2) = c;
};   

void CElement::Initializer(CNode* Node_mast, CNode* Node_slv){
    
    // Associate the nodes object   
    SetNode_1( Node_mast) ;
    SetNode_2( Node_slv);  
    // Calculate element DOFs
    setGlobalDOFs();    
    // set rigid element length
    setLength();
    
    Kprim = MatrixXdDiff::Zero(6,6);    
    
};