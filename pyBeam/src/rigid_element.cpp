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


#include "../include/rigid_element.h"
#include <iostream>

CRBE2::CRBE2(int RBE2_ID) {
    iRBE2 = RBE2_ID;
};

void CRBE2::Initializer(CNode* Node_mast, CNode* Node_slv) {

    // Associate the nodes object
    SetNodeMaster(Node_mast);
    SetNodeSlave(Node_slv);
    // Calculate element DOFs
    setGlobalDOFs();
    // set rigid element length
    setLength();
    // Initialize dimensional axis vector
    InitializeAxisVector();


    InitializeJacobian();
    InitializeHessian();

}



void CRBE2::setGlobalDOFs() {
    
    int nodeIndexM = node_master-> GeID();
    int nodeIndexS = node_slave-> GeID();
    int i = 0;
    for (int i = 1; i <= 6; i++) {
        // This coefficients will be used to asseble full_to_red and red_to_full vectors which have 0 expressing no relation so it is decided that DOFs will be from 1 to 6
        MasterDOFs(i - 1) = (nodeIndexM - 1)*6 + i; // -1;
        SlaveDOFs(i - 1) = (nodeIndexS - 1)*6 + i; // -1;
    };
    
};

void CRBE2::setLength() {
    addouble a = node_slave->GetCoordinate0(0) - node_master->GetCoordinate0(0);
    addouble b = node_slave->GetCoordinate0(1) - node_master->GetCoordinate0(1);
    addouble c = node_slave->GetCoordinate0(2) - node_master->GetCoordinate0(2);
    addouble intermediate = pow(a, 2) + pow(b, 2) + pow(c, 2);
    l_rigid = sqrt(intermediate);
    
};

void CRBE2::InitializeAxisVector() {
    
    addouble a = node_slave->GetCoordinate0(0) - node_master->GetCoordinate0(0);
    addouble b = node_slave->GetCoordinate0(1) - node_master->GetCoordinate0(1);
    addouble c = node_slave->GetCoordinate0(2) - node_master->GetCoordinate0(2);
    
    axis_vector0 = VectorXdDiff::Zero(3);
    axis_vector0(0) = a;
    axis_vector0(1) = b;
    axis_vector0(2) = c;
    
};

void CRBE2::InitializeJacobian( ) {
    
    MStrans = MatrixXdDiff::Zero(3, 3);
    
    
    // L matrix in theory
    /*        MStrans = [        0         axis_vector(3)   -axis_vector(2);
     -axis_vector(3)        0           axis_vector(1);
     axis_vector(2) -axis_vector(1)        0            ];*/
    
    MStrans << 0, axis_vector0(3 - 1), -axis_vector0(2 - 1),
            -axis_vector0(3 - 1), 0, axis_vector0(1 - 1),
            axis_vector0(2 - 1), -axis_vector0(1 - 1), 0;
    
    G.block(0,0, 3 , 3) =  -MatrixXdDiff::Identity(3,3);
    G.block(3,3, 3 , 3) =  -MatrixXdDiff::Identity(3,3);
    
    G.block(0,6, 3 , 3) =  MatrixXdDiff::Identity(3,3);
    G.block(3,9, 3 , 3) =  MatrixXdDiff::Identity(3,3);    
    
    G.block(0,3, 3 , 3) =  -MStrans;  
    

    //std::cout << "G " <<  " = \n" << G  <<std::endl;    
};

void CRBE2::InitializeHessian() {
    
    // H_0
    MatrixXdDiff H_0_red = MatrixXdDiff::Zero(3,3);
    H_0_red << 0, axis_vector0(2 - 1), axis_vector0(3 - 1),
            axis_vector0(2 - 1), -2*axis_vector0(1 - 1),0,
            axis_vector0(3 - 1), 0, -2*axis_vector0(1 - 1);
    
    H_0.block(4-1,4-1,3,3) = H_0_red/2;
    
    // H_1
    MatrixXdDiff H_1_red = MatrixXdDiff::Zero(3,3);
    H_1_red << -2*axis_vector0(2 - 1), axis_vector0(1 - 1), 0,
            axis_vector0(1 - 1), 0,axis_vector0(3 - 1),
            0, axis_vector0(3 - 1), -2*axis_vector0(2 - 1);
    
    H_1.block(4-1,4-1,3,3) = H_1_red/2;    
    
    // H_2
    MatrixXdDiff H_2_red = MatrixXdDiff::Zero(3,3);
    H_2_red << -2*axis_vector0(3 - 1), 0, axis_vector0(1 - 1),
            0, -2*axis_vector0(3 - 1),axis_vector0(2 - 1),
            axis_vector0(1 - 1), axis_vector0(2 - 1), 0;
    
    H_2.block(4-1,4-1,3,3) = H_2_red/2;     
    /*
    std::cout << "H_0 " <<  " = \n" << H_0  <<std::endl;  
    std::cout << "H_1 " <<  " = \n" << H_1  <<std::endl;
    std::cout << "H_2 " <<  " = \n" << H_2  <<std::endl;
    */
};

void CRBE2::EvalConstraintEquation(VectorXdDiff Um, VectorXdDiff Us) {
    
    Vector3dDiff U_M_rot= Vector3dDiff::Zero();
    Matrix3dDiff R= Matrix3dDiff::Zero();
    VectorXdDiff Us_p = VectorXdDiff::Zero(3);
    g = VectorXdDiff::Zero(6);
    U_M_rot = Um.segment(4-1,3);      // Rotation at triad A
    //std::cout << "U_M_rot " <<  " = \n" << U_M_rot  <<std::endl;
    PseudoToRot(U_M_rot, R);
        
    R = R - MatrixXdDiff::Identity(3,3);
    
    g = Us - Um;
    g.segment(0,3) += -R*axis_vector0;
    //Us_p = Um.segment(1-1,3) + R*axis_vector0;
    //std::cout << "g " <<  " = \n" << g  <<std::endl;
    //std::cout << "R " <<  " = \n" << R  <<std::endl;
    //std::cout << "axis_vector0 " <<  " = \n" << axis_vector0  <<std::endl;
    //std::cout << "Us_constraint " <<  " = \n" << Us_p  <<std::endl;
    std::ofstream myfile;
    myfile.open ("g.pyBeam", fstream::in | fstream::out | fstream::app);
    myfile << "g = "; myfile << g ; myfile << "\n ";
    myfile.close();
};

CRBE2::~CRBE2(void) {

};