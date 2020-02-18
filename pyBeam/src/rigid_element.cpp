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
        // This coefficients will be used to assemble full_to_red and red_to_full vectors which have 0 expressing no relation so it is decided that DOFs will be from 1 to 6
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
    
    ax = VectorXdDiff::Zero(3);
    ax(0) = a;
    ax(1) = b;
    ax(2) = c;
    
};

void CRBE2::InitializeJacobian( ) {
    
    MStrans = MatrixXdDiff::Zero(3, 3);
    
    
    // L matrix in theory
    /*        MStrans = [        0         axis_vector(3)   -axis_vector(2);
     -axis_vector(3)        0           axis_vector(1);
     axis_vector(2) -axis_vector(1)        0            ];*/
    
    MStrans << 0, ax(3 - 1), -ax(2 - 1),
            -ax(3 - 1), 0, ax(1 - 1),
            ax(2 - 1), -ax(1 - 1), 0;
    
    G.block(0,0, 3 , 3) =  -MatrixXdDiff::Identity(3,3);
    G.block(3,3, 3 , 3) =  -MatrixXdDiff::Identity(3,3);
    
    G.block(0,6, 3 , 3) =  MatrixXdDiff::Identity(3,3);
    G.block(3,9, 3 , 3) =  MatrixXdDiff::Identity(3,3);    
    
    G.block(0,3, 3 , 3) =  -MStrans;  
    

};

void CRBE2::EvalJacobian_2(VectorXdDiff Um ) {
    
    
    // Jacobian depends on the rotations of the master node
    /*
     * G = | -I  -d (R(U_m_rot)*axis_vector)/dU_m_rot  I   0  |
     *     |  0                I                       0  -I  |
     */
    Vector3dDiff U_M_rot= Vector3dDiff::Zero();
        
    // I call L = -d (R(U_m_rot)*axis_vector)/dU_m_rot 3x3
    MatrixXdDiff L = MatrixXdDiff::Zero(3, 3);
    
    U_M_rot = Um.segment(4-1,3);      // Rotations of master

    // Let's define the ingredients
    
    
    // theta (module ofU_M_rot )
    addouble th = sqrt(  pow(U_M_rot(0),2) + pow(U_M_rot(1),2) + pow(U_M_rot(2),2) );
    if (th < 1.0e-16 ) {th = 1.0e-16; };  //to avoid singularities
    // some relative scalars
    // sin th/th
    addouble sth = sin(th)/th;
    //(1-cos th )/ th^2
    addouble ccsth = (1 - cos(th))/ pow(th,2);
    
    // First order Rodriguez rotation matrix
    MatrixXdDiff K = MatrixXdDiff::Zero(3, 3);
    K << 0            , -U_M_rot(3-1)  , U_M_rot(2-1),
            U_M_rot(3-1) ,      0         , -U_M_rot(1-1),
            -U_M_rot(2-1) ,  U_M_rot(1-1)  ,     0        ;    
    
    // d(K*axis vector)/ d(U_M_rot)
    MatrixXdDiff Kd_prime = MatrixXdDiff::Zero(3, 3);
    Kd_prime << 0, ax(3 - 1), -ax(2 - 1),
            -ax(3 - 1), 0, ax(1 - 1),
            ax(2 - 1), -ax(1 - 1), 0;


   // K*axis_vector
    Vector3dDiff KA = Vector3dDiff::Zero();
    KA = K*ax;
    
    // K*K*axis_vector
    Vector3dDiff KKA = Vector3dDiff::Zero();
    KKA = K*(K*ax);
    
    // d(K*K*axis vector)/ d(U_M_rot)  
    MatrixXdDiff D = MatrixXdDiff::Zero(3, 3);  
    // Careful here
    D <<  U_M_rot(2 -1)*ax(2 - 1) + U_M_rot(3 -1)*ax(3 - 1) , -2*U_M_rot(2 -1)*ax(1 - 1) + U_M_rot(1 -1)*ax(2 - 1), -2*U_M_rot(3 -1)*ax(1 - 1) + U_M_rot(1 -1)*ax(3 - 1),
          U_M_rot(2 -1)*ax(1 - 1) -2*U_M_rot(1 -1)*ax(2 - 1), U_M_rot(1 -1)*ax(1 - 1) + U_M_rot(3 -1)*ax(3 - 1)   , -2*U_M_rot(3 -1)*ax(2 - 1) + U_M_rot(2 -1)*ax(3 - 1),
          U_M_rot(3 -1)*ax(1 - 1) -2*U_M_rot(1 -1)*ax(3 - 1), U_M_rot(3 -1)*ax(2 - 1) -2*U_M_rot(2 -1)*ax(3 - 1)  , U_M_rot(1 -1)*ax(1 - 1) + U_M_rot(2 -1)*ax(2 - 1); 


    // L assembly. So L has 4 contributions:
    // 1
    L = (cos(th)*th - sin(th))/pow(th,3) * KA* U_M_rot.transpose();
    // 2
    L = L + sth*Kd_prime;
    //3
    L = L + (th*sin(th)+2*cos(th) - 2)/pow(th,4) *KKA* U_M_rot.transpose();
    // 4
    L = L + ccsth*D;
    
    G = MatrixXdDiff::Zero(6,12);       // Jacobian of the constraint equations on LHS    
    
    G.block(0,0, 3 , 3) =  -MatrixXdDiff::Identity(3,3);
    G.block(3,3, 3 , 3) =  -MatrixXdDiff::Identity(3,3);
    
    G.block(0,6, 3 , 3) =  MatrixXdDiff::Identity(3,3);
    G.block(3,9, 3 , 3) =  MatrixXdDiff::Identity(3,3); 


    G.block(0,3, 3 , 3) =  -L;  
    
    //std::cout << "G " <<  " = \n" << G  <<std::endl;    
    //std::cout << "G^T*G " <<  " = \n" << G.transpose()*G  <<std::endl;     
    
}
void CRBE2::EvalJacobian(VectorXdDiff Um ) {
    //std::cout << "\nSymbolic\n" << std::endl;
    // G
    G = MatrixXdDiff::Zero(6,12);       // Jacobian of the constraint equations on LHS   
    addouble Umrx = Um(4-1); addouble Umry = Um(5-1); addouble Umrz = Um(6-1);
    //U_M_rot = Um.segment(4-1,3);      // Rotation at triad A
    
    addouble theta = sqrt(pow(Umrx,2) + pow(Umry,2) + pow(Umrz,2)); 
    addouble theta_sq = ( pow(Umrx,2) + pow(Umry,2) + pow(Umrz,2) );
    if (theta < 1.0e-16 ) {theta = 1.0e-16; };  //to avoid singularities 
    if (theta_sq < 1.0e-16 ) {theta_sq = 1.0e-16; };  //to avoid singularities     

       G(0,0) = -1.0;
       G(0,3) = ((ax(1)*Umry+ax(2)*Umrz)*(cos(theta)-1.0))/(theta_sq)+(Umrx*cos(theta)*(ax(1)*Umrz-ax(2)*Umry))/(theta_sq)-Umrx*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)-Umrx*sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umrx*(cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)*2.0;
       G(0,4) = -((ax(0)*Umry*2.0-ax(1)*Umrx)*(cos(theta)-1.0))/(theta_sq)-ax(2)*sin(theta)*1.0/theta+(Umry*cos(theta)*(ax(1)*Umrz-ax(2)*Umry))/(theta_sq)-Umry*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)-Umry*sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umry*(cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)*2.0;
       G(0,5) = -((ax(0)*Umrz*2.0-ax(2)*Umrx)*(cos(theta)-1.0))/(theta_sq)+ax(1)*sin(theta)*1.0/theta+(Umrz*cos(theta)*(ax(1)*Umrz-ax(2)*Umry))/(theta_sq)-Umrz*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)-Umrz*sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umrz*(cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)*2.0;
       G(0,6) = 1.0;
       G(1,1) = -1.0;
       G(1,3) = ((ax(0)*Umry-ax(1)*Umrx*2.0)*(cos(theta)-1.0))/(theta_sq)+ax(2)*sin(theta)*1.0/theta-(Umrx*cos(theta)*(ax(0)*Umrz-ax(2)*Umrx))/(theta_sq)+Umrx*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-Umrx*sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umrx*(cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*2.0;
       G(1,4) = ((ax(0)*Umrx+ax(2)*Umrz)*(cos(theta)-1.0))/(theta_sq)-(Umry*cos(theta)*(ax(0)*Umrz-ax(2)*Umrx))/(theta_sq)+Umry*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-Umry*sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umry*(cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*2.0;
       G(1,5) = -((ax(1)*Umrz*2.0-ax(2)*Umry)*(cos(theta)-1.0))/(theta_sq)-ax(0)*sin(theta)*1.0/theta-(Umrz*cos(theta)*(ax(0)*Umrz-ax(2)*Umrx))/(theta_sq)+Umrz*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-Umrz*sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umrz*(cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*2.0;
       G(1,7) = 1.0;
       G(2,2) = -1.0;
       G(2,3) = ((ax(0)*Umrz-ax(2)*Umrx*2.0)*(cos(theta)-1.0))/(theta_sq)-ax(1)*sin(theta)*1.0/theta+(Umrx*cos(theta)*(ax(0)*Umry-ax(1)*Umrx))/(theta_sq)-Umrx*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-Umrx*sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umrx*(cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*2.0;
       G(2,4) = ((ax(1)*Umrz-ax(2)*Umry*2.0)*(cos(theta)-1.0))/(theta_sq)+ax(0)*sin(theta)*1.0/theta+(Umry*cos(theta)*(ax(0)*Umry-ax(1)*Umrx))/(theta_sq)-Umry*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-Umry*sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umry*(cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*2.0;
       G(2,5) = ((ax(0)*Umrx+ax(1)*Umry)*(cos(theta)-1.0))/(theta_sq)+(Umrz*cos(theta)*(ax(0)*Umry-ax(1)*Umrx))/(theta_sq)-Umrz*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-Umrz*sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umrz*(cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*2.0;
       G(2,8) = 1.0; 
       
       G.block(3,3, 3 , 3) =  -MatrixXdDiff::Identity(3,3);
       G.block(3,9, 3 , 3) =  MatrixXdDiff::Identity(3,3); 
       
    //std::cout << "G " <<  " = \n" << G  <<std::endl;    
    //std::cout << "G^T*G " <<  " = \n" << G.transpose()*G  <<std::endl;        
}
void CRBE2::EvalHessian(VectorXdDiff Um){
 
    // H_0
    H_0 = MatrixXdDiff::Zero(12,12);
    H_1 = MatrixXdDiff::Zero(12,12);
    H_2 = MatrixXdDiff::Zero(12,12);
    addouble Umrx = Um(4-1); addouble Umry = Um(5-1); addouble Umrz = Um(6-1);
    //U_M_rot = Um.segment(4-1,3);      // Rotation at triad A
    
    addouble theta = sqrt(pow(Umrx,2) + pow(Umry,2) + pow(Umrz,2)); 
    addouble theta_sq = ( pow(Umrx,2) + pow(Umry,2) + pow(Umrz,2) );
    if (theta < 1.0e-16 ) {theta = 1.0e-16; };  //to avoid singularities 
    if (theta_sq < 1.0e-16 ) {theta_sq = 1.0e-16; };  //to avoid singularities 
    
    H_0(3,3) = (cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)*-2.0+(cos(theta)*(ax(1)*Umrz-ax(2)*Umry))/(theta_sq)-sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)-sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umrx*sin(theta)*(ax(1)*Umry+ax(2)*Umrz)*1.0/pow(theta_sq,3.0/2.0)*2.0-Umrx*(ax(1)*Umry+ax(2)*Umrz)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*4.0-(Umrx*Umrx)*cos(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,2.0)*3.0-(Umrx*Umrx)*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)+(Umrx*Umrx)*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,5.0/2.0)*3.0-(Umrx*Umrx)*cos(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)+(Umrx*Umrx)*sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+(Umrx*Umrx)*(cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_0(3,4) = (ax(1)*(cos(theta)-1.0))/(theta_sq)-(ax(2)*Umrx*cos(theta))/(theta_sq)+ax(2)*Umrx*sin(theta)*1.0/pow(theta_sq,3.0/2.0)+Umrx*sin(theta)*(ax(0)*Umry*2.0-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-Umry*sin(theta)*(ax(1)*Umry+ax(2)*Umrz)*1.0/pow(theta_sq,3.0/2.0)+Umrx*(ax(0)*Umry*2.0-ax(1)*Umrx)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umry*(ax(1)*Umry+ax(2)*Umrz)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umrx*Umry*cos(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,2.0)*3.0-Umrx*Umry*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)+Umrx*Umry*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,5.0/2.0)*3.0-Umrx*Umry*cos(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)+Umrx*Umry*sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+Umrx*Umry*(cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_0(3,5) = (ax(2)*(cos(theta)-1.0))/(theta_sq)+(ax(1)*Umrx*cos(theta))/(theta_sq)-ax(1)*Umrx*sin(theta)*1.0/pow(theta_sq,3.0/2.0)+Umrx*sin(theta)*(ax(0)*Umrz*2.0-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-Umrz*sin(theta)*(ax(1)*Umry+ax(2)*Umrz)*1.0/pow(theta_sq,3.0/2.0)+Umrx*(ax(0)*Umrz*2.0-ax(2)*Umrx)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umrz*(ax(1)*Umry+ax(2)*Umrz)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umrx*Umrz*cos(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,2.0)*3.0-Umrx*Umrz*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)+Umrx*Umrz*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,5.0/2.0)*3.0-Umrx*Umrz*cos(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)+Umrx*Umrz*sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+Umrx*Umrz*(cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_0(5,3) = H_0(3,5);
    H_0(4,3) = H_0(3,4);
    H_0(4,4) = (ax(0)*(cos(theta)-1.0)*-2.0)/(theta_sq)-(cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)*2.0+(cos(theta)*(ax(1)*Umrz-ax(2)*Umry))/(theta_sq)-sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)-sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0/2.0)-(ax(2)*Umry*cos(theta)*2.0)/(theta_sq)+ax(2)*Umry*sin(theta)*1.0/pow(theta_sq,3.0/2.0)*2.0+Umry*sin(theta)*(ax(0)*Umry*2.0-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)*2.0+Umry*(ax(0)*Umry*2.0-ax(1)*Umrx)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*4.0-(Umry*Umry)*cos(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,2.0)*3.0-(Umry*Umry)*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)+(Umry*Umry)*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,5.0/2.0)*3.0-(Umry*Umry)*cos(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)+(Umry*Umry)*sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+(Umry*Umry)*(cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_0(4,5) = (ax(1)*Umry*cos(theta))/(theta_sq)-(ax(2)*Umrz*cos(theta))/(theta_sq)-ax(1)*Umry*sin(theta)*1.0/pow(theta_sq,3.0/2.0)+ax(2)*Umrz*sin(theta)*1.0/pow(theta_sq,3.0/2.0)+Umrz*sin(theta)*(ax(0)*Umry*2.0-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)+Umry*sin(theta)*(ax(0)*Umrz*2.0-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)+Umrz*(ax(0)*Umry*2.0-ax(1)*Umrx)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0+Umry*(ax(0)*Umrz*2.0-ax(2)*Umrx)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umry*Umrz*cos(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,2.0)*3.0-Umry*Umrz*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)+Umry*Umrz*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,5.0/2.0)*3.0-Umry*Umrz*cos(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)+Umry*Umrz*sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+Umry*Umrz*(cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_0(5,4) = H_0(4,5);
    H_0(5,5) = (ax(0)*(cos(theta)-1.0)*-2.0)/(theta_sq)-(cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)*2.0+(cos(theta)*(ax(1)*Umrz-ax(2)*Umry))/(theta_sq)-sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)-sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0/2.0)+(ax(1)*Umrz*cos(theta)*2.0)/(theta_sq)-ax(1)*Umrz*sin(theta)*1.0/pow(theta_sq,3.0/2.0)*2.0+Umrz*sin(theta)*(ax(0)*Umrz*2.0-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)*2.0+Umrz*(ax(0)*Umrz*2.0-ax(2)*Umrx)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*4.0-(Umrz*Umrz)*cos(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,2.0)*3.0-(Umrz*Umrz)*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)+(Umrz*Umrz)*sin(theta)*(ax(1)*Umrz-ax(2)*Umry)*1.0/pow(theta_sq,5.0/2.0)*3.0-(Umrz*Umrz)*cos(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,2.0)+(Umrz*Umrz)*sin(theta)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+(Umrz*Umrz)*(cos(theta)-1.0)*(-ax(0)*(Umry*Umry+Umrz*Umrz)+ax(1)*Umrx*Umry+ax(2)*Umrx*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    
    
    H_1(3,3) = (ax(1)*(cos(theta)-1.0)*-2.0)/(theta_sq)-(cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*2.0-(cos(theta)*(ax(0)*Umrz-ax(2)*Umrx))/(theta_sq)+sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)+(ax(2)*Umrx*cos(theta)*2.0)/(theta_sq)-ax(2)*Umrx*sin(theta)*1.0/pow(theta_sq,3.0/2.0)*2.0-Umrx*sin(theta)*(ax(0)*Umry-ax(1)*Umrx*2.0)*1.0/pow(theta_sq,3.0/2.0)*2.0-Umrx*(ax(0)*Umry-ax(1)*Umrx*2.0)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*4.0+(Umrx*Umrx)*cos(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,2.0)*3.0+(Umrx*Umrx)*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-(Umrx*Umrx)*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-(Umrx*Umrx)*cos(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+(Umrx*Umrx)*sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+(Umrx*Umrx)*(cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_1(3,4) = (ax(0)*(cos(theta)-1.0))/(theta_sq)+(ax(2)*Umry*cos(theta))/(theta_sq)-ax(2)*Umry*sin(theta)*1.0/pow(theta_sq,3.0/2.0)-Umry*sin(theta)*(ax(0)*Umry-ax(1)*Umrx*2.0)*1.0/pow(theta_sq,3.0/2.0)-Umrx*sin(theta)*(ax(0)*Umrx+ax(2)*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umry*(ax(0)*Umry-ax(1)*Umrx*2.0)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umrx*(ax(0)*Umrx+ax(2)*Umrz)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0+Umrx*Umry*cos(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,2.0)*3.0+Umrx*Umry*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-Umrx*Umry*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-Umrx*Umry*cos(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+Umrx*Umry*sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+Umrx*Umry*(cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_1(3,5) = -(ax(0)*Umrx*cos(theta))/(theta_sq)+(ax(2)*Umrz*cos(theta))/(theta_sq)+ax(0)*Umrx*sin(theta)*1.0/pow(theta_sq,3.0/2.0)-ax(2)*Umrz*sin(theta)*1.0/pow(theta_sq,3.0/2.0)-Umrz*sin(theta)*(ax(0)*Umry-ax(1)*Umrx*2.0)*1.0/pow(theta_sq,3.0/2.0)+Umrx*sin(theta)*(ax(1)*Umrz*2.0-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)-Umrz*(ax(0)*Umry-ax(1)*Umrx*2.0)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0+Umrx*(ax(1)*Umrz*2.0-ax(2)*Umry)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0+Umrx*Umrz*cos(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,2.0)*3.0+Umrx*Umrz*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-Umrx*Umrz*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-Umrx*Umrz*cos(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+Umrx*Umrz*sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+Umrx*Umrz*(cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_1(4,3) = H_1(3,4);
    H_1(4,4) = (cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*-2.0-(cos(theta)*(ax(0)*Umrz-ax(2)*Umrx))/(theta_sq)+sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umry*sin(theta)*(ax(0)*Umrx+ax(2)*Umrz)*1.0/pow(theta_sq,3.0/2.0)*2.0-Umry*(ax(0)*Umrx+ax(2)*Umrz)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*4.0+(Umry*Umry)*cos(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,2.0)*3.0+(Umry*Umry)*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-(Umry*Umry)*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-(Umry*Umry)*cos(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+(Umry*Umry)*sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+(Umry*Umry)*(cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_1(4,5) = (ax(2)*(cos(theta)-1.0))/(theta_sq)-(ax(0)*Umry*cos(theta))/(theta_sq)+ax(0)*Umry*sin(theta)*1.0/pow(theta_sq,3.0/2.0)-Umrz*sin(theta)*(ax(0)*Umrx+ax(2)*Umrz)*1.0/pow(theta_sq,3.0/2.0)+Umry*sin(theta)*(ax(1)*Umrz*2.0-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)-Umrz*(ax(0)*Umrx+ax(2)*Umrz)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0+Umry*(ax(1)*Umrz*2.0-ax(2)*Umry)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0+Umry*Umrz*cos(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,2.0)*3.0+Umry*Umrz*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-Umry*Umrz*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-Umry*Umrz*cos(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+Umry*Umrz*sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+Umry*Umrz*(cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_1(5,3) = H_1(3,5);
    H_1(5,4) = H_1(4,5);
    H_1(5,5) = (ax(1)*(cos(theta)-1.0)*-2.0)/(theta_sq)-(cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*2.0-(cos(theta)*(ax(0)*Umrz-ax(2)*Umrx))/(theta_sq)+sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)-(ax(0)*Umrz*cos(theta)*2.0)/(theta_sq)+ax(0)*Umrz*sin(theta)*1.0/pow(theta_sq,3.0/2.0)*2.0+Umrz*sin(theta)*(ax(1)*Umrz*2.0-ax(2)*Umry)*1.0/pow(theta_sq,3.0/2.0)*2.0+(Umrz*Umrz)*cos(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,2.0)*3.0+Umrz*(ax(1)*Umrz*2.0-ax(2)*Umry)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*4.0+(Umrz*Umrz)*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-(Umrz*Umrz)*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-(Umrz*Umrz)*cos(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+(Umrz*Umrz)*sin(theta)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+(Umrz*Umrz)*(cos(theta)-1.0)*(-ax(1)*(Umrx*Umrx+Umrz*Umrz)+ax(0)*Umrx*Umry+ax(2)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;    
    
    
    H_2(3,3) = (ax(2)*(cos(theta)-1.0)*-2.0)/(theta_sq)-(cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*2.0+(cos(theta)*(ax(0)*Umry-ax(1)*Umrx))/(theta_sq)-sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)-(ax(1)*Umrx*cos(theta)*2.0)/(theta_sq)+ax(1)*Umrx*sin(theta)*1.0/pow(theta_sq,3.0/2.0)*2.0-Umrx*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx*2.0)*1.0/pow(theta_sq,3.0/2.0)*2.0-(Umrx*Umrx)*cos(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,2.0)*3.0-Umrx*(ax(0)*Umrz-ax(2)*Umrx*2.0)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*4.0-(Umrx*Umrx)*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)+(Umrx*Umrx)*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-(Umrx*Umrx)*cos(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+(Umrx*Umrx)*sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+(Umrx*Umrx)*(cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_2(3,4) = (ax(0)*Umrx*cos(theta))/(theta_sq)-(ax(1)*Umry*cos(theta))/(theta_sq)-ax(0)*Umrx*sin(theta)*1.0/pow(theta_sq,3.0/2.0)+ax(1)*Umry*sin(theta)*1.0/pow(theta_sq,3.0/2.0)-Umry*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx*2.0)*1.0/pow(theta_sq,3.0/2.0)-Umrx*sin(theta)*(ax(1)*Umrz-ax(2)*Umry*2.0)*1.0/pow(theta_sq,3.0/2.0)-Umry*(ax(0)*Umrz-ax(2)*Umrx*2.0)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umrx*(ax(1)*Umrz-ax(2)*Umry*2.0)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umrx*Umry*cos(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,2.0)*3.0-Umrx*Umry*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)+Umrx*Umry*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-Umrx*Umry*cos(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+Umrx*Umry*sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+Umrx*Umry*(cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_2(3,5) = (ax(0)*(cos(theta)-1.0))/(theta_sq)-(ax(1)*Umrz*cos(theta))/(theta_sq)+ax(1)*Umrz*sin(theta)*1.0/pow(theta_sq,3.0/2.0)-Umrx*sin(theta)*(ax(0)*Umrx+ax(1)*Umry)*1.0/pow(theta_sq,3.0/2.0)-Umrz*sin(theta)*(ax(0)*Umrz-ax(2)*Umrx*2.0)*1.0/pow(theta_sq,3.0/2.0)-Umrx*(ax(0)*Umrx+ax(1)*Umry)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umrz*(ax(0)*Umrz-ax(2)*Umrx*2.0)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umrx*Umrz*cos(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,2.0)*3.0-Umrx*Umrz*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)+Umrx*Umrz*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-Umrx*Umrz*cos(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+Umrx*Umrz*sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+Umrx*Umrz*(cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_2(4,3) = H_2(3,4);
    H_2(4,4) = (ax(2)*(cos(theta)-1.0)*-2.0)/(theta_sq)-(cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*2.0+(cos(theta)*(ax(0)*Umry-ax(1)*Umrx))/(theta_sq)-sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)+(ax(0)*Umry*cos(theta)*2.0)/(theta_sq)-ax(0)*Umry*sin(theta)*1.0/pow(theta_sq,3.0/2.0)*2.0-Umry*sin(theta)*(ax(1)*Umrz-ax(2)*Umry*2.0)*1.0/pow(theta_sq,3.0/2.0)*2.0-(Umry*Umry)*cos(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,2.0)*3.0-Umry*(ax(1)*Umrz-ax(2)*Umry*2.0)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*4.0-(Umry*Umry)*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)+(Umry*Umry)*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-(Umry*Umry)*cos(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+(Umry*Umry)*sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+(Umry*Umry)*(cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_2(4,5) = (ax(1)*(cos(theta)-1.0))/(theta_sq)+(ax(0)*Umrz*cos(theta))/(theta_sq)-ax(0)*Umrz*sin(theta)*1.0/pow(theta_sq,3.0/2.0)-Umry*sin(theta)*(ax(0)*Umrx+ax(1)*Umry)*1.0/pow(theta_sq,3.0/2.0)-Umrz*sin(theta)*(ax(1)*Umrz-ax(2)*Umry*2.0)*1.0/pow(theta_sq,3.0/2.0)-Umry*(ax(0)*Umrx+ax(1)*Umry)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umrz*(ax(1)*Umrz-ax(2)*Umry*2.0)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*2.0-Umry*Umrz*cos(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,2.0)*3.0-Umry*Umrz*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)+Umry*Umrz*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-Umry*Umrz*cos(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+Umry*Umrz*sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+Umry*Umrz*(cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;
    H_2(5,3) = H_2(3,5);
    H_2(5,4) = H_2(4,5);
    H_2(5,5) = (cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)*-2.0+(cos(theta)*(ax(0)*Umry-ax(1)*Umrx))/(theta_sq)-sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)-sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0/2.0)-Umrz*sin(theta)*(ax(0)*Umrx+ax(1)*Umry)*1.0/pow(theta_sq,3.0/2.0)*2.0-Umrz*(ax(0)*Umrx+ax(1)*Umry)*(cos(theta)-1.0)*1.0/pow(theta_sq,2.0)*4.0-(Umrz*Umrz)*cos(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,2.0)*3.0-(Umrz*Umrz)*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,3.0/2.0)+(Umrz*Umrz)*sin(theta)*(ax(0)*Umry-ax(1)*Umrx)*1.0/pow(theta_sq,5.0/2.0)*3.0-(Umrz*Umrz)*cos(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,2.0)+(Umrz*Umrz)*sin(theta)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,5.0/2.0)*5.0+(Umrz*Umrz)*(cos(theta)-1.0)*(-ax(2)*(Umrx*Umrx+Umry*Umry)+ax(0)*Umrx*Umrz+ax(1)*Umry*Umrz)*1.0/pow(theta_sq,3.0)*8.0;    
    
}


void CRBE2::InitializeHessian() {
    
    // H_0
    MatrixXdDiff H_0_red = MatrixXdDiff::Zero(3,3);
    /*
     H_0_red << 0, ax(2 - 1), ax(3 - 1),
     ax(2 - 1), -2*ax(1 - 1),0,
     ax(3 - 1), 0, -2*ax(1 - 1);
     */
    H_0_red << 0, ax(2 - 1), ax(3 - 1),
            ax(2 - 1), 0  ,0,
            ax(3 - 1), 0, 0     ;    
    
    
    H_0.block(4-1,4-1,3,3) = - H_0_red/2;
    
    // H_1
    MatrixXdDiff H_1_red = MatrixXdDiff::Zero(3,3);
    /*
     H_1_red << 0,ax(2 - 1), ax(3 - 1),
     ax(2 - 1),   0,  0,
     ax(2 - 1), 0 , 0;
     */
    H_1_red << 0,ax(1 - 1), 0,
            ax(1 - 1),   0,  ax(3 - 1),
            0, ax(3 - 1) , 0;
    
    H_1.block(4-1,4-1,3,3) = - H_1_red/2;    
    
    // H_2
    MatrixXdDiff H_2_red = MatrixXdDiff::Zero(3,3);
    /*
     H_2_red << -2*ax(3 - 1), 0, ax(1 - 1),
     0, -2*ax(3 - 1),ax(2 - 1),
     ax(1 - 1), ax(2 - 1), 0;
     */
    H_2_red << 0, 0  , ax(1 - 1), 
            0,   0,  ax(2 - 1),
            ax(1 - 1), ax(2 - 1) , 0;
    
    H_2.block(4-1,4-1,3,3) = - H_2_red/2;     
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
    U_M_rot = Um.segment(4-1,3);      // Rotation at master
    //std::cout << "U_M_rot " <<  " = \n" << U_M_rot  <<std::endl;
    PseudoToRot(U_M_rot, R);
        
    R = R - MatrixXdDiff::Identity(3,3);
    g = Us - Um;
    g.segment(0,3) = g.segment(0,3) -R*ax;
    /*
    Us_p = Um.segment(1-1,3) + R*ax;
    std::cout << "U_M_rot " <<  " = \n" << U_M_rot  <<std::endl;
    std::cout << "g " <<  " = \n" << g  <<std::endl;
    std::cout << "R " <<  " = \n" << R  <<std::endl;
    std::cout << "ax " <<  " = \n" << ax  <<std::endl;
    std::cout << "Us_constraint " <<  " = \n" << Us_p  <<std::endl;
    
    std::ofstream myfile;
    myfile.open ("g.pyBeam", fstream::in | fstream::out | fstream::app);
    myfile << "g = "; myfile << g ; myfile << "\n ";
    myfile.close();
    */
};

CRBE2::~CRBE2(void) {

};