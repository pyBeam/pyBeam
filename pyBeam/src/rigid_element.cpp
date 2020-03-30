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
/*
void CRBE2::EvalJacobian_2(VectorXdDiff Um ) {
    
    
    // Jacobian depends on the rotations of the master node
    
     * G = | -I  -d (R(U_m_rot)*axis_vector)/dU_m_rot  I   0  |
     *     |  0                I                       0  -I  |
     
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
    
    MatrixXdDiff G_old = MatrixXdDiff::Zero(6,12);       // Jacobian of the constraint equations on LHS    
    
    G_old.block(0,0, 3 , 3) =  -MatrixXdDiff::Identity(3,3);
    G_old.block(3,3, 3 , 3) =  -MatrixXdDiff::Identity(3,3);
    
    G_old.block(0,6, 3 , 3) =  MatrixXdDiff::Identity(3,3);
    G_old.block(3,9, 3 , 3) =  MatrixXdDiff::Identity(3,3); 


    G_old.block(0,3, 3 , 3) =  -L;  
    
    std::cout << "G_old " <<  " = \n" << G_old  <<std::endl;    
    //std::cout << "G^T*G " <<  " = \n" << G.transpose()*G  <<std::endl;     
    
}
*/
void CRBE2::EvalJacobian(VectorXdDiff Um ) {
    // Jacobian depends on the rotations of the master node
    /*
     * G = | -I  -d (R(U_m_rot)*axis_vector)/dU_m_rot  I   0  |
     *     |  0                I                       0  -I  |
     */
    G = MatrixXdDiff::Zero(6,12);       // Jacobian of the constraint equations on LHS   
    addouble Dx = ax(0); addouble Dy = ax(1); addouble Dz = ax(2);
    //addouble Umtx = Um(1-1); addouble Umty = Um(2-1); addouble Umtz = Um(3-1);
    addouble Umrx = Um(4-1); addouble Umry = Um(5-1); addouble Umrz = Um(6-1);

    if (abs(Umrx) < 1.0e-19 and abs(Umry) < 1.0e-19 and abs(Umrz) < 1.0e-19 ) {Umrz = 1.0e-19; };  //to avoid singularities 
    addouble theta = Umrx*Umrx+Umry*Umry+Umrz*Umrz;
    G(0,3) = -Dy*(-(Umry*(cos(sqrt(theta))-1.0))/(theta)-(Umrx*Umrz*cos(sqrt(theta)))/(theta)+Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-Dz*(-(Umrz*(cos(sqrt(theta))-1.0))/(theta)+(Umrx*Umry*cos(sqrt(theta)))/(theta)-Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)+Dx*Umrx*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dx*Umrx*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0/2.0);
    G(0,4) = -Dy*(-(Umrx*(cos(sqrt(theta))-1.0))/(theta)-(Umry*Umrz*cos(sqrt(theta)))/(theta)+Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-Dz*(sin(sqrt(theta))*1.0/sqrt(theta)+((Umry*Umry)*cos(sqrt(theta)))/(theta)-(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-(Dx*Umry*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dx*Umry*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dx*Umry*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0/2.0);
    G(0,5) = -Dz*(-(Umrx*(cos(sqrt(theta))-1.0))/(theta)+(Umry*Umrz*cos(sqrt(theta)))/(theta)-Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-Dy*(-sin(sqrt(theta))*1.0/sqrt(theta)-((Umrz*Umrz)*cos(sqrt(theta)))/(theta)+(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-(Dx*Umrz*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dx*Umrz*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dx*Umrz*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0/2.0);
    G(1,3) = -Dx*(-(Umry*(cos(sqrt(theta))-1.0))/(theta)+(Umrx*Umrz*cos(sqrt(theta)))/(theta)-Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-Dz*(-sin(sqrt(theta))*1.0/sqrt(theta)-((Umrx*Umrx)*cos(sqrt(theta)))/(theta)+(Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-(Dy*Umrx*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dy*Umrx*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dy*Umrx*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0/2.0);
    G(1,4) = -Dx*(-(Umrx*(cos(sqrt(theta))-1.0))/(theta)+(Umry*Umrz*cos(sqrt(theta)))/(theta)-Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-Dz*(-(Umrz*(cos(sqrt(theta))-1.0))/(theta)-(Umrx*Umry*cos(sqrt(theta)))/(theta)+Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umry*Umry)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)+Dy*Umry*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dy*Umry*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0/2.0);
    G(1,5) = -Dz*(-(Umry*(cos(sqrt(theta))-1.0))/(theta)-(Umrx*Umrz*cos(sqrt(theta)))/(theta)+Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-Dx*(sin(sqrt(theta))*1.0/sqrt(theta)+((Umrz*Umrz)*cos(sqrt(theta)))/(theta)-(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-(Dy*Umrz*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dy*Umrz*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dy*Umrz*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0/2.0);
    G(2,3) = -Dx*(-(Umrz*(cos(sqrt(theta))-1.0))/(theta)-(Umrx*Umry*cos(sqrt(theta)))/(theta)+Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-Dy*(sin(sqrt(theta))*1.0/sqrt(theta)+((Umrx*Umrx)*cos(sqrt(theta)))/(theta)-(Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-(Dz*Umrx*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dz*Umrx*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)*2.0+Dz*Umrx*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0/2.0);
    G(2,4) = -Dy*(-(Umrz*(cos(sqrt(theta))-1.0))/(theta)+(Umrx*Umry*cos(sqrt(theta)))/(theta)-Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umry*Umry)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-Dx*(-sin(sqrt(theta))*1.0/sqrt(theta)-((Umry*Umry)*cos(sqrt(theta)))/(theta)+(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-(Dz*Umry*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dz*Umry*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)*2.0+Dz*Umry*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0/2.0);
    G(2,5) = -Dx*(-(Umrx*(cos(sqrt(theta))-1.0))/(theta)-(Umry*Umrz*cos(sqrt(theta)))/(theta)+Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)-Dy*(-(Umry*(cos(sqrt(theta))-1.0))/(theta)+(Umrx*Umrz*cos(sqrt(theta)))/(theta)-Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0)+Dz*Umrz*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)*2.0+Dz*Umrz*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0/2.0);
   
    
    G.block(0,0, 3 , 3) =  -MatrixXdDiff::Identity(3,3);
    G.block(3,3, 3 , 3) =  -MatrixXdDiff::Identity(3,3);
    
    G.block(0,6, 3 , 3) =  MatrixXdDiff::Identity(3,3);
    G.block(3,9, 3 , 3) =  MatrixXdDiff::Identity(3,3);
      
}
void CRBE2::EvalHessian(VectorXdDiff Um){
 
    // H_0
    H_0 = MatrixXdDiff::Zero(12,12);
    H_1 = MatrixXdDiff::Zero(12,12);
    H_2 = MatrixXdDiff::Zero(12,12);
    addouble Dx = ax(0); addouble Dy = ax(1); addouble Dz = ax(2);
    addouble Umtx = Um(1-1); addouble Umty = Um(2-1); addouble Umtz = Um(3-1);
    addouble Umrx = Um(4-1); addouble Umry = Um(5-1); addouble Umrz = Um(6-1);
    
     if (abs(Umrx) < 1.0e-19 and abs(Umry) < 1.0e-19 and abs(Umrz) < 1.0e-19 ) {Umrx = 1.0e-19; };  //to avoid singularities 
    //if (abs(Umry) < 1.0e-19 ) {Umry = 1.0e-19; };  //to avoid singularities 
    //if (abs(Umrz) < 1.0e-19 ) {Umrz = 1.0e-19; };  //to avoid singularities
    addouble theta = Umrx*Umrx+Umry*Umry+Umrz*Umrz;
    H_0(3,3) = -Dy*(-(Umrz*cos(sqrt(theta)))/(theta)+Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0+(Umrx*Umrx)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umrx*Umrx*Umrx)*Umry*cos(sqrt(theta))*1.0/pow(theta,2.0)+(Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-(Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrx*Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-(Umrx*Umrx*Umrx)*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)-Dz*((Umry*cos(sqrt(theta)))/(theta)-Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0-(Umrx*Umrx)*Umry*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umrx*Umrx*Umrx)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrx*Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-(Umrx*Umrx*Umrx)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)+Dx*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0/2.0)+Dx*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dx*(Umrx*Umrx)*cos(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)-Dx*(Umrx*Umrx)*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dx*(Umrx*Umrx)*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    H_0(3,4) = -Dz*((Umrx*cos(sqrt(theta)))/(theta)-Umrx*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-Umrx*(Umry*Umry)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0-Umrx*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrx*Umrx)*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+(Umrx*Umrx)*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umrx*Umrx)*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dy*(-(cos(sqrt(theta))-1.0)/(theta)+(Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+(Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-(Umrx*Umrx)*(Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+(Umrx*Umrx)*(Umry*Umry)*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umrx*Umrx)*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0+Umrx*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0)+Dx*Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*2.0+Dx*Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*4.0+Dx*Umrx*Umry*cos(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)-Dx*Umrx*Umry*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dx*Umrx*Umry*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    H_0(3,5) = -Dy*(-(Umrx*cos(sqrt(theta)))/(theta)+Umrx*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+Umrx*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+Umrx*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrx*Umrx)*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+(Umrx*Umrx)*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umrx*Umrx)*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dz*(-(cos(sqrt(theta))-1.0)/(theta)+(Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-(Umrx*Umrx)*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+(Umrx*Umrx)*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umrx*Umrx)*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-Umrx*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0-Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0)+Dx*Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*2.0+Dx*Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*4.0+Dx*Umrx*Umrz*cos(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)-Dx*Umrx*Umrz*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dx*Umrx*Umrz*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    H_0(4,3) = H_0(3,4);
    H_0(4,4) = -Dy*(-(Umrz*cos(sqrt(theta)))/(theta)+Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0+Umrx*(Umry*Umry*Umry)*cos(sqrt(theta))*1.0/pow(theta,2.0)+(Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*(Umry*Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*(Umry*Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)+Dz*((Umry*Umry*Umry)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umry*Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-(Umry*Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umry*cos(sqrt(theta))*3.0)/(theta)+Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0-Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+Umrx*(Umry*Umry)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0-Umrx*(Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)+Umrx*(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-(Dx*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dx*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*4.0+Dx*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0/2.0)+Dx*(Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*8.0+Dx*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dx*(Umry*Umry)*cos(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)-Dx*(Umry*Umry)*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dx*(Umry*Umry)*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    H_0(4,5) = -Dy*(-(Umry*cos(sqrt(theta)))/(theta)+Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+Umry*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*(Umry*Umry)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+Umrx*(Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dz*((Umrz*cos(sqrt(theta)))/(theta)-Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-(Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0-(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*Umry*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+Umrx*Umry*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)+Dx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*4.0+Dx*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*8.0+Dx*Umry*Umrz*cos(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)-Dx*Umry*Umrz*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dx*Umry*Umrz*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    H_0(5,3) = H_0(3,5);
    H_0(5,4) = H_0(4,5);
    H_0(5,5) = -Dy*((Umrz*Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umrz*Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-(Umrz*Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrz*cos(sqrt(theta))*3.0)/(theta)+Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-Umrx*Umry*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+Umrx*Umry*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dz*((Umry*cos(sqrt(theta)))/(theta)-Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0+Umrx*(Umrz*Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umry*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0-Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*(Umrz*Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0+Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*(Umrz*Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)-(Dx*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dx*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*4.0+Dx*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0/2.0)+Dx*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*8.0+Dx*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dx*(Umrz*Umrz)*cos(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,2.0)-Dx*(Umrz*Umrz)*sin(sqrt(theta))*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dx*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    
    H_1(3,3) = -Dx*((Umrz*cos(sqrt(theta)))/(theta)-Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0-(Umrx*Umrx)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umrx*Umrx*Umrx)*Umry*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrx*Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-(Umrx*Umrx*Umrx)*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)-Dz*((Umrx*Umrx*Umrx)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umrx*Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-(Umrx*Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrx*cos(sqrt(theta))*3.0)/(theta)+Umrx*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-(Umrx*Umrx)*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+(Umrx*Umrx)*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umrx*Umrx)*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-(Dy*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dy*(Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*4.0+Dy*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0/2.0)+Dy*(Umrx*Umrx)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*8.0+Dy*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dy*(Umrx*Umrx)*cos(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)-Dy*(Umrx*Umrx)*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dy*(Umrx*Umrx)*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    H_1(3,4) = -Dz*(-(Umry*cos(sqrt(theta)))/(theta)+Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+(Umrx*Umrx)*Umry*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-(Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*(Umry*Umry)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+Umrx*(Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dx*(-(cos(sqrt(theta))-1.0)/(theta)+(Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+(Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-(Umrx*Umrx)*(Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+(Umrx*Umrx)*(Umry*Umry)*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umrx*Umrx)*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-Umrx*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0-Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0)+Dy*Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*2.0+Dy*Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*4.0+Dy*Umrx*Umry*cos(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)-Dy*Umrx*Umry*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dy*Umrx*Umry*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    H_1(3,5) = -Dx*((Umrx*cos(sqrt(theta)))/(theta)-Umrx*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-Umrx*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0-Umrx*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrx*Umrx)*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+(Umrx*Umrx)*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umrx*Umrx)*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dz*(-(Umrz*cos(sqrt(theta)))/(theta)+Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+(Umrx*Umrx)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-(Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*Umry*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+Umrx*Umry*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)+Dy*Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*4.0+Dy*Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*8.0+Dy*Umrx*Umrz*cos(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)-Dy*Umrx*Umrz*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dy*Umrx*Umrz*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    H_1(4,3) = H_1(3,4);
    H_1(4,4) = -Dx*((Umrz*cos(sqrt(theta)))/(theta)-Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0+Umrx*(Umry*Umry*Umry)*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0-(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*(Umry*Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0+(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*(Umry*Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)-Dz*(-(Umrx*cos(sqrt(theta)))/(theta)+Umrx*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0+Umrx*(Umry*Umry)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umry*Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)+Umrx*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umry*Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-(Umry*Umry*Umry)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)+Dy*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0/2.0)+Dy*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dy*(Umry*Umry)*cos(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)-Dy*(Umry*Umry)*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dy*(Umry*Umry)*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    H_1(4,5) = -Dx*((Umry*cos(sqrt(theta)))/(theta)-Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-Umry*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0-Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*(Umry*Umry)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+Umrx*(Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dz*(-(cos(sqrt(theta))-1.0)/(theta)+(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-(Umry*Umry)*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+(Umry*Umry)*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umry*Umry)*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0+Umrx*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0)+Dy*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*2.0+Dy*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*4.0+Dy*Umry*Umrz*cos(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)-Dy*Umry*Umrz*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dy*Umry*Umrz*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    H_1(5,3) = H_1(3,5);
    H_1(5,4) = H_1(4,5);
    H_1(5,5) = Dx*((Umrz*Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umrz*Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-(Umrz*Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrz*cos(sqrt(theta))*3.0)/(theta)+Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0-Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+Umrx*Umry*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0-Umrx*Umry*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)+Umrx*Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dz*(-(Umrx*cos(sqrt(theta)))/(theta)+Umrx*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0+Umrx*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+Umry*(Umrz*Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)+Umrx*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umry*(Umrz*Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-Umry*(Umrz*Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)-(Dy*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dy*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*4.0+Dy*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0/2.0)+Dy*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*8.0+Dy*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)*2.0+Dy*(Umrz*Umrz)*cos(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,2.0)-Dy*(Umrz*Umrz)*sin(sqrt(theta))*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,5.0/2.0)*5.0-Dy*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz)*1.0/pow(theta,3.0)*8.0;
    
    H_2(3,3) = -Dx*(-(Umry*cos(sqrt(theta)))/(theta)+Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0+(Umrx*Umrx)*Umry*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umrx*Umrx*Umrx)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)+(Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-(Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrx*Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-(Umrx*Umrx*Umrx)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)+Dy*((Umrx*Umrx*Umrx)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umrx*Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-(Umrx*Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrx*cos(sqrt(theta))*3.0)/(theta)+Umrx*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0-Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+(Umrx*Umrx)*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0-(Umrx*Umrx)*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)+(Umrx*Umrx)*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-(Dz*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dz*(Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*4.0+Dz*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0/2.0)+Dz*(Umrx*Umrx)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*8.0+Dz*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)*2.0+Dz*(Umrx*Umrx)*cos(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)-Dz*(Umrx*Umrx)*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,5.0/2.0)*5.0-Dz*(Umrx*Umrx)*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0)*8.0;
    H_2(3,4) = -Dx*(-(Umrx*cos(sqrt(theta)))/(theta)+Umrx*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+Umrx*(Umry*Umry)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+Umrx*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umrx*Umrx)*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+(Umrx*Umrx)*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umrx*Umrx)*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dy*((Umry*cos(sqrt(theta)))/(theta)-Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-(Umrx*Umrx)*Umry*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0-(Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umry*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*(Umry*Umry)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+Umrx*(Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)+Dz*Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*4.0+Dz*Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*8.0+Dz*Umrx*Umry*cos(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)-Dz*Umrx*Umry*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,5.0/2.0)*5.0-Dz*Umrx*Umry*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0)*8.0;
    H_2(3,5) = -Dy*((Umrz*cos(sqrt(theta)))/(theta)-Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-(Umrx*Umrx)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0-(Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*Umry*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+Umrx*Umry*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dx*(-(cos(sqrt(theta))-1.0)/(theta)+(Umrx*Umrx)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrx*Umrx)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-(Umrx*Umrx)*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+(Umrx*Umrx)*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umrx*Umrx)*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0+Umrx*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0)+Dz*Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*2.0+Dz*Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*4.0+Dz*Umrx*Umrz*cos(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)-Dz*Umrx*Umrz*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,5.0/2.0)*5.0-Dz*Umrx*Umrz*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0)*8.0;
    H_2(4,3) = H_2(3,4);
    H_2(4,4) = -Dx*((Umry*Umry*Umry)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umry*Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-(Umry*Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umry*cos(sqrt(theta))*3.0)/(theta)+Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-Umrx*(Umry*Umry)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+Umrx*(Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dy*((Umrx*cos(sqrt(theta)))/(theta)-Umrx*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0-Umrx*(Umry*Umry)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umry*Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-(Umry*Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-(Umry*Umry*Umry)*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)-(Dz*(cos(sqrt(theta))-1.0)*2.0)/(theta)+Dz*(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*4.0+Dz*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0/2.0)+Dz*(Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*8.0+Dz*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)*2.0+Dz*(Umry*Umry)*cos(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)-Dz*(Umry*Umry)*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,5.0/2.0)*5.0-Dz*(Umry*Umry)*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0)*8.0;
    H_2(4,5) = -Dx*(-(Umrz*cos(sqrt(theta)))/(theta)+Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+(Umry*Umry)*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-(Umry*Umry)*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*Umry*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+Umrx*Umry*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0)-Dy*(-(cos(sqrt(theta))-1.0)/(theta)+(Umry*Umry)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+(Umry*Umry)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0+(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*2.0-(Umry*Umry)*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0+(Umry*Umry)*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)-(Umry*Umry)*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-Umrx*Umry*Umrz*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0-Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0)+Dz*Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*2.0+Dz*Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*4.0+Dz*Umry*Umrz*cos(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)-Dz*Umry*Umrz*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,5.0/2.0)*5.0-Dz*Umry*Umrz*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0)*8.0;
    H_2(5,3) = H_2(3,5);
    H_2(5,4) = H_2(4,5);
    H_2(5,5) = -Dx*(-(Umry*cos(sqrt(theta)))/(theta)+Umry*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umrx*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0+Umrx*(Umrz*Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)+Umry*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)-Umrx*(Umrz*Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-Umry*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umrx*(Umrz*Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)-Dy*((Umrx*cos(sqrt(theta)))/(theta)-Umrx*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umry*Umrz*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)*3.0+Umry*Umrz*(cos(sqrt(theta))-1.0)*1.0/pow(theta,2.0)*6.0-Umrx*(Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)*3.0+Umry*(Umrz*Umrz*Umrz)*cos(sqrt(theta))*1.0/pow(theta,2.0)-Umrx*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,3.0/2.0)+Umrx*(Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*3.0-Umry*(Umrz*Umrz*Umrz)*sin(sqrt(theta))*1.0/pow(theta,5.0/2.0)*5.0-Umry*(Umrz*Umrz*Umrz)*(cos(sqrt(theta))-1.0)*1.0/pow(theta,3.0)*8.0)+Dz*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0/2.0)+Dz*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)*2.0+Dz*(Umrz*Umrz)*cos(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,2.0)-Dz*(Umrz*Umrz)*sin(sqrt(theta))*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,5.0/2.0)*5.0-Dz*(Umrz*Umrz)*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry)*1.0/pow(theta,3.0)*8.0;
}

void CRBE2::EvalConstraintEquation(VectorXdDiff Um, VectorXdDiff Us) {
    
    addouble Dx = ax(0); addouble Dy = ax(1); addouble Dz = ax(2);
    addouble Umtx = Um(1-1); addouble Umty = Um(2-1); addouble Umtz = Um(3-1);
    addouble Umrx = Um(4-1); addouble Umry = Um(5-1); addouble Umrz = Um(6-1);
    addouble Ustx = Us(1-1); addouble Usty = Us(2-1); addouble Ustz = Us(3-1);
    addouble Usrx = Us(4-1); addouble Usry = Us(5-1); addouble Usrz = Us(6-1);
    //U_M_rot = Um.segment(4-1,3);      // Rotation at triad A
        
    g = VectorXdDiff::Zero(6);
    //VectorXdDiff gb = VectorXdDiff::Zero(6);   
     if (abs(Umrx) < 1.0e-19 and abs(Umry) < 1.0e-19 and abs(Umrz) < 1.0e-19 ) {Umrx = 1.0e-19; };  //to avoid singularities 
    addouble theta = Umrx*Umrx+Umry*Umry+Umrz*Umrz;
    g(0) = -Umtx+Ustx+Dy*(Umrz*sin(sqrt(theta))*1.0/sqrt(theta)+(Umrx*Umry*(cos(sqrt(theta))-1.0))/(theta))-Dz*(Umry*sin(sqrt(theta))*1.0/sqrt(theta)-(Umrx*Umrz*(cos(sqrt(theta))-1.0))/(theta))-(Dx*(cos(sqrt(theta))-1.0)*(Umry*Umry+Umrz*Umrz))/(theta);
    g(1) = -Umty+Usty-Dx*(Umrz*sin(sqrt(theta))*1.0/sqrt(theta)-(Umrx*Umry*(cos(sqrt(theta))-1.0))/(theta))+Dz*(Umrx*sin(sqrt(theta))*1.0/sqrt(theta)+(Umry*Umrz*(cos(sqrt(theta))-1.0))/(theta))-(Dy*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umrz*Umrz))/(theta);
    g(2) = -Umtz+Ustz+Dx*(Umry*sin(sqrt(theta))*1.0/sqrt(theta)+(Umrx*Umrz*(cos(sqrt(theta))-1.0))/(theta))-Dy*(Umrx*sin(sqrt(theta))*1.0/sqrt(theta)-(Umry*Umrz*(cos(sqrt(theta))-1.0))/(theta))-(Dz*(cos(sqrt(theta))-1.0)*(Umrx*Umrx+Umry*Umry))/(theta);
    g(3) = -Umrx+Usrx;
    g(4) = -Umry+Usry;
    g(5) = -Umrz+Usrz;    
    
    
}


CRBE2::~CRBE2(void) {

};