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


#include "../include/element.h"
#include <iostream>

CElement::CElement(int element_ID) {

    iElement = element_ID;

    GlobalDOFs = VectorXdDiff::Zero(12);  // Global DOFs
    aux_vector = Vector3dDiff::Zero(3);   // Vector that orientates the section of the beam

};

CElement::~CElement(void) {};

void CElement::setGlobalDOFs(){

    int nodeIndexA =  nodeA-> GeID();
    int nodeIndexB =  nodeB-> GeID();
    int i = 0;
    for (i = 0; i < 6; i++){
        GlobalDOFs(i) = (nodeIndexA)*6 + i;
        GlobalDOFs(6+i) = (nodeIndexB)*6 + i;
    }

};

void CElement::setLength() {

    addouble a = nodeA->GetCoordinate0(0) - nodeB->GetCoordinate0(0);
    addouble b = nodeA->GetCoordinate0(1) - nodeB->GetCoordinate0(1);
    addouble c = nodeA->GetCoordinate0(2) - nodeB->GetCoordinate0(2);

    addouble intermediate = pow(a, 2) + pow(b, 2) + pow(c, 2);

    l_ini = sqrt(intermediate);
};

void CElement::Initializer(CNode* Node1, CNode* Node2, CProperty* Property, CInput* Input,
                           passivedouble AuxVector_x, passivedouble AuxVector_y, passivedouble AuxVector_z){

    // Associate the nodes object
    SetNode_1(Node1);
    SetNode_2(Node2);

    // Calculate element DOFs
    setGlobalDOFs();

    //  Associate property
    SetProperty(Property);

    // Associate all inputs
    SetInput(Input);

    // Set auxiliary vector to orientate the beam section
    SetAuxVector(AuxVector_x, AuxVector_y, AuxVector_z);

    // Set length and element masses
    setLength();
    setElementMass();

    // Store element properties from the input property object
    J0  = elprop->GetJ0();
    A   = elprop->GetA();
    EIz = input->GetYoungModulus()*elprop->GetIzz();
    EIy = input->GetYoungModulus()*elprop->GetIyy();
    GJ  = input->GetShear()*elprop->GetJt();
    AE  = input->GetYoungModulus()*elprop->GetA();
    Iyy = elprop->GetIyy();
    Izz = elprop->GetIzz();

    int elemdofs = 12;

    //Mass matrix initialization (element level)
    Mfem  = MatrixXdDiff::Zero(elemdofs, elemdofs);

    //Internal forces vector initialization (element level)
    fint = VectorXdDiff::Zero(elemdofs);

    // Initializing rotation matrices for the undeformed element
    InitializeRotMats();

    // Length intialization
    l_curr = l_ini;   // Current length
    l_prev = l_ini;   // Previous length

    // Initialize cumulative strain vector
    eps = VectorXdDiff::Zero(6);

    // Initialize stress vector
    phi = VectorXdDiff::Zero(6);

    // INITIALIZATION of KPRIM  (linear)
    Kprim = MatrixXdDiff::Zero(6,6);

    VectorXdDiff diagonal = VectorXdDiff::Zero(6);
    diagonal << AE/l_ini  , GJ/l_ini  ,  4*EIy/l_ini  ,   4*EIz/l_ini , 4*EIy/l_ini , 4*EIz/l_ini ;

    // Writing the diagonal
    for (unsigned short index = 0; index < 6; index++)
    {
        Kprim(index,index) = diagonal(index);
    }

    Kprim(3-1,5-1) = 2*EIy/l_ini;  Kprim(5-1,3-1) = Kprim(3-1,5-1);
    Kprim(4-1,6-1) = 2*EIz/l_ini;  Kprim(6-1,4-1) = Kprim(4-1,6-1);

}


void CElement::SetDependencies(void){

    // Store element properties from the input property object
    J0  = elprop->GetJ0();
    A   = elprop->GetA();
    EIz = input->GetYoungModulus()*elprop->GetIzz();
    EIy = input->GetYoungModulus()*elprop->GetIyy();
    GJ  = input->GetShear()*elprop->GetJt();
    AE  = input->GetYoungModulus()*elprop->GetA();
    Iyy = elprop->GetIyy();
    Izz = elprop->GetIzz();

    // INITIALIZATION of KPRIM  (linear)
    Kprim = MatrixXdDiff::Zero(6,6);

    VectorXdDiff diagonal = VectorXdDiff::Zero(6);
    diagonal << AE/l_ini  , GJ/l_ini  ,  4*EIy/l_ini  ,   4*EIz/l_ini , 4*EIy/l_ini , 4*EIz/l_ini ;

    // Writing the diagonal
    for (unsigned short index = 0; index < 6; index++) {
        Kprim(index,index) = diagonal(index);
    }

    Kprim(3-1,5-1) = 2*EIy/l_ini;  Kprim(5-1,3-1) = Kprim(3-1,5-1);
    Kprim(4-1,6-1) = 2*EIz/l_ini;  Kprim(6-1,4-1) = Kprim(4-1,6-1);

}


//-----------------------------------------------
// Evaluates FEM element matrix according to Rao
//-----------------------------------------------
void CElement::ElementMass_Rao() {
    // This mass matrix is currently evaluated using the original length of the element
    // The Finite Element Method in Engineering- S.S. Rao
    addouble r = J0/A;

    // Element matrix

    // Needed for storing the diagonal
    VectorXdDiff diagonal(12);
    diagonal << 1.0/3.0, 13.0/35.0, 13.0/35.0, r/3.0, pow(l_ini,2)/105.0, pow(l_ini,2)/105.0,
                1.0/3.0, 13.0/35.0, 13.0/35.0, r/3.0, pow(l_ini,2)/105.0, pow(l_ini,2)/105.0;

    Mfem(1-1,7-1) = 1/6.0;
    Mfem(2-1,6-1) = 11/210.0*l_ini;    Mfem(2-1,8-1)=  9/70.0;  Mfem(2-1,12-1)= -13/420.0*l_ini;
    Mfem(3-1,5-1) = -11/210.0*l_ini;   Mfem(3-1,9-1)=  9/70.0;  Mfem(3-1,11-1)=  13/420.0*l_ini;
    Mfem(4-1,10-1)= r/6.0;
    Mfem(5-1,9-1) =  -13/420.0*l_ini ;  Mfem(5-1,11-1) = -pow(l_ini,2)/140.0;
    Mfem(6-1,8-1) = 13/420.0*l_ini;     Mfem(6-1,12-1) = -pow(l_ini,2)/140.0;
    Mfem(8-1,12-1) =  -11/210.0*l_ini;
    Mfem(9-1,11-1) =   11/210.0*l_ini;

    // Writing the diagonal
    for (int iii=0; iii<12; iii++) {
        Mfem(iii,iii) = diagonal(iii);
    }

    // Symmetrizing the matrix
    for (int iii=0; iii<12; iii++) {
        for (int jjj=iii; jjj<12; jjj++)
            Mfem(jjj,iii) = Mfem(iii,jjj);
    }

    // Rescaling with the element's mass
    Mfem = m_e*Mfem;
}

/***************************************************************
 *
 *         Eval Na Nb
 *
 ************************************************************/
/*
 * This routine, for a given value of the actual beam length (l_curr) evaluates
 * the "kinematic" matrices. Recall that
 * fint = ([Na' ; Nb']*Kprim* [Na, Nb]) u;
 */

void CElement::EvalNaNb(MatrixXdDiff &Na, MatrixXdDiff  &Nb) {

    addouble one_to_l = 1.0/l_curr;

    //-------------    KINEMATIC MATRIX  --------------------------------------
    //    % Na=1/Lcurr* [ -Lcurr   0      0     0       0      0;
    //    %                0     0      0    -Lcurr    0      0;
    //    %                0     0     -1     0       0      0;
    //    %                0     1      0     0       0      0;
    //    %                0     0     -1     0      Lcurr    0;
    //    %                0     1      0     0       0     Lcurr];

    Na(1-1,1-1) =   -1.0;    Na(2-1,4-1) =  -1.0;
    Na(3-1,3-1) =   -one_to_l;
    Na(4-1,2-1) =    one_to_l;
    Na(5-1,3-1) =   -one_to_l;     Na(5-1,5-1) = 1.0;
    Na(6-1,2-1) =    one_to_l;     Na(6-1,6-1) = 1.0;

    //    % Nb=1/Lcurr* [  Lcurr    0       0     0     0      0;
    //    %                0      0       0    Lcurr   0      0;
    //    %                0      0       1     0    Lcurr    0;
    //    %                0     -1       0     0     0    Lcurr;
    //    %                0      0       1     0     0      0;
    //    %                0     -1       0     0     0      0];

    Nb(1-1,1-1) =  1.0;
    Nb(2-1,4-1) =  1.0;
    Nb(3-1,3-1) = one_to_l;     Nb(3-1,5-1) = 1.0;
    Nb(4-1,2-1) = -one_to_l;    Nb(4-1,6-1) = 1.0;
    Nb(5-1,3-1) = one_to_l;
    Nb(6-1,2-1) = -one_to_l;

}

//------------------------------------
// Evaluates FEM element matrix (with update)
//------------------------------------
void CElement::ElementElastic_Rao(MatrixXdDiff &Kel) {

    // The Finite Element Method in Engineering- S.S. Rao

    MatrixXdDiff Na = MatrixXdDiff::Zero(6,6);
    MatrixXdDiff Nb = MatrixXdDiff::Zero(6,6);

    EvalNaNb(Na,  Nb);

    // =================   ELASTIC MATRIX
    // KEL   = [Na'; Nb']*Kprim *[ Na Nb];

    Kel.block(1-1,1-1,6,6) = Na.transpose() * Kprim * Na;
    Kel.block(1-1,7-1,6,6) = Na.transpose() * Kprim * Nb;
    Kel.block(7-1,1-1,6,6) = Nb.transpose() * Kprim * Na;
    Kel.block(7-1,7-1,6,6) = Nb.transpose() * Kprim * Nb;

}

/*##############################################
 *
 *    Evaluates FEM tangent element matrix
 *
 *##############################################*/


void CElement::ElementTang_Rao(int iIter, MatrixXdDiff & Ktang){

    MatrixXdDiff Na = MatrixXdDiff::Zero(6,6);
    MatrixXdDiff Nb = MatrixXdDiff::Zero(6,6);

    EvalNaNb(Na,  Nb);

    VectorXdDiff df_dl =  VectorXdDiff::Zero(12);

    //---------------------------------------------
    //          dKel/dl*uel* dl/du
    //---------------------------------------------
    VectorXdDiff dl_du =  VectorXdDiff::Zero(12);
    dl_du(1-1) = -1.0;    dl_du(7-1) = 1.0;

    MatrixXdDiff Kstretch = MatrixXdDiff::Zero(12,12);
    /*
     *
     */
    addouble onetol = 1.0/(l_curr);

    df_dl(2-1) = -onetol*fint(2-1);
    df_dl(3-1) = -onetol*fint(3-1);
    df_dl(8-1) = -onetol*fint(8-1);
    df_dl(9-1) = -onetol*fint(9-1);
    Kstretch = df_dl*dl_du.transpose();

    MatrixXdDiff Kel = MatrixXdDiff::Zero(12,12);
    ElementElastic_Rao(Kel);

    // Element Level tangent Matrix
    // Still needs to be rotated and also added the rigid contribution
    Ktang = Kstretch + Kel;

}


/***************************************************************
 *
 *         EvalRotMat
 *
 ************************************************************/
/*
 * This routine, given the incremental displacement vector  of the current finite element
 * (a) calculates the new Rotation matrix R
 * (b) Find the Rrig (incremental)
 *
 * IMPORTANT: X need to be the last values,
 */

void CElement::EvalRotMat(VectorXdDiff &dU_AB,  VectorXdDiff  &X_AB) {

    Vector3dDiff pa = Vector3dDiff::Zero();
    Vector3dDiff pb = Vector3dDiff::Zero();
    Vector3dDiff p= Vector3dDiff::Zero();
    Vector3dDiff pseudo= Vector3dDiff::Zero();
    Vector3dDiff pseudoA= Vector3dDiff::Zero();
    Vector3dDiff pseudoB= Vector3dDiff::Zero();
    Matrix3dDiff Rnode= Matrix3dDiff::Zero();

    // New versor in old local coord system
    Vector3dDiff e1 = Vector3dDiff::Zero();
    Vector3dDiff e2 = Vector3dDiff::Zero();
    Vector3dDiff e3 = Vector3dDiff::Zero();

    /*---------------------
     *       e1
     *---------------------*/

    e1 = X_AB.tail(3) - X_AB.head(3);
    e1 = e1/e1.norm();

    /*---------------------
     *      p
     *---------------------*/

    //===> Node A
    Vector3dDiff e2_old = Vector3dDiff::Zero();
    e2_old = R.block(1-1,2-1,3,1);      // This is the old y in global ref
    pseudoA = dU_AB.segment(4-1,3);      // Rotation at triad A
    PseudoToRot(pseudoA, Rnode);
    pa = Rnode*e2_old;
    

    //===> Node B
    pseudoB = dU_AB.segment(10-1,3);     // Rotation at triad B
    PseudoToRot(pseudoB ,  Rnode);
    pb = Rnode*e2_old;

    // Auxiliary Vector for building the new e3
    p = 0.5*(pa + pb);
    //p= p/p.norm();
    

    /*---------------------
     *       e3
     *---------------------*/

    // Find the new e3 versor (in old local CS)
    e3 = e1.cross(p);
    e3 = e3/e3.norm();

    /*---------------------
     *       e2
     *---------------------*/

    // Find the new e2 versor (in old local CS)
    e2 = e3.cross(e1);
    e2 = e2/e2.norm();

    Rprev = R;

    // Update
    R.block(1-1,1-1,3,1) = e1.segment(1-1,3);
    R.block(1-1,2-1,3,1) = e2.segment(1-1,3);
    R.block(1-1,3-1,3,1) = e3.segment(1-1,3);

    R.block(4-1,4-1,3,3) = R.block(1-1,1-1,3,3);

    Rrig = Rprev.transpose() * R;


}


/***************************************************************
 *
 *         EvalRotMat_FP
 *
 ************************************************************/
/*
 * This routine, given the incremental displacement vector  of the current finite element
 * (a) calculates the new Rotation matrix R
 * (b) Find the Rrig (incremental)
 *
 * IMPORTANT: X need to be the last values,
 */

void CElement::EvalRotMat_FP(VectorXdDiff &dU_AB,  VectorXdDiff  &X_AB)
{

    Vector3dDiff pa = Vector3dDiff::Zero();
    Vector3dDiff pb = Vector3dDiff::Zero();
    Vector3dDiff p= Vector3dDiff::Zero();
    Vector3dDiff pseudo= Vector3dDiff::Zero();    
    Vector3dDiff pseudoA= Vector3dDiff::Zero();
    Vector3dDiff pseudoB= Vector3dDiff::Zero();    
    Matrix3dDiff Rnode= Matrix3dDiff::Zero();

    // New versor in old local coord system
    Vector3dDiff e1 = Vector3dDiff::Zero();
    Vector3dDiff e2 = Vector3dDiff::Zero();
    Vector3dDiff e3 = Vector3dDiff::Zero();

    /*---------------------
     *       e1
     *---------------------*/

    e1 = X_AB.tail(3) - X_AB.head(3);
    e1 = e1/e1.norm();

    /*---------------------
     *      p
     *---------------------*/

    //===> Node A

    Vector3dDiff e2_old = Vector3dDiff::Zero();
    e2_old = R0.block(1-1,2-1,3,1);      // This is the old y in global ref
    pseudoA = dU_AB.segment(4-1,3);      // Rotation at triad A  
    PseudoToRot(pseudoA, Rnode);
    pa = Rnode*e2_old;

    //===> Node B

    pseudoB = dU_AB.segment(10-1,3);     // Rotation at triad B   
    PseudoToRot(pseudoB ,  Rnode);
    pb = Rnode*e2_old;

    // Auxiliary Vector for building the new e3
    p = 0.5*(pa + pb);
    //p= p/p.norm();

    /*---------------------
     *       e3
     *---------------------*/

    // Find the new e3 versor (in old local CS)
    e3 = e1.cross(p);
    e3 = e3/e3.norm();
    
    /*---------------------
     *       e2
     *---------------------*/

    // Find the new e2 versor (in old local CS)
    e2 = e3.cross(e1);
    e2 = e2/e2.norm();

    Rprev = R;
    // Update

    R.block(1-1,1-1,3,1) = e1.segment(1-1,3);
    R.block(1-1,2-1,3,1) = e2.segment(1-1,3);
    R.block(1-1,3-1,3,1) = e3.segment(1-1,3);

    R.block(4-1,4-1,3,3) = R.block(1-1,1-1,3,3);

    Rrig = Rprev.transpose() * R;

}


/***************************************************************
 *
 *         EvalInitialRotMat
 *
 ************************************************************/
/*
 * This routine, given the auxiliary vector of the finite element in the undeformed configuration
 * (a) Initializes the Rotation matrix R
 * (b) Initializes the rigid Rotation matrix Rrig
 * (c) Initializes Rprev
 */

void CElement::InitializeRotMats()
{
    
    // Initialization
    Rrig = MatrixXdDiff::Zero(6,6);         // Rotation Matrix
    R     = MatrixXdDiff::Identity(6,6);    // ROtation from global to rigid in current deformed conf
    Rprev = MatrixXdDiff::Identity(6,6);    // ROtation from global to rigid in old deformed conf
    R0     = MatrixXdDiff::Identity(6,6);    // ROtation from global to rigid in  undeformed conf

    // New versor in old local coord system
    Vector3dDiff e1 = Vector3dDiff::Zero();
    Vector3dDiff e2 = Vector3dDiff::Zero();
    Vector3dDiff e3 = Vector3dDiff::Zero();

    /*---------------------
     *       e1  (along element exis)
     *---------------------*/
    e1(1 -1,0) = nodeB->GetCoordinate0(1 - 1) - nodeA->GetCoordinate0(1 - 1);
    e1(2 -1,0) = nodeB->GetCoordinate0(2 - 1) - nodeA->GetCoordinate0(2 - 1);
    e1(3 -1,0) = nodeB->GetCoordinate0(3 - 1) - nodeA->GetCoordinate0(3 - 1);
    e1 = e1/e1.norm();

    /*---------------------
     *       e3 (local z axis)
     *---------------------*/

    // Find the new e3 versor (in old local CS)
    e3 = e1.cross(aux_vector);
    e3 = e3/e3.norm();

    /*---------------------
     *       e2
     *---------------------*/

    // Find the  e2 versor (in old local CS)
    e2 = e3.cross(e1);

    // Update

    R.block(1 - 1, 1 - 1, 3, 1) = e1.segment(1 - 1, 3);
    R.block(1 - 1, 2 - 1, 3, 1) = e2.segment(1 - 1, 3);
    R.block(1 - 1, 3 - 1, 3, 1) = e3.segment(1 - 1, 3);

    R.block(4 - 1, 4 - 1, 3, 3) = R.block(1 - 1, 1 - 1, 3, 3);

    Rprev = R;
    R0 = R;

    Rrig = Rprev.transpose() * R;

}

