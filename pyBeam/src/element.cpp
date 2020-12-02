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

    nodeA = nullptr;
    nodeB= nullptr;
    elprop = nullptr;
    input = nullptr;

};

CElement::~CElement(void) {

};

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
    // cout<<"slenderness"<<pow(J0/A,0.5)<<endl;
    EIz = input->GetYoungModulus()*elprop->GetIzz();
    EIy = input->GetYoungModulus()*elprop->GetIyy();
    GJ  = input->GetShear()*elprop->GetJt();
    AE  = input->GetYoungModulus()*elprop->GetA();
    Iyy = elprop->GetIyy();
    Izz = elprop->GetIzz();
    //
    Iyy_b= elprop->GetIyy_b();
    Izz_b= elprop->GetIzz_b();
    A_b  = elprop->GetA_b();
    //
    h      = elprop->Geth();
    C_wb   = elprop->GetC_wb();
    t_sk   = elprop->Gett_sk();
    t_sp   = elprop->Gett_sp();
    n_stiff= elprop->Getn_stiff();
    A_stiff= elprop->GetA_stiff();
    A_fl   = elprop->GetA_fl();    
    
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
    

    if (elprop->GetisWBDV() == 1){
        C_wb = elprop->GetC_wb(); 
        h    = elprop->Geth();
        t_sk = elprop->Gett_sk();
        t_sp = elprop->Gett_sp();
        A_fl = elprop->GetA_fl();    
        A_stiff = elprop->GetA_stiff();  
        A_b= elprop->GetA_b();     
        Iyy_b= elprop->GetIyy_b();
        Izz_b= elprop->GetIzz_b();    }

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
    //    %                0       0      0    -Lcurr    0      0;
    //    %                0       0     -1     0       0      0;
    //    %                0       1      0     0       0      0;
    //    %                0       0     -1     0      Lcurr    0;
    //    %                0       1      0     0       0     Lcurr];

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

//------------------------------------
// Evaluates FEM element matrix (with update)
//------------------------------------
void CElement::ElementElastic_DBG(MatrixXdDiff &Kel) {


    Kel = MatrixXdDiff::Zero(12,12);
    VectorXdDiff diagonal(12);
    diagonal << AE/l_ini  ,  12*EIz/pow(l_ini,3) , 12*EIy/pow(l_ini,3) , GJ/l_ini ,
                       4*EIy/l_ini ,  4*EIz/l_ini , 
               AE/l_ini  ,  12*EIz/pow(l_ini,3) , 12*EIy/pow(l_ini,3) , GJ/l_ini ,
                       4*EIy/l_ini ,  4*EIz/l_ini ;
            
    // Writing the diagonal
    for (unsigned short index = 0; index < 12; index++) {
        Kel(index,index) = diagonal(index);
    }        
    
    Kel(2-1,6-1) =   (6*EIz)/pow(l_ini,2); 
    Kel(3-1,5-1) =  -(6*EIy)/pow(l_ini,2);  
    // SYM
    Kel(5-1,3-1) =  -(6*EIy)/pow(l_ini,2);
    Kel(6-1,2-1) =   (6*EIz)/pow(l_ini,2);
    
    VectorXdDiff shortdiagonal(6);
    shortdiagonal << -AE/l_ini  ,  -12*EIz/pow(l_ini,3) , -12*EIy/pow(l_ini,3) , GJ/l_ini ,
                       2*EIy/l_ini ,  2*EIz/l_ini  ;

    // BLOCK 1,6 -> 7,12
    // Writing the second diagonal
    for (unsigned short index = 0; index < 6; index++) {
        Kel(index,index+6) = shortdiagonal(index);
    } 
    
    Kel(2-1,12-1) =   Kel(2-1,6-1);
    Kel(3-1,11-1) =   Kel(3-1,5-1);        
    Kel(5-1,9-1)  =  -Kel(5-1,3-1);
    Kel(6-1,8-1)  =  -Kel(6-1,2-1); 
    
    // SYM 
    // Writing the second diagonal
    for (unsigned short index = 0; index < 6; index++) {
        Kel(index+6,index) = shortdiagonal(index);
    } 
    Kel(12-1,2-1) = Kel(2-1,12-1); 
    Kel(11-1,3-1) = Kel(3-1,11-1) ;        
    Kel(9-1,5-1) = Kel(5-1,9-1);  
    Kel(8-1,6-1) = Kel(6-1,8-1);      
    
   // BLOCK 7,12 -> 7,12 out of diagonal
    Kel(8-1,12-1) = - Kel(2-1,6-1) ;
    Kel(9-1,11-1) = - Kel(3-1,5-1) ;
    //SYM 
    Kel(12-1,8-1) = Kel(8-1,12-1) ;
    Kel(11-1,9-1) = Kel(9-1,11-1);        
                
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
 *         EvalRotMat   (to be removed in future u[date of the code)
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

    //Rrig = Rprev.transpose() * R;


}


/***************************************************************
 *
 *         EvalRotMat_FP
 *
 ************************************************************/
/*
 * This routine, given the cumulative displacement vector  of the current finite element
 *  (a) calculates the new Rotation matrix R
 *  (b) Initializes Rprev
 *
 * IMPORTANT: X need to be the last values,
 */

void CElement::EvalRotMat_FP(VectorXdDiff &U_AB,  VectorXdDiff  &X_AB)
{

    Vector3dDiff pa = Vector3dDiff::Zero();
    Vector3dDiff pb = Vector3dDiff::Zero();
    Vector3dDiff p= Vector3dDiff::Zero();
    Vector3dDiff pseudo= Vector3dDiff::Zero();    
    Vector3dDiff pseudoA= Vector3dDiff::Zero();
    Vector3dDiff pseudoB= Vector3dDiff::Zero();    
    Matrix3dDiff Rnode= Matrix3dDiff::Zero();

    // New versors
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
    e2_old = R0.block(1-1,2-1,3,1);      // This is the original y in global ref
    pseudoA = U_AB.segment(4-1,3);      // Rotation at triad A  
    PseudoToRot(pseudoA, Rnode);
    pa = Rnode*e2_old;

    //===> Node B

    pseudoB = U_AB.segment(10-1,3);     // Rotation at triad B   
    PseudoToRot(pseudoB ,  Rnode);
    pb = Rnode*e2_old;

    // Auxiliary Vector for building the new e3
    p = 0.5*(pa + pb);
    //p= p/p.norm();

    /*---------------------
     *       e3
     *---------------------*/

    // Find the new e3 versor 
    e3 = e1.cross(p);
    e3 = e3/e3.norm();
    
    /*---------------------
     *       e2
     *---------------------*/

    // Find the new e2 versor 
    e2 = e3.cross(e1);
    e2 = e2/e2.norm();

    Rprev = R;
    // Update

    R.block(1-1,1-1,3,1) = e1.segment(1-1,3);
    R.block(1-1,2-1,3,1) = e2.segment(1-1,3);
    R.block(1-1,3-1,3,1) = e3.segment(1-1,3);

    R.block(4-1,4-1,3,3) = R.block(1-1,1-1,3,3);

    //Rrig = Rprev.transpose() * R;

}


/***************************************************************
 *
 *         EvalInitialRotMat
 *
 ************************************************************/
/*
 * This routine, given the auxiliary vector of the finite element in the undeformed configuration
 * (a) Initializes the Rotation matrix R
 * (b) Initializes Rprev
 */

void CElement::InitializeRotMats()
{
    
    // Initialization
    Rrig = MatrixXdDiff::Zero(6,6);         // Rotation Matrix
    R     = MatrixXdDiff::Identity(6,6);    // ROtation from global to rigid in current deformed conf
    Rprev = MatrixXdDiff::Identity(6,6);    // ROtation from global to rigid in old deformed conf
    R0     = MatrixXdDiff::Identity(6,6);    // ROtation from global to rigid in  undeformed conf

    // Versor in the initial coord system:  
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

    // Find the e3 versor 
    e3 = e1.cross(aux_vector);
    e3 = e3/e3.norm();

    /*---------------------
     *       e2
     *---------------------*/

    // Find the  e2 versor 
    e2 = e3.cross(e1);

    // Update

    R.block(1 - 1, 1 - 1, 3, 1) = e1.segment(1 - 1, 3);
    R.block(1 - 1, 2 - 1, 3, 1) = e2.segment(1 - 1, 3);
    R.block(1 - 1, 3 - 1, 3, 1) = e3.segment(1 - 1, 3);

    R.block(4 - 1, 4 - 1, 3, 3) = R.block(1 - 1, 1 - 1, 3, 3);

    Rprev = R;
    R0 = R;

    //Rrig = Rprev.transpose() * R; 

}



void  CElement::StressRetrieving()
{ 
    
    int n_tot = 4 + n_stiff;  // n_stiff + 4 flanges  
    
    addouble b=(C_wb)/(((n_tot)/2)-1);
    
    
    addouble L_Qxy = -h/2;
    addouble L_Qxz = -C_wb/2; 
    addouble Edim =input->GetYoungModulus_dimensional();
    addouble N    =  Edim*fint(7-1);  
    addouble Qxy  = Edim*fint(8-1);                    // shear y end B and A
    addouble Qxz  = Edim*fint(9-1);                    // shear z end B and A
    addouble Mt   = Edim*fint(10-1);     
    addouble My   = Edim*fint(11-1)-Qxz*(l_curr/2.0);    // moment y, mid section
    addouble Mz   = Edim*fint(12-1)+Qxy*(l_curr/2.0);    // moment z, mid section
    //cout<<"My ="<< My<<endl;
    //cout<<"MZ ="<< Mz<<endl;
    //cout<<"N = "<< N<<endl;
    //cout<<"Qxy = "<< Qxy<<endl;
    //cout<<"Qxz = "<< Qxz<<endl;
      
   
    
    
    VectorXdDiff dsigma_dx   = VectorXdDiff::Zero(n_tot);
    VectorXdDiff axial_load  = VectorXdDiff::Zero(n_tot);
    
    sigma_booms = VectorXdDiff::Zero(n_tot);
    
    /// Calculation of Normal stress absorbed by booms (Navier Formula) 
    int r= ((n_stiff)/2)%2;
    
    if (n_stiff == 0 ){
        sigma_booms(1-1)=(N/A_b) - (Mz/Izz_b)*C_wb*0.5 +(My/Iyy_b)*(h/2.0);
        sigma_booms(2-1)=(N/A_b) + (Mz/Izz_b)*C_wb*0.5 +(My/Iyy_b)*(h/2.0);
        sigma_booms(3-1)=(N/A_b) + (Mz/Izz_b)*C_wb*0.5 -(My/Iyy_b)*(h/2.0);
        sigma_booms(4-1)=(N/A_b) - (Mz/Izz_b)*C_wb*0.5 -(My/Iyy_b)*(h/2.0);

        dsigma_dx(1-1)= -A_fl*(Qxy/Izz_b)*C_wb*0.5  + A_fl*(Qxz/Iyy_b)*(h/2.0);
        dsigma_dx(2-1)=  A_fl*(Qxy/Izz_b)*C_wb*0.5  + A_fl*(Qxz/Iyy_b)*(h/2.0);
        dsigma_dx(3-1)=  A_fl*(Qxy/Izz_b)*C_wb*0.5  - A_fl*(Qxz/Iyy_b)*(h/2.0);
        dsigma_dx(4-1)= -A_fl*(Qxy/Izz_b)*C_wb*0.5  - A_fl*(Qxz/Iyy_b)*(h/2.0);

        axial_load  = sigma_booms*A_fl;} 
    
    else if (r==0) // Even number 
    {
        for (int i=1-1  ; i<= ((n_tot)/4 - 1) ; i += 1){
            sigma_booms(i)                  = (N/A_b) - (Mz/Izz_b)*b*((n_tot/4)-1-i+(1/2)) + (My/Iyy_b)*(h/2);
            sigma_booms( ((n_tot)/4)+i )    = (N/A_b) + (Mz/Izz_b)*b*(i + (1/2))           + (My/Iyy_b)*(h/2);
            sigma_booms( ((n_tot)/2)+i )    = (N/A_b) + (Mz/Izz_b)*b*((n_tot/4)-1-i+(1/2)) - (My/Iyy_b)*(h/2);
            sigma_booms( (((n_tot)*3/4))+i )= (N/A_b) - (Mz/Izz_b)*b*(i + (1/2))           - (My/Iyy_b)*(h/2);
    
            dsigma_dx(i)                       =  -A_stiff*(Qxy/Izz_b)*b*((n_tot/4)-1-i+(1/2))  + A_stiff*(Qxz/Iyy_b)*(h/2);
            dsigma_dx( ((n_tot)/4)+i )         =  A_stiff*(Qxy/Izz_b)*b*(i + (1/2))             + A_stiff*(Qxz/Iyy_b)*(h/2);
            dsigma_dx( ((n_tot)/2)+i )         =  A_stiff*(Qxy/Izz_b)*b*((n_tot/4)-1-i+(1/2))   - A_stiff*(Qxz/Iyy_b)*(h/2);
            dsigma_dx( ((n_tot)*3/4)+i)        =  -A_stiff*(Qxy/Izz_b)*b*(i + (1/2))            - A_stiff*(Qxz/Iyy_b)*(h/2);
        }
        //Take into account the different Spars' Area in the corners
        dsigma_dx(1-1)             =  dsigma_dx(1-1)*(A_fl/A_stiff); 
        dsigma_dx((n_tot/2) - 1)   =  dsigma_dx((n_tot/2) - 1)*(A_fl/A_stiff); 
        dsigma_dx((n_tot/2+1) - 1) =  dsigma_dx((n_tot/2+1) - 1)*(A_fl/A_stiff); 
        dsigma_dx(n_tot - 1)       =  dsigma_dx(n_tot - 1)*(A_fl/A_stiff); 

        axial_load                  = sigma_booms*A_stiff;
        axial_load(1-1)             = axial_load (1 -1)*(A_fl/A_stiff);
        axial_load((n_tot/2) - 1)   = axial_load((n_tot/2) - 1)*(A_fl/A_stiff);
        axial_load((n_tot/2+1) - 1) = axial_load((n_tot/2+1) - 1)*(A_fl/A_stiff); 
        axial_load(n_tot - 1)       = axial_load(n_tot - 1)*(A_fl/A_stiff);
    }
    else{
    //odd number
        for (int i=1-1 ; i<= ((n_tot-2)/4 - 1) ; i += 1){
            sigma_booms(i)                      = (N/A_b) - (Mz/Izz_b)*b*(((n_tot-2)/4)-i)  + (My/Iyy_b)*(h/2);
            sigma_booms((((n_tot-2)/4)+1) +i)   = (N/A_b) + (Mz/Izz_b)*b*(i + 1)            + (My/Iyy_b)*(h/2);
            sigma_booms((((n_tot-2)/2)+1) +i)   = (N/A_b) + (Mz/Izz_b)*b*(((n_tot-2)/4)-i)  - (My/Iyy_b)*(h/2);
            sigma_booms((((n_tot-2)*3/4)+2) +i) = (N/A_b) - (Mz/Izz_b)*b* (i + 1)           - (My/Iyy_b)*(h/2);
     
            dsigma_dx(i)                       =  -A_stiff*(Qxy/Izz_b)*b*(((n_tot-2)/4)-i)  + A_stiff*(Qxz/Iyy_b)*(h/2);
            dsigma_dx((((n_tot-2)/4)+1) +i)    =   A_stiff*(Qxy/Izz_b)*b*(i + 1)            + A_stiff*(Qxz/Iyy_b)*(h/2);
            dsigma_dx((((n_tot-2)/2)+1) +i)    =   A_stiff*(Qxy/Izz_b)*b*(((n_tot-2)/4)-i)  - A_stiff*(Qxz/Iyy_b)*(h/2);
            dsigma_dx((((n_tot-2)*3/4)+2) +i)  =  -A_stiff*(Qxy/Izz_b)*b*(i + 1)            - A_stiff*(Qxz/Iyy_b)*(h/2);
        }
        
        //Take into account the different Spars' Area in the corners
        dsigma_dx(1-1)             =  dsigma_dx(1-1)*(A_fl/A_stiff); 
        dsigma_dx((n_tot/2) - 1)   =  dsigma_dx((n_tot/2) - 1)*(A_fl/A_stiff); 
        dsigma_dx((n_tot/2+1) - 1) =  dsigma_dx((n_tot/2+1) - 1)*(A_fl/A_stiff); 
        dsigma_dx(n_tot - 1)       =  dsigma_dx(n_tot - 1)*(A_fl/A_stiff); 

        sigma_booms((n_tot-2)/4)           = (N/A_b) + (My/Iyy_b)*(h/2);    //upper stiffener on Z-axis
        sigma_booms(((n_tot-2)*3/4)+1)     = (N/A_b) - (My/Iyy_b)*(h/2);   //lower stiffener on Z-axis

        dsigma_dx((n_tot-2)/4)           =  A_stiff*(Qxz/Iyy_b)*(h/2);    //upper stiffener on Z-axis
        dsigma_dx(((n_tot-2)*3/4)+1)     = -A_stiff*(Qxz/Iyy_b)*(h/2);   //lower stiffener on Z-axis
      
        axial_load                  = sigma_booms*A_stiff;
        axial_load(1-1)             = axial_load (1 -1)*(A_fl/A_stiff);
        axial_load((n_tot/2) - 1)   = axial_load((n_tot/2) - 1)*(A_fl/A_stiff);
        axial_load((n_tot/2+1) - 1) = axial_load((n_tot/2+1) - 1)*(A_fl/A_stiff); 
        axial_load(n_tot - 1)       = axial_load(n_tot - 1)*(A_fl/A_stiff);
      }
      
 
    
    /// Shear Flux calculation

    //Solve the equation :
     
    // tau_coeff*tau + dsigma_dx=0 ---> tau= -dsigma_dx*(tau_coeff)^-1
     
    //    % tau_coeff = ( 1    0       0     0      0 ...     -1;                dsigma_dx=( dsigma/dx (1st )
    //    %                -1     1      0     0      0 ...     0;                            dsigma/dx (2nd)
    //    %                0      -1      1     0     0         0;                               .
    //    %                0      0       -1    1     0         0;                               .
    //    %                0      0       0     -1    1         0;                            dsigma/dx (ntot-1 )
    //    %               bh     bh ...  bh     0     0 ...     0];                              M_Q]     
    
    MatrixXdDiff tau_coeff= MatrixXdDiff::Zero(n_tot,n_tot);
    
    VectorXdDiff q  = VectorXdDiff::Zero(n_tot);
   
    addouble M_Q = Mt + Qxz*L_Qxz - Qxy*L_Qxy ;   // total Torque in the section due to Qxz and Qxy
    
    dsigma_dx(n_tot-1)= M_Q;    //index start from 0   
                
    
      
    // fill tau_coeff matrix 
    for (int j=1 -1 ; j<= (n_tot-1) -1; j+=1){
        tau_coeff(j,j)=1;}           // Diagonal
    
  
    for (int jj=1 -1; jj<= (n_tot-2) -1; jj+=1){
        tau_coeff(jj+1,jj)= -1;}      // sub-diagonal   
   
    for (int jjj=1 -1 ; jjj<= (n_tot/2) -1 ; jjj+=1){
        tau_coeff(n_tot-1,jjj)= b*h;}   // last row
   
        
    tau_coeff(1-1,(n_tot)-1)= -1;  // up right corner
       
    // System resolution 
    q = (tau_coeff).fullPivHouseholderQr().solve(-dsigma_dx);
    
    // Tau retrieving
    tau = VectorXdDiff::Zero(n_tot);
    tau.segment(1-1,n_tot/2 -1 )           = q.segment(1-1, n_tot/2 -1) /t_sk;
    tau(n_tot/2 -1)                        = q(n_tot/2 -1) /t_sp;
    tau.segment(n_tot/2+1 -1, n_tot/2 -1 ) = q.segment(n_tot/2+1 -1, n_tot/2 -1 )/t_sk;
    tau(n_tot -1)                          = q(n_tot -1) /t_sp;
          
    
    
    
    
    //--- Section Verification ---------
   //N resultant in the section 
    addouble N_sec =0;
    
    for (int i = 1-1 ;i<=(n_tot) -1 ; i=i+1){
        N_sec=N_sec+axial_load(i);}
    
    //cout<<"N"<<N_sec<<endl;
 /* 
     
    // Tz resultant in the section 
    addouble Tz_sec= -tau( (n_tot/2)-1)*h + tau(n_tot-1)*h ; 
      
    
    //Ty resultant in the section 
    VectorXdDiff Ty_vect = VectorXdDiff::Zero((n_tot/2)-1);   
    Ty_vect.segment(1-1,(n_tot/2)-1)=tau.segment(1-1, (n_tot/2)-1)*b - tau.segment((n_tot/2), (n_tot/2)-1)*b;
   
    addouble Ty_sec=0;
    for (int iy = 1-1 ;iy<=((n_tot/2)-1) -1 ; iy=iy+1){
        Ty_sec=Ty_sec + Ty_vect(iy);            
    }   
  */ 
    //cout <<"sigma_booms"<<sigma_booms<<endl;
} 



void CElement:: VonMises(){
 // Von mises criteria (skin and booms ) ----> (sigma_e/sigma_all) -1 <=0
          
     int n_tot = n_stiff+4;                     // n_stiff + 4 flanges
     g_element = VectorXdDiff::Zero(2*n_tot);  //element constraint  equation Von mises    dim = (2*n_tot) 
     
     SF = 1.5;                                   // Safety factor
     sigma_y = 468.5;                          // yielding stress Alluminum 7075
     addouble sigma_adm = sigma_y/SF;         // Allowable Stress
     
    // Constraints
     
    
    for (int i= 1 -1 ; i<= n_tot -1 ; i=i+1)
    {
      g_element(i)=(fabs(sigma_booms(i))/sigma_adm)-1;                     // Normal stress state (Booms)     ---> g=(sigma_x/sigma_all) -1 
      g_element((n_tot+1 - 1) +i )=( sqrt(3.0)*fabs(tau(i))/sigma_adm)-1;  // Pure shear state (Spar and skin)---> g= (sqrt(3)*tau /sigma_all) -1 
    
    }
    // cout<<"__________________________"<<endl;
}



void CElement:: BoomsBuckling(){
    
 // Formula of Vallat for the calculation of sigma_cr of buckling for the concentrated areas 
    
    // has been considered " L " section booms , 
          
     int n_tot = n_stiff+4;                     // n_stiff + 4 flanges
     addouble t_b =1;
     addouble h_stiff = (A_stiff + pow(t_b,2))/(2*t_b);
     addouble h_fl =  (A_fl + pow(t_b,2))/(2*t_b);
     
     
     addouble K = 8.5;   // constant for open stiffeners section 
     addouble beta_fl = h_fl/t_b;
     addouble beta_stiff = h_stiff/t_b;
     
      //cout << "H fl"<< h_fl<<endl;
      //cout << "H stiff"<< h_stiff<<endl;
      
     addouble Edim =input->GetYoungModulus_dimensional();   
     
     
     sigma_y = 468.5;
     
     addouble  sigma_buckl_fl =  sigma_y/(1 + K*beta_fl*sigma_y/Edim );
     addouble  sigma_buckl_stiff =  sigma_y/(1 + K*beta_stiff*sigma_y/Edim );
     
     //cout << "Sigma Critica Buckling FLANGE"<< sigma_buckl_fl<<endl;
     //cout << "Sigma Critica Buckling STIFF"<< sigma_buckl_stiff<<endl;
     
    // Constraints
         
    int n_neg=0;
  
    for (int i= 1 -1 ; i<= n_tot -1 ; i=i+1)
    {
      if (sigma_booms(i) < 0){
          n_neg++ ;
      }
    }
      
    g_buckl_element = VectorXdDiff::Zero(n_neg);
    
    int j=0;
    for (int i= 1 -1 ; i<= n_tot -1 ; i=i+1)
    {
      if (sigma_booms(i) < 0){
    
        if (i == 0 or i== (n_tot/2) - 1 or i==(n_tot/2+1) - 1 or i ==n_tot - 1 ){
             g_buckl_element(j++)=(fabs(sigma_booms(i)/sigma_buckl_fl))-1;
            
        }else
         {
            
      g_buckl_element(j++)=(fabs(sigma_booms(i)/sigma_buckl_stiff))-1;   
     
         } 
      }
    }
          
    
   // cout<<"g_buckling = "<<g_buckl_element<<endl;
    
}

 

addouble CElement:: RetrieveNint(){
    std::cout << fint << std::endl;
    return fint[12-1];
}
