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


#include "../include/element.h"
#include <iostream>


CElement::CElement(int element_ID) { iElement = element_ID; };

void CElement::setGlobalDOFs(){
    
    int nodeIndexA =  nodeA-> GeID();
    int nodeIndexB =  nodeB-> GeID();
    int i = 0;
    for (int i=1; i<=6; i++){
        GlobalDOFs(i-1) = (nodeIndexA-1)*6 + i -1;      
        GlobalDOFs(6+i-1) = (nodeIndexB-1)*6 + i -1;       
    };
    
};  

void CElement::setLength() {
    addouble a = nodeA->GetCoordinate0(0) - nodeB->GetCoordinate0(0);
    addouble b = nodeA->GetCoordinate0(1) - nodeB->GetCoordinate0(1);
    addouble c = nodeA->GetCoordinate0(2) - nodeB->GetCoordinate0(2);
    addouble intermediate = pow(a ,2) + pow(b,2) + pow( c ,2) ;
    l_ini =  sqrt(intermediate );
    //std::cout<< "nodeA-> GeID() = " << nodeA-> GeID() << std::endl;
    //std::cout<< "nodeB-> GeID() = " << nodeB-> GeID() << std::endl;
    //std::cout<< "nodeA->GetCoordinate0(0) = " << nodeA->GetCoordinate0(0) << std::endl;
    //std::cout<< "nodeA->GetCoordinate0(1) = " << nodeA->GetCoordinate0(1) << std::endl;
    //std::cout<< "nodeA->GetCoordinate0(2) = " << nodeA->GetCoordinate0(2) << std::endl;
    //std::cout<< "nodeB->GetCoordinate0(0) = " << nodeB->GetCoordinate0(0) << std::endl;
    //std::cout<< "nodeB->GetCoordinate0(1) = " << nodeB->GetCoordinate0(1) << std::endl;
    //std::cout<< "nodeB->GetCoordinate0(2) = " << nodeB->GetCoordinate0(2) << std::endl;
    //std::cout<< "a = " << a << std::endl;
    //std::cout<< "b = " << b << std::endl;
    //std::cout<< "c = " << c << std::endl;
};

void CElement::setElementMass(){
    // Still don't get why I need two variables for that
    //referred to the undeformed configuration
    m = property->GetA()*l_ini* input->GetDensity();
    m_e = property->GetA()*l_ini* input->GetDensity();
}; 


void CElement::Initializer(CNode* Node1, CNode* Node2, CProperty* Property, CInput* Input, passivedouble AuxVector_x, passivedouble AuxVector_y, passivedouble AuxVector_z){
    
    // Associate the nodes object   
    SetNode_1( Node1) ;
    SetNode_2( Node2);  
    // Calculate element DOFs
    setGlobalDOFs();
    //  Associate property
    SetProperty(Property);
    // Associate all inputs
    SetInput(Input);        
    // set auxiliary vector        
    SetAuxVector(AuxVector_x, AuxVector_y, AuxVector_z);
    
    
    setLength();
    setElementMass();
    J0 = property->GetJ0();
    A = property->GetA();
    EIz = input->GetYoungModulus()*property->GetIzz();
    EIy = input->GetYoungModulus()*property->GetIyy();
    GJ = input->GetShear()*property->GetJt();
    AE = input->GetYoungModulus()*property->GetA();
    Iyy = property->GetIyy();
    Izz = property->GetIzz();
        
    
    //Mass matrix initialization (element level)
    Mfem  = MatrixXdDiff::Zero(elemdofs,elemdofs);   
    //Internal forces vector initialization (element level)
    fint = VectorXdDiff::Zero(elemdofs);    

    // Initializing rotation matrices for the undeformed element
    InitializeRotMats();  
 
    // Length intialization
    l_act  = l_ini;  // This is the current length
    //l_ini  = le;  // This is the original
    l_prev = l_ini;  // This is the previous length  
    
    eps  = VectorXdDiff::Zero(6);   // Elastic Cumulative deformation
    phi  = VectorXdDiff::Zero(6);   // Elastic Cumulative tension
 
        
    // INITIALIZATION of KPRIM  (linear)
    
    VectorXdDiff diagonale = VectorXdDiff::Zero(6);
    Kprim = MatrixXdDiff::Zero(6,6);
    
    diagonale << AE/l_ini  , GJ/l_ini  ,  4*EIy/l_ini  ,   4*EIz/l_ini , 4*EIy/l_ini , 4*EIz/l_ini ;
    
    // Writing the diagonal
    for (unsigned short index=0; index < 6; index++)
    {
        Kprim(index,index) = diagonale(index);
    }
    
    Kprim(3-1,5-1) = 2*EIy/l_ini;  Kprim(5-1,3-1) = Kprim(3-1,5-1);
    Kprim(4-1,6-1) = 2*EIz/l_ini;  Kprim(6-1,4-1) = Kprim(4-1,6-1);
    
    //std::cout<< "iElement = \n" << iElement << std::endl;
    //std::cout<< "l_ini = \n" << l_ini << std::endl;
    /*if (iElement == 2){
    cout.precision(17);
    std::cout<< "Jt = \n" << property->GetJt() << std::endl;
    std::cout<< "Izz = \n" << property->GetIzz() << std::endl;
    std::cout<< "Iyy = \n" << property->GetIyy() << std::endl; 
    std::cout<< "A = \n" << property->GetA() << std::endl; } 
    //std::cout<< "Kprim = \n" << Kprim << std::endl; */
     
    
}

//-----------------------------------------------
// Evaluates FEM element matrix according to Rao
//-----------------------------------------------
void CElement::ElementMass_Rao()
{
    // This mass matrix is currently evaluated using the original length of the element
    // The Finite Element Method in Engineering- S.S. Rao
    addouble r=J0/A;
    
    // Element matrix
    
    // Needed for storing the diagonal
    VectorXdDiff diagonale(12);
    diagonale << 1.0/3.0, 13.0/35.0, 13.0/35.0, r/3.0, pow(l_ini,2)/105.0, pow(l_ini,2)/105.0,
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
    for (int iii=0; iii<12; iii++)
    {
    	Mfem(iii,iii) = diagonale(iii);
    }
    
    
    // Symmetrizing the matrix
    for (int iii=0; iii<12; iii++)
    {
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
 * This routine, for a given value of the actual beam length (l_act) evaluates
 * the "kinematic" matrices. Recall that
 * fint = ([Na' ; Nb']*Kprim* [Na, Nb]) u;
 */

void CElement::EvalNaNb(MatrixXdDiff &Na,  MatrixXdDiff  &Nb)

{
    
    addouble one_to_l = 1.0/l_act;
    
    //-------------    KINEMATIC MATRIX  --------------------------------------
    //    % Na=1/Lact* [ -Lact   0      0     0       0      0;
    //    %                0     0      0    -Lact    0      0;
    //    %                0     0     -1     0       0      0;
    //    %                0     1      0     0       0      0;
    //    %                0     0     -1     0      Lact    0;
    //    %                0     1      0     0       0     Lact];
    
    
    
    Na(1-1,1-1) =   -1.0;    Na(2-1,4-1) =  -1.0;
    Na(3-1,3-1) =   -one_to_l;
    Na(4-1,2-1) =    one_to_l;
    Na(5-1,3-1) =   -one_to_l;     Na(5-1,5-1) = 1.0;
    Na(6-1,2-1) =    one_to_l;     Na(6-1,6-1) = 1.0;
    
    //    % Nb=1/Lact* [  Lact    0       0     0     0      0;
    //    %                0      0       0    Lact   0      0;
    //    %                0      0       1     0    Lact    0;
    //    %                0     -1       0     0     0    Lact;
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
void CElement::ElementElastic_Rao(MatrixXdDiff &Kel)

{
    
    // The Finite Element Method in Engineering- S.S. Rao
    
    MatrixXdDiff Na = MatrixXdDiff::Zero(6,6);
    MatrixXdDiff Nb = MatrixXdDiff::Zero(6,6);
    
    EvalNaNb(Na,  Nb);
    
    // =================   ELASTIC MATRIX
    // KEL   = [Na'; Nb']*Kprim *[ Na Nb];
    
    Kel.block(1-1,1-1,6,6) = Na.transpose()*Kprim*Na;
    Kel.block(1-1,7-1,6,6) = Na.transpose()*Kprim*Nb;
    Kel.block(7-1,1-1,6,6) = Nb.transpose()*Kprim*Na;
    Kel.block(7-1,7-1,6,6) = Nb.transpose()*Kprim*Nb;
    
    //std::cout << "Kel = \n" << Kel <<std::endl;
    //std::cout << "Kprim = \n" << Kprim <<std::endl;
    //std::cout << "Na = \n" << Na <<std::endl;
    //std::cout << "Nb = \n" << Nb <<std::endl;    
}



/*##############################################
 *
 *    Evaluates FEM tangent element matrix
 *
 *##############################################*/


void CElement::ElementTang_Rao(int iIter, MatrixXdDiff & Ktang)
{
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
    addouble onetol = 1.0/(l_act);
    
    df_dl(2-1) = -onetol*fint(2-1);
    df_dl(3-1) = -onetol*fint(3-1);
    df_dl(8-1) = -onetol*fint(8-1);
    df_dl(9-1) = -onetol*fint(9-1);
    Kstretch = df_dl*dl_du.transpose();
    
    MatrixXdDiff Kel = MatrixXdDiff::Zero(12,12);
    ElementElastic_Rao(Kel);
    
    Ktang = Kstretch + Kel;  // Element Level tangent Matrix (not all! Needs to be rotated and also added the rigid contribution)
    //std::cout << "Kstretch = \n" << Kstretch <<std::endl;
    //std::cout << "Kel = \n" << Kel <<std::endl;

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

void CElement::EvalRotMat(VectorXdDiff &dU_AB,  VectorXdDiff  &X_AB)
{
    
    
    Vector3dDiff pa = Vector3dDiff::Zero();
    Vector3dDiff pb = Vector3dDiff::Zero();
    Vector3dDiff p= Vector3dDiff::Zero();
    Vector3dDiff pseudo= Vector3dDiff::Zero();
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
    e2_old = R.block(1-1,2-1,3,1);    // This is the old y in global ref
    pseudo = dU_AB.segment(4-1,3);      // Rotation at triad A
    PseudoToRot(pseudo ,  Rnode);
    pa = Rnode*e2_old;
    
    //===> Node B
    
    pseudo = dU_AB.segment(10-1,3);      // Rotation at triad A
    PseudoToRot(pseudo ,  Rnode);
    pb = Rnode*e2_old;
    
    // Auxiliary Vector for building the new e3
    p = 0.5*(pa + pb);
    
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
    
    Rprev = R;
    // Update
    
    R.block(1-1,1-1,3,1) = e1.segment(1-1,3);
    R.block(1-1,2-1,3,1) = e2.segment(1-1,3);
    R.block(1-1,3-1,3,1) = e3.segment(1-1,3);
    
    R.block(4-1,4-1,3,3) = R.block(1-1,1-1,3,3);
    
    Rrig= Rprev.transpose()*R;
        
}


void CElement::EvalRotMatDEBUG(VectorXdDiff &dU_AB , VectorXdDiff &X_AB , MatrixXdDiff &R)
{
    
    Vector3dDiff pa = Vector3dDiff::Zero();
    Vector3dDiff pb = Vector3dDiff::Zero();
    Vector3dDiff p= Vector3dDiff::Zero();
    Vector3dDiff pseudo= Vector3dDiff::Zero();
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
    e2_old = R.block(1-1,2-1,3,1);    // This is the old y in global ref  (this may be source of error)
    pseudo = dU_AB.segment(4-1,3);      // Rotation at triad A
    PseudoToRot(pseudo ,  Rnode);
    pa = Rnode*e2_old;
    
    //===> Node B
    
    pseudo = dU_AB.segment(10-1,3);      // Rotation at triad A
    PseudoToRot(pseudo ,  Rnode);
    pb = Rnode*e2_old;
    
    // Auxiliary Vector for building the new e3
    p = 0.5*(pa + pb);
    
    
    /*---------------------
     *       e3
     *---------------------*/
    
    
    // Find the new e3 versor (in old local CS)
    e3 = e1.cross(p);
    e3 = e3/e3.norm();
    
    e3(1-1) = 0.0;  e3(2-1) = 0.0;    e3(3-1) = 1.0;
    
    
    /*---------------------
     *       e2
     *---------------------*/
    
    // Find the new e2 versor (in old local CS)
    e2 = e3.cross(e1);
    
    
    // Update
    R = MatrixXdDiff::Zero(6,6);
    
    
    R.block(1-1,1-1,3,1) = e1.segment(1-1,3);
    R.block(1-1,2-1,3,1) = e2.segment(1-1,3);
    R.block(1-1,3-1,3,1) = e3.segment(1-1,3);
    
    R.block(4-1,4-1,3,3) = R.block(1-1,1-1,3,3);
    
    
    
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
    
    R.block(1-1,1-1,3,1) = e1.segment(1-1,3);
    R.block(1-1,2-1,3,1) = e2.segment(1-1,3);
    R.block(1-1,3-1,3,1) = e3.segment(1-1,3);
    
    R.block(4-1,4-1,3,3) = R.block(1-1,1-1,3,3);
    
    Rprev = R;    
    
    Rrig= Rprev.transpose()*R;

}


CElement::~CElement(void) {};
