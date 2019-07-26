/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rocco Bombardieri, Rauno Cavallaro 
 *
 * File developers: Rocco Bombardieri (Carlos III University Madrid)
 *                  Rauno Cavallaro (Carlos III University Madrid)
 *                  Ruben Sanchez (SciComp, TU Kaiserslautern)
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

#include <iostream>
#include <fstream>
#include <chrono>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "../include/beam.h"

using namespace std;
using namespace Eigen;

CBeamSolver::CBeamSolver(void) {
    
    //==============================================================
    //      Input initialization
    //==============================================================
    
    register_loads = false;
    objective_function = 0.0;
    
}

CBeamSolver::~CBeamSolver(void) {
    
}

void CBeamSolver::InitializeInput(CInput* py_input){   // insert node class and connectivity
    
    // I'm memorizing as a member variable the object input passed from outside
    input = py_input;

    if (input->GetDiscreteAdjoint()){
        input->RegisterInput_E();
        input->RegisterInput_Nu();
    }

    input->SetParameters();

    nDOF = input->Get_nDOF();
    nFEM = input->Get_nFEM();
    nTotalDOF = input->Get_nNodes() * input->Get_nDOF();
    nRBE2 = input->Get_nRBE2();
    
    //==============================================================
    //      Load Vector initialization
    //==============================================================
    
    loadVector = new addouble[nTotalDOF];
    for (int iLoad = 0; iLoad < nTotalDOF; iLoad++)
        loadVector[iLoad] = 0.0;
    
    //==============================================================
    //      Node vector initialization
    //==============================================================
    
    cout << "--> Node Structure Initialization... ";
    unsigned long nNodes = input->Get_nNodes();   //substitute from mesh file
    node = new CNode*[nNodes];
    std::cout << "ok!"  <<std::endl;
    
    //==============================================================
    //      Finite Element vector initialization
    //==============================================================
    
    cout << "--> Finite Element Initialization... ";
    element = new CElement*[nFEM];
    std::cout << "ok!"  <<std::endl;
    
    //==============================================================
    //      RBE2 initialization
    //==============================================================
    if (nRBE2 != 0){  
        cout << "--> RBE2 initialization... ";
        RBE2 = new CRBE2*[nRBE2];
        std::cout << "ok!"  <<std::endl;
        std::cout << "--> Warning: the code works if slave nodes are not connected to beams! "  << std::endl;
    }
    else {RBE2 = NULL; }     
    
    
    //===============================================
    //  Initialize structural solver
    //===============================================
    structure = NULL;
    
}

void CBeamSolver::Solve(int FSIIter = 0){
    
    // Beam total length
    addouble TotalLength = 0;
    for  ( unsigned long iFEM = 0; iFEM < nFEM; iFEM++) {
        TotalLength += element[iFEM]->GetInitial_Length();
    }

    std::cout << "--> Setting External Forces" << std::endl;
    structure->ReadForces(nTotalDOF, loadVector);
    
    if (nRBE2 != 0){ 
        std::cout << "--> Setting RBE2 Matrix for Rigid Constraints" << std::endl;
        structure->AddRBE2(input, RBE2);        
    };
        
    //===============================================
    // LOAD STEPPING
    //===============================================
    
    std::cout << "--> Starting Load Stepping" << std::endl << std::endl;

    addouble  lambda = 1.0;
    addouble dlambda =  1.0/input->Get_LoadSteps() ;
    addouble initResNorm   =  1.0;
    addouble initDispNorm  =  1.0;
    unsigned long iIter;
    unsigned long totalIter = 0;
    unsigned long loadStep = 1;
    cout.setf(ios::fixed, ios::floatfield);
    
    // This function set the current initial coordinates and memorizes them as the old one before the converging procedure starts
    structure->InitialCoord();

    for  ( loadStep = 0; loadStep < input->Get_LoadSteps(); loadStep++) {

        lambda = dlambda*(loadStep+1);
        std::cout << "===========================================================================" << std::endl;
        std::cout << "==       LOAD STEP     " << loadStep << std::endl;
        std::cout.precision(8);
        std::cout << "==       Lambda        " << lambda << std::endl;
        std::cout << "===========================================================================" << std::endl;

        std::cout.width(6); std::cout << "Iter";
        std::cout.width(17); std::cout << "Log10(Norm_Res)";
        std::cout.width(17); std::cout << "Log10(Lin_Sol)";
        std::cout.width(17); std::cout << "Log10(Norm_Disp)";
        std::cout.width(17); std::cout << "Log10(Disp_Fact)" << std::endl;
        
        //===============================================
        //               ITERATIVE SEQUENCE
        //===============================================
        bool converged = false;
        
        for (iIter = 0; iIter < input->Get_nIter(); iIter++) {

            std::cout.width(6); std::cout << iIter;

            //std::cout << "   ----- ITERATION  -----" << iIter << std::endl;
            
            /*--------------------------------------------------
             *   Updates  Fext, Residual,
             *----------------------------------------------------*/
            
            // Update the External Forces with the loadStep
            structure->UpdateExtForces(lambda);
                        
            // Evaluate the Residual
            structure->EvalResidual(input->Get_RigidCriteria());

            if(iIter == 0){initResNorm = structure->Residual.norm();}
            std::cout.width(17); std::cout << log10(structure->Residual.norm() / initResNorm);

            /*--------------------------------------------------
             *   Assembly Ktang, Solve System
             *----------------------------------------------------*/
            
            // Reassembling Stiffness Matrix + Applying Boundary Conditions
            structure->AssemblyTang(iIter);
            
            // Solve Linear System   K*dU = Res = Fext - Fin
            if (nRBE2 != 0 and input->Get_RigidCriteria() == 0) {
              std::cout << "-->  Update KRBE matrix "  << std::endl;
              structure->AssemblyRigidConstr();
              structure->SolveLinearStaticSystem_RBE2(iIter);
            }
            else if (nRBE2 != 0 and input->Get_RigidCriteria() == 1) {
              std::cout << "-->  Update penalty matrix for RBEs "  << std::endl;
              structure->AssemblyRigidPenalty(input->GetPenalty());
              structure->SolveLinearStaticSystem_RBE2_penalty(iIter);
            }
            else {
              structure->SolveLinearStaticSystem(iIter);
             }

            if(iIter == 0){initDispNorm = structure->dU.norm();}
            std::cout.width(17); std::cout << log10(structure->dU.norm() / initDispNorm);

            /*--------------------------------------------------
             *   Updates Coordinates, Updates Rotation Matrices
             *----------------------------------------------------*/
            
            structure->UpdateCoord();
            
            // Now only X is updated
            //structure->UpdateRotationMatrix();  // based on the rotational displacements
            structure->UpdateRotationMatrix_FP();  // based on the rotational displacements            
            structure->UpdateLength();          // Updating length, important
                     
            /*--------------------------------------------------
             *   Update Internal Forces
             *----------------------------------------------------*/
            // Now, X, R, l are updated
            //structure->UpdateInternalForces();
            structure-> UpdateInternalForces_FP();            
            /*--------------------------------------------------
             *   Update Penalty Forces
             *----------------------------------------------------*/
            
            if (nRBE2 != 0 and input->Get_RigidCriteria() == 1)
            {
            //structure->UpdateAxvector_RBE2();
            //structure->EvalPenaltyForces(input->GetPenalty());
            //structure->UpdateRigidConstr(iIter); 
            }
            /*--------------------------------------------------
             *    Check Convergence
             *----------------------------------------------------*/
            
            addouble disp_factor =   structure->dU.norm()/TotalLength;

            UpdateDisplacements();

            std::cout.width(17); std::cout << log10(disp_factor);
            std::cout << std::endl;
            
            if (disp_factor <= input->Get_ConvCriteria()) {
                converged = true;
                totalIter += iIter;
                break;
            }
        }
        
    }
    
    std::cout << "===========================================================================" << std::endl;
    std::cout << std::endl << "--> Exiting Iterative Sequence." << std::endl;
    
    if (input->Get_WriteRestartFlag() ==1)
    {
        WriteRestart();
    }
    
    
}

void CBeamSolver::Restart(int FSIIter = 0){
    
    // Beam total length
    addouble TotalLength = 0;
    for  ( unsigned long iFEM = 0; iFEM < nFEM; iFEM++) {
        TotalLength += element[iFEM]->GetInitial_Length();
    }
    
    std::cout << "--> Setting External Forces" << std::endl;
    structure->ReadForces(nTotalDOF, loadVector);
    
    if (nRBE2 != 0){ 
        std::cout << "--> Setting RBE2 Matrix for Rigid Constraints" << std::endl;
        structure->AddRBE2(input, RBE2);        
    };
    
    addouble  lambda = 1.0;
    addouble dlambda =  1.0/input->Get_LoadSteps() ;
    addouble initResNorm   =  1.0;
    addouble initDispNorm  =  1.0;
    unsigned long iIter;
    unsigned long totalIter = 0;
    unsigned long loadStep = 1;
    cout.setf(ios::fixed, ios::floatfield);

    // This function set the current initial coordinates and memorizes them as the old one before the converging procedure starts
    structure->InitialCoord();
    
    ReadRestart();

    std::cout << "--> Initializing from restart file" << std::endl; 
    structure-> RestartCoord();
    structure-> UpdateLength();
    structure-> InitializeInternalForces();
    
        std::cout << "--> Starting Restart Sequence" << std::endl; 
    std::cout << "===========================================================================" << std::endl;
    
    std::cout.width(6); std::cout << "Iter";
    std::cout.width(17); std::cout << "Log10(Norm_Res)";
    std::cout.width(17); std::cout << "Log10(Lin_Sol)";
    std::cout.width(17); std::cout << "Log10(Norm_Disp)";
    std::cout.width(17); std::cout << "Log10(Disp_Fact)" << std::endl;
    
    int nIter = 1; //  input->Get_nIter()
    
    //===============================================
    //               ITERATIVE SEQUENCE
    //===============================================
    bool converged = false;
    
    for (iIter = 0; iIter < nIter; iIter++) {
        
        std::cout.width(6); std::cout << iIter;
        
        //std::cout << "   ----- ITERATION  -----" << iIter << std::endl;
        
        /*--------------------------------------------------
         *   Updates  Fext, Residual,
         *----------------------------------------------------*/
        
        // Update the External Forces with the loadStep
        structure->UpdateExtForces(lambda);
        
        // Evaluate the Residual
        structure->EvalResidual(input->Get_RigidCriteria());
        
        if(iIter == 0){initResNorm = structure->Residual.norm();}
        std::cout.width(17); std::cout << log10(structure->Residual.norm() / initResNorm);
        
        /*--------------------------------------------------
         *   Assembly Ktang, Solve System
         *----------------------------------------------------*/
        
        // Reassembling Stiffness Matrix + Applying Boundary Conditions
        structure->AssemblyTang(iIter);
        
        // Solve Linear System   K*dU = Res = Fext - Fin
        if (nRBE2 != 0 and input->Get_RigidCriteria() == 0) {
            std::cout << "-->  Update KRBE matrix "  << std::endl;
            structure->AssemblyRigidConstr();
            structure->SolveLinearStaticSystem_RBE2(iIter);
        }
        else if (nRBE2 != 0 and input->Get_RigidCriteria() == 1) {
            std::cout << "-->  Update penalty matrix for RBEs "  << std::endl;
            structure->AssemblyRigidPenalty(input->GetPenalty());
            structure->SolveLinearStaticSystem_RBE2_penalty(iIter);
        }
        else {
            structure->SolveLinearStaticSystem(iIter);
        }
        
        if(iIter == 0){initDispNorm = structure->dU.norm();}
        std::cout.width(17); std::cout << log10(structure->dU.norm() / initDispNorm);
        
        /*--------------------------------------------------
         *   Updates Coordinates, Updates Rotation Matrices
         *----------------------------------------------------*/
        
        structure->UpdateCoord();
        
        // Now only X is updated
        structure->UpdateRotationMatrix();  // based on the rotational displacements
        structure->UpdateLength();          // Updating length, important
        
        /*--------------------------------------------------
         *   Update Internal Forces
         *----------------------------------------------------*/
        // Now, X, R, l are updated
        structure->UpdateInternalForces();
        
        /*--------------------------------------------------
         *   Update Penalty Forces
         *----------------------------------------------------*/
        
        if (nRBE2 != 0 and input->Get_RigidCriteria() == 1)
        {
            //structure->UpdateAxvector_RBE2();
            //structure->EvalPenaltyForces(input->GetPenalty());
            //structure->UpdateRigidConstr(iIter); 
        }
        /*--------------------------------------------------
         *    Check Convergence
         *----------------------------------------------------*/
        
        addouble disp_factor =   structure->dU.norm()/TotalLength;

        UpdateDisplacements();
        
        std::cout.width(17); std::cout << log10(disp_factor);
        std::cout << std::endl;
        
        if (disp_factor <= input->Get_ConvCriteria()) {
            converged = true;
            totalIter += iIter;
            break;
        }
    }    
    
    std::cout << "===========================================================================" << std::endl;
    std::cout << std::endl << "--> Exiting Restart Sequence." << std::endl;    
    
}

passivedouble CBeamSolver::OF_NodeDisplacement(int iNode){
    
    addouble pos1, pos2, pos3;
    
    pos1 = structure->GetDisplacement(iNode, 0);
    pos2 = structure->GetDisplacement(iNode, 1);
    pos3 = structure->GetDisplacement(iNode, 2);
    
    objective_function = sqrt(pow(pos1, 2.0) + pow(pos2, 2.0) + pow(pos3, 2.0));
    
    std::cout.width(20); std::cout << "objective_function = " << objective_function << std::endl;
    
    return AD::GetValue(objective_function);
    
}

void CBeamSolver::SetDependencies(void){

    addouble E, Nu, G;
    unsigned long iFEM, iLoad;

    input->RegisterInput_E();
    input->RegisterInput_Nu();

    E = input->GetYoungModulus();
    Nu = input->GetPoisson();
    G = E/(2*(1+Nu));

    input->SetShear(G);

    for (iFEM = 0; iFEM < nFEM; iFEM++) {
        element[iFEM]->SetDependencies();
    }

    /*--- Initialize vector to store the gradient ---*/
    loadGradient = new passivedouble[nTotalDOF];

    for (iLoad = 0; iLoad < nTotalDOF; iLoad++){

        /*--- Register the load vector ---*/
        AD::RegisterInput(loadVector[iLoad]);

        /*--- Initialize the load gradient to 0 ---*/
        loadGradient[iLoad] = 0.0;
    }

}

void CBeamSolver::ComputeAdjoint(void){
    
    AD::SetDerivative(objective_function, 1.0);

    unsigned long iNode, iLoad;
    unsigned short iDim;
    for (iNode = 0; iNode <  input->Get_nNodes(); iNode++){
      for (iDim =0; iDim < 3; iDim++){
       structure->SetDisplacementAdjoint(iNode, iDim);
      }
    }

    AD::ComputeAdjoint();

    E_grad = input->GetGradient_E();
    Nu_grad = input->GetGradient_Nu();

    for (iLoad = 0; iLoad < nTotalDOF; iLoad++){
        loadGradient[iLoad] = AD::GetValue(AD::GetDerivative(loadVector[iLoad]));
    }

}

void CBeamSolver::UpdateDisplacements(void){

    unsigned long iNode;
    unsigned short iDim;
    for (iNode = 0; iNode <  input->Get_nNodes(); iNode++){
      for (iDim =0; iDim < 3; iDim++){
         structure->SetDisplacement(iNode, iDim);
      }
    }

}


void CBeamSolver::StopRecording(void) {

 AD::RegisterOutput(objective_function);

  unsigned long iNode;
  unsigned short iDim;
  for (iNode = 0; iNode <  input->Get_nNodes(); iNode++){
    for (iDim =0; iDim < 3; iDim++){
       structure->RegisterDisplacement(iNode, iDim);
    }
  }

 AD::StopRecording();

}

void CBeamSolver::StoreDisplacementAdjoint(int iNode, int iDim, passivedouble val_adj){
   structure->StoreDisplacementAdjoint(iNode, iDim, val_adj);
}

void CBeamSolver::WriteRestart(){
    std::ofstream myfile;
    myfile.open ("restart_structure.dat");
        //==== Writing Nodes info
    myfile << "Node ID          ";
    myfile << "Coord. X         ";
    myfile << "Coord. Y         ";
    myfile << "Coord. Z        \n";     
    int posX = 1;    // current  position in the X array
    for (int id_node=1; id_node<= input->Get_nNodes() ; id_node++)   {
        myfile << id_node  << "           ";
        for (int iDim=0; iDim < 3; iDim++) {
            myfile << std::scientific << setprecision(17) << structure->X(posX+iDim-1) << "           " ;
        }
        myfile << "\n";
        posX += 3;
    }    
    //==== Writing Elements info   
    myfile << "Element ID              "; myfile << "Strain 1                "; myfile << "Strain 2                "; myfile << "Strain 3                ";
    myfile << "Strain 4                "; myfile << "Strain 5                "; myfile << "Strain 6            //\n";    
    myfile << "Element ID              "; myfile << "Reference e1(1)         "; myfile << "Reference e1(2)         "; myfile << "Reference e1(3)        ";
    myfile << "Reference e2(1)        "; myfile << "Reference e2(2)        "; myfile << "Reference e2(3)        ";
    myfile << "Reference e3(1)        "; myfile << "Reference e3(2)        "; myfile << "Reference e3(3)    //\n";        
    for (int id_fe=1;     id_fe <= input->Get_nFEM() ; id_fe++){
        myfile << id_fe  << "           ";   
        //strains
        for (int i=1;     i <= 6 ; i++){
            myfile << std::scientific << setprecision(17) << element[id_fe-1]->eps(i-1) << "           " ;            
        }
        myfile << "\n";
        //e1
        myfile << id_fe  << "           ";
        for (int j=1;  j <= 3 ; j++){
            myfile << std::scientific << setprecision(17) << element[id_fe-1]->R(j-1,1-1) << "           " ;                        
        }
        //e2
        for (int j=1;  j <= 3 ; j++){
            myfile << std::scientific << setprecision(17) << element[id_fe-1]->R(j-1,2-1) << "           " ;                        
        }
        //e3
        for (int j=1;  j <= 3 ; j++){
            myfile << std::scientific << setprecision(17) << element[id_fe-1]->R(j-1,3-1) << "           " ;                        
        }        
        myfile << "\n";   
    }
    
    
    myfile.close();
}


void CBeamSolver::ReadRestart(){
    int nNode; double x;double y;double z;
    int nElem; double eps1; double eps2; double eps3; double eps4; double eps5; double eps6; double e11; double e12; double e13; double e21; double e22; double e23; double e31; double e32; double e33;
    int posX = 1;    // current  position in the X array
    string line;
    Vector3dDiff e1; Vector3dDiff e2; Vector3dDiff e3;
    
    ifstream myfile ("restart_structure.dat");
    if (myfile.is_open()){
        getline (myfile,line); //Line of comments for Nodes
        for (int id_node=1; id_node<= input->Get_nNodes() ; id_node++)   {
            getline (myfile,line);   
            CoordExtract( line ,  nNode,  x, y, z);
            structure->X(posX+0-1) = x; structure->X(posX+1-1) = y; structure->X(posX+2-1) = z; 
            posX += 3;
            
        }
        getline (myfile,line); //Line of comments for Elements (1)
        getline (myfile,line); //Line of comments for Elements (2)      
        for (int id_fe=1;     id_fe <= input->Get_nFEM() ; id_fe++)
        {
            getline (myfile,line);   // line of the strain  
            ElemStrainExtract( line , nElem, eps1,eps2, eps3, eps4, eps5, eps6);
            element[id_fe-1]->eps(1-1) = eps1; element[id_fe-1]->eps(2-1) = eps2; element[id_fe-1]->eps(3-1) = eps3;
            element[id_fe-1]->eps(4-1) = eps4; element[id_fe-1]->eps(5-1) = eps5; element[id_fe-1]->eps(6-1) = eps6;     
            
            getline (myfile,line);   // line of the ref system 
            ElemRefExtract( line , nElem, e11, e12, e13, e21, e22, e23, e31, e32, e33);
            e1(1-1) = e11;  e1(2-1) = e12; e1(3-1) = e13;
            e2(1-1) = e21;  e2(2-1) = e22; e2(3-1) = e23;
            e3(1-1) = e31;  e3(2-1) = e32; e3(3-1) = e33;
            
            element[id_fe-1]->R.block(1-1,1-1,3,1) = e1.segment(1-1,3);
            element[id_fe-1]->R.block(1-1,2-1,3,1) = e2.segment(1-1,3);
            element[id_fe-1]->R.block(1-1,3-1,3,1) = e3.segment(1-1,3);
            
            element[id_fe-1]->R.block(4-1,4-1,3,3) = element[id_fe-1]->R.block(1-1,1-1,3,3);   
     
        }
    }
}

void CBeamSolver::CoordExtract(std::string line , int &nNode, double &x,double &y,double &z)
{
        std::istringstream is( line );
        is >> nNode >> x >> y >> z;
}


void CBeamSolver::Debug_Print(int iElement){
    /*
    std::cout << "For element " <<  iElement  << ":" << std::endl;
    std::cout << "Auxiliary vector    = \n" <<  element[iElement]->aux_vector  << std::endl;        
    std::cout << "Kprim    = \n" <<  element[iElement]->Kprim  << std::endl;    
    std::cout << "Length   = " <<  element[iElement]->GetInitial_Length()  << std::endl;
    std::cout << "Rotation matrix = " <<  element[iElement]->Rprev  << std::endl;    
    std::cout << "Property: A   = " <<  element[iElement]->property->GetA()  << std::endl;
    std::cout << "Property: Iyy = " <<  element[iElement]->property->GetIyy()  << std::endl;
    std::cout << "Property: Izz = " <<  element[iElement]->property->GetIzz()  << std::endl;
    std::cout << "Property: Jt  = " <<  element[iElement]->property->GetJt()  << std::endl;  
    std::cout << "Input: E      = " <<  input->GetYoungModulus()  << std::endl;      
    std::cout << "Input: ni     = " <<  input->GetPoisson()  << std::endl;    
    std::cout << "Input: G      = " <<  input->GetShear()  << std::endl;  
    std::cout << "Stiffness matrix:       = \n" <<  structure->Ksys.block(0,0,12,12)  << std::endl;     
    */
}
void CBeamSolver::ElemStrainExtract(std::string line , int &nElem, double &eps1, double &eps2, double &eps3, double &eps4, double &eps5, double &eps6)
{
        std::istringstream is( line );
        is >> nElem >> eps1 >> eps2 >> eps3 >> eps4 >> eps5 >> eps6 ;
}

void CBeamSolver::ElemRefExtract(std::string line , int &nElem, double &e11, double &e12, double &e13, double &e21, double &e22, double &e23, double &e31, double &e32, double &e33)
{
        std::istringstream is( line );
        is >> nElem >> e11 >> e12 >> e13 >> e21 >> e22 >> e23 >> e31 >> e32 >> e33 ;

}
