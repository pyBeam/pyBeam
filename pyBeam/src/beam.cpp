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

    totalIter = 0;
    initResNorm   =  1.0;
    initDispNorm  =  1.0;
    
}

CBeamSolver::~CBeamSolver(void) {
    
}

void CBeamSolver::InitializeInput(CInput* py_input){   // insert node class and connectivity
    
    // I'm memorizing as a member variable the object input passed from outside
    input = py_input;

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
    unsigned long iIter;

    unsigned long loadStep = 1;
    cout.setf(ios::fixed, ios::floatfield);
    
    // This function set the current initial coordinates and memorizes them as the old one before the converging procedure starts
    structure->InitialCoord();

    totalIter = 0;
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

            std::cout.width(17); std::cout << log10(disp_factor);
            std::cout << std::endl;
            
            if (disp_factor <= input->Get_ConvCriteria()) {
                converged = true;
                break;
            }
            totalIter++;
        }
        
    }
    
    std::cout << "===========================================================================" << std::endl;
    std::cout << std::endl << "--> Exiting Iterative Sequence." << std::endl;
       
    if (input->Get_WriteRestartFlag() ==1)
    {
        WriteRestart();
    }
    
    
}

void CBeamSolver::RunRestart(int FSIIter = 0){

    // This function set the current initial coordinates and memorizes them as the old one before the converging procedure starts

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


    /*--- Restart the internal forces ---*/
    std::cout << "--> Initializing from restart file" << std::endl;
    structure->InitialCoord();
    structure->RestartCoord();
    structure->UpdateLength();
    structure->UpdateRotationMatrix_FP();  // based on the rotational displacements
    structure->UpdateInternalForces_FP();

    std::cout << "--> Starting Restart Sequence" << std::endl;
    std::cout << "===========================================================================" << std::endl;

    std::cout.width(8); std::cout << "Iter";
    std::cout.width(16); std::cout << "Log10(Res)";
    std::cout.width(17); std::cout << "Log10(Lin_Sol)";
    std::cout.width(16); std::cout << "Log10(Disp)";
    std::cout.width(17); std::cout << "Log10(Disp_Fact)" << std::endl;

    int nIter = 1; //  input->Get_nIter()

    //===============================================
    //               ITERATIVE SEQUENCE
    //===============================================
    bool converged = false;

    std::cout.width(8); std::cout << "RESTART";

    /*--------------------------------------------------
    *   Updates  Fext, Residual,
    *----------------------------------------------------*/

    // Update the External Forces with the loadStep
    structure->UpdateExtForces(1);

    // Evaluate the Residual
    structure->EvalResidual(input->Get_RigidCriteria());

    std::cout.width(16); std::cout << log10(structure->Residual.norm());

    /*--------------------------------------------------
    *   Assembly Ktang, Solve System
    *----------------------------------------------------*/

    // Reassembling Stiffness Matrix + Applying Boundary Conditions
    structure->AssemblyTang(1);

    // Solve Linear System   K*dU = Res = Fext - Fin
    if (nRBE2 != 0 and input->Get_RigidCriteria() == 0) {
        std::cout << "-->  Update KRBE matrix "  << std::endl;
        structure->AssemblyRigidConstr();
        structure->SolveLinearStaticSystem_RBE2(1);
    }
    else if (nRBE2 != 0 and input->Get_RigidCriteria() == 1) {
        std::cout << "-->  Update penalty matrix for RBEs "  << std::endl;
        structure->AssemblyRigidPenalty(input->GetPenalty());
        structure->SolveLinearStaticSystem_RBE2_penalty(1);
    }
    else {
        structure->SolveLinearStaticSystem(1);
    }

    std::cout.width(16); std::cout << log10(structure->dU.norm());

    /*--------------------------------------------------
    *   Updates Coordinates, Updates Rotation Matrices
    *----------------------------------------------------*/

    structure->UpdateCoord();

    /*--------------------------------------------------
    *    Check Convergence
    *----------------------------------------------------*/

    addouble disp_factor =   structure->dU.norm()/TotalLength;

    std::cout.width(17); std::cout << log10(disp_factor);
    std::cout << std::endl;


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

    /** Register the solution as input **/
    structure->RegisterSolutionInput();

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

    unsigned long iLoad;

    for (unsigned short iTer = 0; iTer < input->Get_nIter(); iTer++){
    
        AD::SetDerivative(objective_function, 1.0);

        structure->SetSolutionAdjoint();

        AD::ComputeAdjoint();

        structure->ExtractSolutionAdjoint();

        E_grad = input->GetGradient_E();
        Nu_grad = input->GetGradient_Nu();

        for (iLoad = 0; iLoad < nTotalDOF; iLoad++){
            loadGradient[iLoad] = AD::GetValue(AD::GetDerivative(loadVector[iLoad]));
        }

        AD::ClearAdjoints();

   }

}

void CBeamSolver::StopRecording(void) {

  AD::RegisterOutput(objective_function);

  /** Register the solution as output **/
  structure->RegisterSolutionOutput();

  AD::StopRecording();

}

void CBeamSolver::StoreDisplacementAdjoint(int iNode, int iDim, passivedouble val_adj){
   structure->StoreDisplacementAdjoint(iNode, iDim, val_adj);
}

void CBeamSolver::WriteRestart(){
    std::ofstream myfile;
    myfile.open ("restart_structure.dat");
    int posX = 1;
    //==== Writing Nodes info
    myfile << "Node ID              "; myfile << "U 1                "; myfile << "U 2                "; myfile << "U 3                ";
    myfile << "U 4                "; myfile << "U 5                "; myfile << "U 6            //\n";          
    for (int id_node=1; id_node<= input->Get_nNodes() ; id_node++)   {
        myfile << id_node  << "           ";
        for (int iDim=0; iDim < 6; iDim++) {
            myfile << std::scientific << setprecision(17) << structure->U(posX+iDim-1) << "           " ;
        }
        myfile << "\n";
        posX += 6;
    }
    
    myfile.close();
}


void CBeamSolver::ReadRestart(){
    int nNode; double Ux;double Uy;double Uz; double Urx;double Ury;double Urz;
    //int nElem; double eps1; double eps2; double eps3; double eps4; double eps5; double eps6; double e11; double e12; double e13; double e21; double e22; double e23; double e31; double e32; double e33;
    int posX = 1;    // current  position in the X array
    string line;


    ifstream myfile ("restart_structure.dat");
    if (myfile.is_open()){
        getline (myfile,line); //Line of comments for Nodes
        for (int id_node=1; id_node<= input->Get_nNodes() ; id_node++)   {
            getline (myfile,line);   
            UExtract( line , nNode, Ux, Uy, Uz, Urx, Ury, Urz);
            structure->U(posX+0-1) = Ux; structure->U(posX+1-1) = Uy; structure->U(posX+2-1) = Uz;
            structure->U(posX+3-1) = Urx; structure->U(posX+4-1) = Ury; structure->U(posX+5-1) = Urz;
            posX += 6;
            
        }
    } 
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

void CBeamSolver::UExtract(std::string line ,int &nNode, double &Ux, double &Uy, double &Uz, double &Urx, double &Ury, double &Urz)
{
        std::istringstream is( line );
        is >> nNode >> Ux >> Uy >> Uz >> Urx >> Ury >> Urz ;
}

