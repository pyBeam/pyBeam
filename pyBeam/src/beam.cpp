/*
 * pyBeam, an open-source Beam Solver
 *
 * Copyright (C) 2019 by the authors
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
    
    if (verbose){cout << "--> Node Structure Initialization... ";}
    unsigned long nNodes = input->Get_nNodes();   //substitute from mesh file
    node = new CNode*[nNodes];
    if (verbose){std::cout << "ok!"  <<std::endl;}
    
    //==============================================================
    //      Finite Element vector initialization
    //==============================================================
    
    if (verbose){cout << "--> Finite Element Initialization... ";}
    element = new CElement*[nFEM];
    if (verbose){std::cout << "ok!"  <<std::endl;}
    
    //==============================================================
    //      RBE2 initialization
    //==============================================================
    if (nRBE2 != 0){
        if (verbose){cout << "--> RBE2 initialization... ";}
        RBE2 = new CRBE2*[nRBE2];
        if (verbose){std::cout << "ok!"  <<std::endl;}
        if (verbose){std::cout << "--> Warning: the code works if slave nodes are not connected to beams! "
                               << std::endl;}
        if (verbose){std::cout << "--> Remember to add a check about this..."  <<std::endl;}
    }
    else {RBE2 = NULL; }
    
    
    //===============================================
    //  Initialize structural solver
    //===============================================
    structure = NULL;
    
}

/**
 * @todo Remove unused parameters
 * @body Remove unused FSIIter from the arguments of the function, and also from the python wrapper and library.
 */


void CBeamSolver::Solve(int FSIIter = 0){
    
    std::ofstream history;
    history.open ("History.pyBeam");
    
    // Beam total length
    addouble TotalLength = 0;
    for  ( unsigned long iFEM = 0; iFEM < nFEM; iFEM++) {
        TotalLength += element[iFEM]->GetInitial_Length();
    }

    if (verbose){std::cout << "--> Setting External Forces" << std::endl;}
    structure->ReadForces(nTotalDOF, loadVector);
    
    if (nRBE2 != 0){
        structure->AddRBE2(input, RBE2);
    }

    //===============================================
    // LOAD STEPPING
    //===============================================
    
    if (verbose){std::cout << "--> Starting Load Stepping" << std::endl << std::endl;}

    addouble  lambda = 1.0;
    addouble dlambda =  1.0/input->Get_LoadSteps() ;
    unsigned long iIter;

    unsigned long loadStep = 1;
    cout.setf(ios::fixed, ios::floatfield);
    history.setf(ios::fixed, ios::floatfield);
    
    // This function set the current initial coordinates and memorizes them as the old one 
    // before the converging procedure starts (necessary in case of FSI)
    structure->InitialCoord();
    structure->RestartCoord();
    structure->UpdateRotationMatrix_FP();  // based on the rotational displacements
    structure->UpdateLength();
    structure->UpdateInternalForces_FP();
    
    totalIter = 0;
    for  ( loadStep = 0; loadStep < input->Get_LoadSteps(); loadStep++) {
        
        lambda = dlambda*(loadStep+1);
        if (verbose){
            std::cout << "===========================================================================" << std::endl;
            std::cout << "==       LOAD STEP     " << loadStep << std::endl;
            std::cout.precision(8);
            std::cout << "==       Lambda        " << lambda << std::endl;
            std::cout << "===========================================================================" << std::endl;
            
            std::cout.width(6); std::cout << "Iter";
            std::cout.width(17); std::cout << "Log10(Norm_Res)";
            if (nRBE2 != 0 ) {
                std::cout.width(17); std::cout << " Log10(Norm[Res+Constr])";    
            }
            std::cout.width(17); std::cout << "Log10(Lin_Sol)";
            std::cout.width(17); std::cout << "Log10(Norm_Disp)";
            std::cout.width(17); std::cout << "Log10(Disp_Fact)" << std::endl;
            
            // WRITE HISTORY FILE
            history << "===========================================================================" << std::endl;
            history << "==       LOAD STEP     " << loadStep << std::endl ;
            history.precision(8);
            history << "==       Lambda        " << lambda << std::endl ;
            history << "===========================================================================" << std::endl;
            
            history.width(6); history << "Iter";
            history.width(17); history << "Log10(Norm_Res)";
            if (nRBE2 != 0 ) {
                history.width(17); history << " Log10(Norm[Res+Constr])";    
            }
            history.width(17); history << "Log10(Lin_Sol)";
            history.width(17); history << "Log10(Norm_Disp)";
            history.width(17); history << "Log10(Disp_Fact)" << std::endl;        
        }
        
        //===============================================
        //               ITERATIVE SEQUENCE
        //===============================================
        bool converged = false;
        
        for (iIter = 0; iIter < input->Get_nIter(); iIter++) {

            if (verbose){std::cout.width(6); std::cout << iIter;
                         history.width(6); history << iIter;}

            //std::cout << "   ----- ITERATION  -----" << iIter << std::endl;
            
            /*--------------------------------------------------
             *   Updates  Fext, Residual,
             *----------------------------------------------------*/
            
            // Update the External Forces with the loadStep
            structure->UpdateExtForces(lambda);

            // Evaluate the Residual
            structure->EvalResidual();

            
            if (verbose){
                std::cout.width(17); std::cout << log10(structure->Residual.norm());
                history.width(17); history << log10(structure->Residual.norm());
            }

            /*--------------------------------------------------
             *   Assembly Ktang, Solve System
             *----------------------------------------------------*/
            
            // Reassembling Stiffness Matrix + Applying Boundary Conditions
            structure->AssemblyTang(1);
            
            
            if (nRBE2 != 0 ) {
                if (iRigid == 0 ) {
                structure->SetPenalty();
                structure->RigidResidual();
                //structure->RigidResidual_FD();
                structure->AssemblyRigidPenalty();
                //structure->AssemblyRigidPenalty_FD();
                }
                else{
                structure->SetRigidLagrangeDimensions(); 
                structure->RigidResidualLagrange();
                structure->AssemblyRigidLagrange();
                }
                if (verbose){
                    std::cout.width(17); std::cout << log10(structure->Residual.norm());
                    history.width(17); history << log10(structure->Residual.norm());
                }                
            }
            
            if (nRBE2 != 0 and iRigid == 1){
            structure->ImposeBC_RigidLagrangian(); 
            // Solve Linear System   Ksys_lam*dU_lam = Res_lam 
            structure->SolveLinearStaticSystem_RigidLagrangian(iIter,history,1);                
            }
            else{
            structure->ImposeBC(); 
            // Solve Linear System   Ksys*dU = Res =
            structure->SolveLinearStaticSystem(iIter,history,1);
            }

            if (verbose){
                std::cout.width(17); std::cout << log10(structure->dU.norm());
                history.width(17); history << log10(structure->dU.norm());
            }

            /*--------------------------------------------------
             *   Updates Coordinates, Updates Rotation Matrices
             *----------------------------------------------------*/
            
            structure->UpdateCoord(nRBE2,iRigid);
            
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
             *    Check Convergence
             *----------------------------------------------------*/
            
            addouble disp_factor =   structure->dU.norm()/TotalLength;
            
            if (verbose){
                std::cout.width(17); std::cout << log10(disp_factor);
                history.width(17); history << log10(disp_factor);
            }
            if (verbose){std::cout << std::endl;  history << std::endl;}
            
            if (disp_factor <= input->Get_ConvCriteria()) {
                converged = true;
                break;
            }
            totalIter++;
        }
        
    }
    
    if (verbose){
        std::cout << "===========================================================================" << std::endl;
        history << "==========================================================================="<< std::endl;
        std::cout << std::endl << "--> Writing Restart file (restart.pyBeam)." << std::endl;
    }
    WriteRestart();
    if (verbose){std::cout << std::endl << "--> Exiting Iterative Sequence." << std::endl;}
    
    history.close();
    
}

void CBeamSolver::RunRestart(int FSIIter = 0){
    
    std::ofstream history;
    
    // This function set the current initial coordinates and memorizes them as the old one before the converging procedure starts
    
    // Beam total length
    addouble TotalLength = 0;
    for  ( unsigned long iFEM = 0; iFEM < nFEM; iFEM++) {
        TotalLength += element[iFEM]->GetInitial_Length();
    }
    
    if (verbose){std::cout << "--> Setting External Forces" << std::endl;}
    structure->ReadForces(nTotalDOF, loadVector);
    
    if (nRBE2 != 0){
        if (verbose){std::cout << "--> Setting RBE2 Matrix for Rigid Constraints" << std::endl;}
        structure->AddRBE2(input, RBE2);
    }
    
    /*--- Restart the internal forces ---*/
    if (verbose){std::cout << "--> Initializing from restart file" << std::endl;}
    structure->InitialCoord();
    structure->RestartCoord();
    structure->UpdateLength();
    structure->UpdateRotationMatrix_FP();  // based on the rotational displacements
    structure->UpdateInternalForces_FP();
    
    if (verbose){
        std::cout << "--> Starting Restart Sequence" << std::endl;
        std::cout << "===========================================================================" << std::endl;
        
        cout.setf(ios::fixed, ios::floatfield);
        history.setf(ios::fixed, ios::floatfield);
        
        std::cout.width(8); std::cout << "Iter";
        std::cout.width(16); std::cout << "Log10(Res)";
        std::cout.width(17); std::cout << "Log10(Lin_Sol)";
        std::cout.width(16); std::cout << "Log10(Disp)";
        std::cout.width(17); std::cout << "Log10(Disp_Fact)" << std::endl;
    }
    
    //===============================================
    //               RESTART SEQUENCE
    //===============================================
    
    if (verbose){std::cout.width(8); std::cout << "RESTART";}
    
    /*--------------------------------------------------
     *   Updates  Fext, Residual,
     *----------------------------------------------------*/
    
    // Update the External Forces with the loadStep
    structure->UpdateExtForces(1);
    
    // Evaluate the Residual
    structure->EvalResidual();
    
    if (verbose){std::cout.width(16); std::cout << log10(structure->Residual.norm());}
    
    /*--------------------------------------------------
     *   Assembly Ktang, Solve System
     *----------------------------------------------------*/
    
    // Reassembling Stiffness Matrix + Applying Boundary Conditions
    structure->AssemblyTang(1);
    
    if (nRBE2 != 0 ) {
        if (iRigid == 0 ) {
            structure->SetPenalty();
            structure->RigidResidual();
            //structure->RigidResidual_FD();
            structure->AssemblyRigidPenalty();
            //structure->AssemblyRigidPenalty_FD();
        }
        else{
            structure->SetRigidLagrangeDimensions(); 
            structure->RigidResidualLagrange();
            structure->AssemblyRigidLagrange();
        }
        if (verbose){
            std::cout.width(17); std::cout << log10(structure->Residual.norm());
            history.width(17); history << log10(structure->Residual.norm());
        }                
    }
    
    if (nRBE2 != 0 and iRigid == 1){
        structure->ImposeBC_RigidLagrangian(); 
        // Solve Linear System   Ksys_lam*dU_lam = Res_lam 
        structure->SolveLinearStaticSystem_RigidLagrangian(0,history,1);                
    }
    else{
        structure->ImposeBC(); 
        // Solve Linear System   Ksys*dU = Res =
        structure->SolveLinearStaticSystem(0,history,1);
    }
    
    if (verbose){
        std::cout.width(17); std::cout << log10(structure->dU.norm());
        history.width(17); history << log10(structure->dU.norm());
    }
    
    /*--------------------------------------------------
     *   Updates Coordinates, Updates Rotation Matrices
     *----------------------------------------------------*/
    
    structure->UpdateCoord(nRBE2,iRigid);
    
    
    /*--------------------------------------------------
     *    Check Convergence
     *----------------------------------------------------*/

    addouble disp_factor =   structure->dU.norm()/TotalLength;

    if (verbose){
    std::cout.width(17); std::cout << log10(disp_factor);
    std::cout << std::endl;


    std::cout << "===========================================================================" << std::endl;
    std::cout << std::endl << "--> Exiting Restart Sequence." << std::endl;
    }
    
}

passivedouble CBeamSolver::OF_NodeDisplacement(int iNode){

    addouble pos1, pos2, pos3;

    pos1 = structure->GetDisplacement(iNode, 0);
    pos2 = structure->GetDisplacement(iNode, 1);
    pos3 = structure->GetDisplacement(iNode, 2);
    
    objective_function = sqrt(pow(pos1, 2.0) + pow(pos2, 2.0) + pow(pos3, 2.0));

    return AD::GetValue(objective_function);

}

void CBeamSolver::SetDependencies(void){

    /** Register the solution as input **/
    structure->RegisterSolutionInput();

    addouble E, E_dim, Nu, G;
    unsigned long iFEM, iLoad;

    input->RegisterInput_E();
    input->RegisterInput_Nu();

    E_dim = input->GetYoungModulus_dimensional();
    structure->SetDimensionalYoungModulus(E_dim);

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
    myfile.open ("restart.pyBeam");
    int posX = 1;
    //==== Writing Nodes info
    myfile << "Node ID              ";
    myfile << "U 1                "; myfile << "U 2                ";
    myfile << "U 3                "; myfile << "U 4                ";
    myfile << "U 5                "; myfile << "U 6            //\n";
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
    int posX = 1;    // current  position in the X array
    string line;

    ifstream myfile ("solution.pyBeam");
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
    if (myfile.fail()){
        cout << "Error opening solution file (solution.pyBeam)." << endl;
        exit (EXIT_FAILURE);
    }
}

void CBeamSolver::UExtract(std::string line ,int &nNode, double &Ux, double &Uy,
                           double &Uz, double &Urx, double &Ury, double &Urz)
{
    std::istringstream is( line );
    is >> nNode >> Ux >> Uy >> Uz >> Urx >> Ury >> Urz ;
}

