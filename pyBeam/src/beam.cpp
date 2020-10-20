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
    resp_weight = 0.0;
    resp_KS = 0.0;

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
    nProp=input->Get_nProp();    
    
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
    //      Property vector initialization
    //==============================================================
    
    if (verbose){cout << "--> Property Initialization... ";}
   
    Prop = new CProperty*[nProp];
    
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


void CBeamSolver::InitializeStructure(void) {
    
    structure = new CStructure(input, element, node); structure->SetCoord0();
    // If there are RBE and Lagrange multiplier method is used it's important to define already the dimension as AD recording needs it
    if (nRBE2 != 0){ 
        if (verbose){std::cout << "--> Setting RBE2 for Rigid Constraints" << std::endl;}        
            structure->AddRBE2(input, RBE2);}
            // Sizing rigid support matrix. Necessary in case of SPARSE approach
            // Should be also done here for SPARSE approach and rigid lagraignan (in the future)
        if (nRBE2 != 0 and iRigid == 0){
#ifdef DENSE
            // Resizes and zeros the K matrices
            structure->K_penal.resize(input->Get_nNodes()*6,input->Get_nNodes()*6);
            structure->K_penal = MatrixXdDiff::Zero(input->Get_nNodes()*6,input->Get_nNodes()*6);
#else
            // Resizes the sparse Ksys matriX
            structure->K_penal.resize( (input->Get_nNodes())*6,(input->Get_nNodes())*6);    
            structure->tripletListRBEPenalty.reserve((input->Get_nNodes())*6*100);     
#endif  
    }
    if (nRBE2 != 0 and iRigid == 1){
    structure->SetRigidLagrangeDimensions();}
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
    structure->InitialCoord(); //sets the Coordinates 0 of the configuration (as given in the mesh)
    structure->RestartCoord(); //sets the current Coordinates of the configuration (as given in restart or from the previous iter of FSI)
    structure->UpdateLength();
//    structure->UpdateRotationMatrix_FP();  // based on the rotational displacements
//    structure->UpdateInternalForces_FP();
       
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
                 if (iIter == 0 ) {structure->SetPenalty();}
                structure->RigidResidual();
                structure->AssemblyRigidPenalty();
                }
                else{
                
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
    
    // Resetting Fnom in case Solve is called again
    ResetLoads(); 
    
}



void CBeamSolver::SolveLin(int FSIIter = 0){
    
    std::ofstream history;
    history.open ("History.pyBeam");
    
    // Beam total length
    addouble TotalLength = 0;
    for  ( unsigned long iFEM = 0; iFEM < nFEM; iFEM++) {
        TotalLength += element[iFEM]->GetInitial_Length();
    }

    if (verbose){std::cout << "--> Setting External Forces" << std::endl;}
    structure->ReadForces(nTotalDOF, loadVector);
    
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
    structure->InitialCoord(); //sets the Coordinates 0 of the configuration (as given in the mesh)
    structure->RestartCoord(); //sets the current Coordinates of the configuration (as given in restart or from the previous iter of FSI)
    // N.B. Internal forces are the one found at the previosu iteration
    structure-> UpdateInternalForcesLinear ();  // For consistency is better to update them
        
        
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
   
        if (verbose){std::cout.width(6); std::cout << iIter;
                     history.width(6); history << iIter;}

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
         *   Assembly Kelastic tangent, Solve System
         *----------------------------------------------------*/
        // Assembling Stiffness Matrix + Applying Boundary Conditions

        structure->AssemblyElasticStiffness();


        if (nRBE2 != 0 ) {
            if (iRigid == 0 ) {
             if (iIter == 0 ) {structure->SetPenalty();}
            structure->RigidResidual();
            structure->AssemblyRigidPenalty();
            }
            else{

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

        structure->UpdateCoordLIN(nRBE2,iRigid);

        structure-> UpdateInternalForcesLinear ();       
        

        if (verbose){std::cout << std::endl;  history << std::endl;}
   
    }
    
    if (verbose){
        std::cout << "===========================================================================" << std::endl;
        history << "==========================================================================="<< std::endl;
        std::cout << std::endl << "--> Writing Restart file (restart.pyBeam)." << std::endl;
    }
    WriteRestart();
    if (verbose){std::cout << std::endl << "--> Exiting Iterative Sequence." << std::endl;}
    
    history.close();
    
    // Resetting Fnom in case Solve is called again
    ResetLoads(); 
    
}



void CBeamSolver::RunRestart(int FSIIter = 0){
 
    std::ofstream history;
    history.open ("History_restart.pyBeam");    
    
    // This function set the current initial coordinates and memorizes them as the old one before the converging procedure starts
    
    // Beam total length
    addouble TotalLength = 0;
    for  ( unsigned long iFEM = 0; iFEM < nFEM; iFEM++) {
        TotalLength += element[iFEM]->GetInitial_Length();
    }
    
    if (verbose){std::cout << "--> Setting External Forces" << std::endl;}
    structure->ReadForces(nTotalDOF, loadVector);
    
    
    /*--- Restart the internal forces ---*/
    if (verbose){std::cout << "--> Initializing from restart file" << std::endl;}
    structure->InitialCoord(); //sets the Coordinates 0 of the configuration (as given in the mesh)
    structure->RestartCoord(); //sets the current Coordinates of the configuration (as given in restart or from the previous iter of FSI)
    structure->UpdateLength();
    structure->UpdateRotationMatrix_FP();  // based on the rotational displacements
    structure->UpdateInternalForces_FP();
    
    if (verbose){
        std::cout << "--> Starting Restart Sequence" << std::endl;
        std::cout << "===========================================================================" << std::endl;
        
        history << "===========================================================================" << std::endl;
        
        cout.setf(ios::fixed, ios::floatfield);
        history.setf(ios::fixed, ios::floatfield);
        
        std::cout.width(8); std::cout << "Iter";
        std::cout.width(16); std::cout << "Log10(Res)";
        if (nRBE2 != 0 ) {
            std::cout.width(17); std::cout << " Log10(Norm[Res+Constr])";    
        }        
        std::cout.width(17); std::cout << "Log10(Lin_Sol)";
        std::cout.width(16); std::cout << "Log10(Disp)";
        std::cout.width(17); std::cout << "Log10(Disp_Fact)" << std::endl;
        
        // WRITE HISTORY FILE
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
    //               RESTART SEQUENCE
    //===============================================
    
    if (verbose){std::cout.width(8); std::cout << "RESTART";
                 history.width(8); history << "RESTART";}
    
    /*--------------------------------------------------
     *   Updates  Fext, Residual,
     *----------------------------------------------------*/
    
    // Update the External Forces with the loadStep
    structure->UpdateExtForces(1);
    
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
            structure->AssemblyRigidPenalty();
        }
        else{
            
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
                history.width(17); history << log10(disp_factor);
                std::cout << std::endl;
                history << std::endl;


        std::cout << "===========================================================================" << std::endl;
        std::cout << std::endl << "--> Exiting Restart Sequence." << std::endl;
        history << "===========================================================================" << std::endl;
        history << std::endl << "--> Exiting Restart Sequence." << std::endl;
    }

    // Reset of load is not used here as restart is usually run before Adjoint which requires the loadvector array


}


void CBeamSolver::RunRestartLin(int FSIIter = 0){
 
    std::ofstream history;
    history.open ("History_restart.pyBeam");    
    
    // This function set the current initial coordinates and memorizes them as the old one before the converging procedure starts
    
    // Beam total length
    addouble TotalLength = 0;
    for  ( unsigned long iFEM = 0; iFEM < nFEM; iFEM++) {
        TotalLength += element[iFEM]->GetInitial_Length();
    }
    
    if (verbose){std::cout << "--> Setting External Forces" << std::endl;}
    structure->ReadForces(nTotalDOF, loadVector);
    
    
    /*--- Restart the internal forces ---*/
    if (verbose){std::cout << "--> Initializing from restart file" << std::endl;}
    structure->InitialCoord(); //sets the Coordinates 0 of the configuration (as given in the mesh)
    structure->RestartCoord(); //sets the current Coordinates of the configuration (as given in restart or from the previous iter of FSI)
    structure-> UpdateInternalForcesLinear ();
     
    
    if (verbose){
        std::cout << "--> Starting Restart Sequence" << std::endl;
        std::cout << "===========================================================================" << std::endl;        
        history << "===========================================================================" << std::endl;
        
        cout.setf(ios::fixed, ios::floatfield);
        history.setf(ios::fixed, ios::floatfield);
        
        std::cout.width(8); std::cout << "Iter";
        std::cout.width(16); std::cout << "Log10(Res)";
        if (nRBE2 != 0 ) {
            std::cout.width(17); std::cout << " Log10(Norm[Res+Constr])";    
        }        
        std::cout.width(17); std::cout << "Log10(Lin_Sol)";
        std::cout.width(16); std::cout << "Log10(Disp)";
        std::cout.width(17); std::cout << "Log10(Disp_Fact)" << std::endl;
        
        // WRITE HISTORY FILE
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
    //               RESTART SEQUENCE
    //===============================================
    
    if (verbose){std::cout.width(8); std::cout << "RESTART";
        history.width(8); history << "RESTART";}
    
    /*--------------------------------------------------
     *   Updates  Fext, Residual,
     *----------------------------------------------------*/
       
    // Update the External Forces with the loadStep
    structure->UpdateExtForces(1);
    
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
    structure->AssemblyElasticStiffness();
    if (nRBE2 != 0 ) {
        if (iRigid == 0 ) {
            structure->SetPenalty();
            structure->RigidResidual();
            structure->AssemblyRigidPenalty();
        }
        else{
            
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
//    std::cout << " ==== RESIDUAL ==== " <<std::endl;
//    std::cout << structure->Residual.segment(60*6,3) <<std::endl;
//    std::cout << " ================== END RESIDUAL " <<std::endl;
//    std::cout << "  ==== FEXT==== " <<std::endl;
//    std::cout << structure->Fext.segment(60*6,3) <<std::endl;
    
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
    
    structure->UpdateCoordLIN(nRBE2,iRigid);
//    structure->U=structure->dU;
//    std:cout << "\n\n NOW WRITING U" << structure->U<<std::endl;
     structure-> UpdateInternalForcesLinear ();
    
    /*--------------------------------------------------
     *    Check Convergence
     *----------------------------------------------------*/

    addouble disp_factor =   structure->dU.norm()/TotalLength;

    if (verbose){
                std::cout.width(17); std::cout << log10(disp_factor);
                history.width(17); history << log10(disp_factor);
                std::cout << std::endl;
                history << std::endl;


        std::cout << "===========================================================================" << std::endl;
        std::cout << std::endl << "--> Exiting Restart Sequence." << std::endl;
        history << "===========================================================================" << std::endl;
        history << std::endl << "--> Exiting Restart Sequence." << std::endl;
    }

    // Reset of load is not used here as restart is usually run before Adjoint which requires the loadvector array


}




passivedouble CBeamSolver::OF_NodeDisplacement(int iNode){

    addouble pos1, pos2, pos3;

    pos1 = structure->GetDisplacement(iNode, 0);
    pos2 = structure->GetDisplacement(iNode, 1);
    pos3 = structure->GetDisplacement(iNode, 2);
    
    objective_function = sqrt(pow(pos1, 2.0) + pow(pos2, 2.0) + pow(pos3, 2.0));

    return AD::GetValue(objective_function);
}

passivedouble CBeamSolver::EvalWeight(){
    resp_weight = structure->EvaluateWeight();
    return AD::GetValue(resp_weight);     };
   
    
passivedouble CBeamSolver::EvalKSStress(){
    resp_KS = structure->Evaluate_no_AdaptiveKSstresses();
    return AD::GetValue(resp_KS);     };
 
    
//////  DEBUG 
passivedouble CBeamSolver::EvalSigmaBoom(){
//    resp_sigmaboom = element[0]->GetSB(); //sigma_booms(0);
//    resp_sigmaboom = element[0]->Getdsigma_dx();   
//    resp_sigmaboom = element[0]->Gettau();
//    resp_sigmaboom = element[0]->Gettau();
    resp_sigmaboom = element[0]->Getg(); 
  
    return AD::GetValue(resp_sigmaboom);     };
    
//////  DEBUG     
passivedouble CBeamSolver::EvalEA(){
    resp_EA = element[0]->GetEA();
    return AD::GetValue(resp_EA);     };
//////  DEBUG     
passivedouble CBeamSolver::EvalIzz_b(){
    resp_Izz_b = element[0]->GetIzz_b();
    return AD::GetValue(resp_Izz_b);     };
//////  DEBUG  
passivedouble CBeamSolver::EvalNint(){
    resp_Nint = element[0]->RetrieveNint();
    return AD::GetValue(resp_Nint);     };

    
   
/** Core function that exposes the dependencies*/    
/** Order d/d  U, E, Nu, Props, A,Loads      
 */
void CBeamSolver::SetDependencies(void){

    /** Register the solution as input **/
    structure->RegisterSolutionInput(iRigid);

    addouble E, E_dim, Nu, G;
    unsigned long iFEM,iRBE2, iLoad, iP;

    input->RegisterInput_E();
    input->RegisterInput_Nu();

    E_dim = input->GetYoungModulus_dimensional();
    structure->SetDimensionalYoungModulus(E_dim);

    E = input->GetYoungModulus();
    Nu = input->GetPoisson();
    G = E/(2*(1+Nu));

    input->SetShear(G);
    
    /* Registering wing-box sizes (or inertias) as inputs.
     Need to be done before exposing dependencies of FE from Properties*/
    
    // Counting total numbers of prop DVs 
    nPropDVs = 0;
    std::vector<int> index; 
    index.push_back(0);
    for (iP= 0; iP<nProp; iP ++){
        //Prop[iP]->RegisterInput_WB();
        if (Prop[iP]->GetisWBDV() == 0){
            nPropDVs += 4; }
        else if (Prop[iP]->GetisWBDV() == 1){
            nPropDVs += 6;}    
        else {std::cout << "======= ERROR  ======" <<  Prop[iP]->GetisWBDV() << std::endl; }
        index.push_back(nPropDVs);
        }
    /* Creating vector storing all WB/Property DVs*/  

    propDVsVector = new addouble[nPropDVs]; 
    
    for (int iPDVs= 0; iPDVs<nPropDVs; iPDVs ++){propDVsVector[iPDVs]=0.0;}         

    for (iP= 0; iP<nProp; iP ++){
        Prop[iP]->InitializePropDVsVec(propDVsVector,index[iP]);  }
 
        /// FOR DEBUGGING
//    Prop[0]->RegisterInput_A();
    
    /* Registering as inputs WB/Property DVs*/    
    for (int iPDV= 0; iPDV<nPropDVs; iPDV ++){
        AD::RegisterInput(propDVsVector[iPDV]);
    }
    

    
    /* Exposing dependencies  DVs*/
    for (iP= 0; iP<nProp; iP ++){
        Prop[iP]->SetDependencyfromDVVec(propDVsVector,index[iP]);         
        }
    

    /*--- Initialize vector to store the gradient wrt to propDV vector ---*/
    propGradient = new passivedouble[nPropDVs]; 
    for (int iPDV= 0; iPDV<nPropDVs; iPDV ++){propGradient[iPDV] = 0.0;}

    /* Exposing dependencies of FEMs*/
    for (iFEM = 0; iFEM < nFEM; iFEM++) {
        element[iFEM]->SetDependencies();
    }
    
    /* Exposing dependencies of nRBE2*/
    if (nRBE2 != 0){
        for (iRBE2 = 0; iRBE2 < nRBE2; iRBE2++) {    
            RBE2[iRBE2]->SetDependencies();    
        }
    }
    
    /*--- Initialize vector to store the gradient wrt to LOADS ---*/
    /*--- and REGISTER LOADS as INPUTS ---*/    
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
        
        
        //In the adjoint system of the field  and objective dual equations sets the
        // initial values of objective adjoint to 1 
        AD::SetDerivative(objective_function, 1.0);
        
        //In the adjoint system of the field  and objective dual equations sets the
        // initial values of state adjoints to the previous iteration values + adds source term
        
        structure->SetSolutionAdjoint(iRigid);  // takes care of the source term and adjoint to disp 
        
        //          
        AD::ComputeAdjoint();
        
        // Extracts new value of the adjoint to U
        structure->ExtractSolutionAdjoint(iRigid);  // extract gradient wrt to state variables
        
        // Extracts the sensitivities wrt to E_grad and Nu_grad 
        E_grad = input->GetGradient_E();
        Nu_grad = input->GetGradient_Nu();
        

        for (int iPDV = 0; iPDV < nPropDVs; iPDV++){
            propGradient[iPDV] = AD::GetValue(AD::GetDerivative(propDVsVector[iPDV]));
        }   
        
        for (iLoad = 0; iLoad < nTotalDOF; iLoad++){
            loadGradient[iLoad] = AD::GetValue(AD::GetDerivative(loadVector[iLoad]));
        }        
        
        AD::ClearAdjoints();
    }
}




/////DEBUG
void CBeamSolver::ComputeAdjointNint(void){

    unsigned long iLoad;

    for (unsigned short iTer = 0; iTer < input->Get_nIter(); iTer++){ 
        
        //In the adjoint system of the field  and objective dual equations sets the
        // initial values of objective adjoint to 1 
        AD::SetDerivative(resp_Nint, 1.0);
        
        //In the adjoint system of the field  and objective dual equations sets the
        // initial values of state adjoints to the previous iteration values + adds source term
        
        structure->SetSolutionAdjoint(iRigid);  // takes care of the source term and adjoint to disp 
        
        //          
        AD::ComputeAdjoint();
        
        // Extracts new value of the adjoint to U
        structure->ExtractSolutionAdjoint(iRigid);  // extract gradient wrt to state variables
        
        // Extracts the sensitivities wrt to E_grad and Nu_grad 
        E_grad = input->GetGradient_E();
        Nu_grad = input->GetGradient_Nu();
        
        for (int iPDV = 0; iPDV < nPropDVs; iPDV++){
            propGradient[iPDV] = AD::GetValue(AD::GetDerivative(propDVsVector[iPDV]));
        }           
        for (iLoad = 0; iLoad < nTotalDOF; iLoad++){
            loadGradient[iLoad] = AD::GetValue(AD::GetDerivative(loadVector[iLoad]));
        }                
        AD::ClearAdjoints();
    }
}


void CBeamSolver::ComputeAdjointWeight(void){

    unsigned long iLoad;

    AD::SetDerivative(resp_weight, 1.0);

//    structure->SetSolutionAdjoint(iRigid);

    AD::ComputeAdjoint();
    structure->ExtractSolutionAdjoint(iRigid);

    E_grad = input->GetGradient_E();
    Nu_grad = input->GetGradient_Nu();
    
//    A_grad = Prop[0]->GetGradient_A();
    
    for (int iPDV = 0; iPDV < nPropDVs; iPDV++){
        propGradient[iPDV] = AD::GetValue(AD::GetDerivative(propDVsVector[iPDV]));
    }      
    
    for (iLoad = 0; iLoad < nTotalDOF; iLoad++){
        loadGradient[iLoad] = AD::GetValue(AD::GetDerivative(loadVector[iLoad]));
    }
 
    AD::ClearAdjoints();    
}



void CBeamSolver::ComputeAdjointKS(void){

    unsigned long iLoad;
    
    for (unsigned short iTer = 0; iTer < input->Get_nIter(); iTer++){ 

        AD::SetDerivative(resp_KS, 1.0);

        structure->SetSolutionAdjoint(iRigid);

        AD::ComputeAdjoint();
        structure->ExtractSolutionAdjoint(iRigid);

        E_grad = input->GetGradient_E();
        Nu_grad = input->GetGradient_Nu();


        for (int iPDV = 0; iPDV < nPropDVs; iPDV++){
            propGradient[iPDV] = AD::GetValue(AD::GetDerivative(propDVsVector[iPDV]));
        }      

        for (iLoad = 0; iLoad < nTotalDOF; iLoad++){
            loadGradient[iLoad] = AD::GetValue(AD::GetDerivative(loadVector[iLoad]));
        }

        AD::ClearAdjoints();  
    }
}


void CBeamSolver::ComputeAdjointSigmaBoom(void){
    
    unsigned long iLoad;
    
    for (unsigned short iTer = 0; iTer < input->Get_nIter(); iTer++){ 

        AD::SetDerivative(resp_sigmaboom, 1.0);

        structure->SetSolutionAdjoint(iRigid);

        AD::ComputeAdjoint();
        structure->ExtractSolutionAdjoint(iRigid);

        E_grad = input->GetGradient_E();
        Nu_grad = input->GetGradient_Nu();


        for (int iPDV = 0; iPDV < nPropDVs; iPDV++){
            propGradient[iPDV] = AD::GetValue(AD::GetDerivative(propDVsVector[iPDV]));
        }      

        for (iLoad = 0; iLoad < nTotalDOF; iLoad++){
            loadGradient[iLoad] = AD::GetValue(AD::GetDerivative(loadVector[iLoad]));
        }

        AD::ClearAdjoints();  
    }
}

void CBeamSolver::ComputeAdjointIzz_b(void){    
    unsigned long iLoad;    
    for (unsigned short iTer = 0; iTer < input->Get_nIter(); iTer++){ 
        AD::SetDerivative(resp_Izz_b, 1.0);
        structure->SetSolutionAdjoint(iRigid);
        AD::ComputeAdjoint();
        structure->ExtractSolutionAdjoint(iRigid);
        E_grad = input->GetGradient_E();
        Nu_grad = input->GetGradient_Nu();
        for (int iPDV = 0; iPDV < nPropDVs; iPDV++){
            propGradient[iPDV] = AD::GetValue(AD::GetDerivative(propDVsVector[iPDV]));
        }      
        for (iLoad = 0; iLoad < nTotalDOF; iLoad++){
            loadGradient[iLoad] = AD::GetValue(AD::GetDerivative(loadVector[iLoad]));
        }
        AD::ClearAdjoints();  
    }
}


//* Order d/d  U, E, Nu, Props, A,Loads  
void CBeamSolver::ComputeAdjointEA(void){

    unsigned long iLoad;

    AD::SetDerivative(resp_EA, 1.0);

    structure->SetSolutionAdjoint(iRigid);

    AD::ComputeAdjoint();
    structure->ExtractSolutionAdjoint(iRigid);
        
    E_grad = input->GetGradient_E();
    Nu_grad = input->GetGradient_Nu();
        
    
//    for (int iPDV = 0; iPDV < nPropDVs; iPDV++){
//        propGradient[iPDV] = AD::GetValue(AD::GetDerivative(propDVsVector[iPDV]));}       
        
    A_grad = Prop[0]->GetGradient_A();
    
    for (iLoad = 0; iLoad < nTotalDOF; iLoad++){
        loadGradient[iLoad] = AD::GetValue(AD::GetDerivative(loadVector[iLoad]));}
    
    AD::ClearAdjoints();
   
}




void CBeamSolver::StopRecording(void) {

    AD::RegisterOutput(objective_function);
    /** Register the solution as output **/
    structure->RegisterSolutionOutput(iRigid);
    AD::StopRecording();}

void CBeamSolver::StopRecordingNint(void) {

    AD::RegisterOutput(resp_Nint);
    /** Register the solution as output **/
    structure->RegisterSolutionOutput(iRigid);
    AD::StopRecording();}


void CBeamSolver::StopRecordingWeight(void) {

    AD::RegisterOutput(resp_weight);
    /** Register the solution as output **/
    structure->RegisterSolutionOutput(iRigid);
    AD::StopRecording();}


void CBeamSolver::StopRecordingKS(void) {

    AD::RegisterOutput(resp_KS);
    /** Register the solution as output **/
    structure->RegisterSolutionOutput(iRigid);
    AD::StopRecording();}

void CBeamSolver::StopRecordingSigmaBoom(void) {

    AD::RegisterOutput(resp_sigmaboom);
    /** Register the solution as output **/
    structure->RegisterSolutionOutput(iRigid);
    AD::StopRecording();}

void CBeamSolver::StopRecordingIzz_b(void) {
    AD::RegisterOutput(resp_Izz_b);
    /** Register the solution as output **/
    structure->RegisterSolutionOutput(iRigid);
    AD::StopRecording();}


void CBeamSolver::StopRecordingEA(void) {

    AD::RegisterOutput(resp_EA);
    /** Register the solution as output **/
    structure->RegisterSolutionOutput(iRigid);
    AD::StopRecording();}



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

