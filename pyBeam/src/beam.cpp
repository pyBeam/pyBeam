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
    
    //input = new CInput();
    
    register_loads = false;
    objective_function = 0.0;
    
}

CBeamSolver::~CBeamSolver(void) {
    
}

void CBeamSolver::InitializeInput(CInput* py_input){   // insert node class and connectivity

  // I'm memorizing as a member variable the object input passed from outside
  input = py_input;

  nDOF = input->Get_nDOF();
  nFEM = input->Get_nFEM();
  nTotalDOF = input->Get_nNodes() * input->Get_nDOF();

  //==============================================================
  //      Load Vector initialization
  //==============================================================

  loadVector = new addouble[nTotalDOF];
  for (int iLoad = 0; iLoad < nTotalDOF; iLoad++)
    loadVector[iLoad] = 0.0;

  //==============================================================
  //      Node vector initialization
  //==============================================================

  cout << "=========  Node Structure Initialization  ====" << std::endl;
  unsigned long nNodes = input->Get_nNodes();   //substitute from mesh file
  node = new CNode*[nNodes];

  //==============================================================
  //      Finite Element vector initialization
  //==============================================================

  cout << "=========  Finite Element Initialization  ====" << std::endl;
  unsigned long nFEM = input->Get_nFEM();   //substitute from mesh file
  element = new CElement*[nFEM];

  //===============================================
  //  Initialize structural solver
  //===============================================
  structure = NULL;

}


void CBeamSolver::InitializeNode(CNode *py_node, unsigned long iNode) {node[iNode] = py_node;}

void CBeamSolver::InitializeElement(CElement* py_element,unsigned long iFEM) {element[iFEM] = py_element;}

void CBeamSolver::InitializeStructure(void){structure = new CStructure(input, element, node);}

void CBeamSolver::Solve(void){

  // Beam total length
  addouble TotalLength = 0;
  for  ( unsigned long iFEM = 0; iFEM < nFEM; iFEM++)   
  {
      TotalLength += element[iFEM]->l_ini;
  }
    
    
  std::cout << "#####    SETTING EXTERNAL FORCES   #####" << std::endl;
  structure->ReadForces(nTotalDOF, loadVector); // to be interfaced with the aerodynamic part

	//===============================================
	// LOAD STEPPING
	//===============================================
	
	std::cout << "#####    STARTING LOAD STEPPING   #####" << std::endl;

  addouble  lambda = 1.0;
  addouble dlambda =  1.0/input->Get_LoadSteps() ;
	unsigned long iIter;
	unsigned long totalIter = 0;
	unsigned long loadStep = 1;

	structure->InitialCoord();

	for  ( loadStep = 0; loadStep < input->Get_LoadSteps(); loadStep++)
	{
		std::cout << "==============================" << std::endl;
		std::cout << "==       LOAD STEP     " << loadStep << std::endl;
		std::cout << "==============================" << std::endl;

		lambda = dlambda*(loadStep+1);
		std::cout << "lambda = " << lambda << std::endl;

		//===============================================
		//               ITERATIVE SEQUENCE
		//===============================================
		bool converged = false;

		for (iIter = 0; iIter < input->Get_nIter(); iIter++) {
			
			std::cout << "   ----- ITERATION  -----" << iIter << std::endl;

			/*--------------------------------------------------
			 *   Updates  Fext, Residual,
			 *----------------------------------------------------*/

			// Update External Forces (if they are follower )
			structure->UpdateExtForces(lambda);
			structure->EchoFext();

			// Evaluate the Residual
			structure->EvalResidual();

			std::cout << " Residual = "  <<  structure->Residual.norm() << std::endl;
			structure->EchoRes();

			/*--------------------------------------------------
			 *   Assembly Ktang, Solve System
			*----------------------------------------------------*/

			// Reassembling Stiffness Matrix + Applying Boundary Conditions
			structure->AssemblyTang();

			// Solve Linear System   K*dU = Res = Fext - Fin
			structure->SolveLinearStaticSystem();

			/*--------------------------------------------------
			 *   Updates Coordinates, Updates Rotation Matrices
			*----------------------------------------------------*/

			structure->UpdateCoord();

			structure->EchoDisp();   structure->EchoCoord();

			// Now only X is updated!

			structure->UpdateRotationMatrix();  // based on the rotational displacements
			structure->UpdateLength();          // Updating length, important


			/*--------------------------------------------------
			 *   Update Internal Forces
			*----------------------------------------------------*/
            // Now, X, R, l are updated
            
            structure->UpdateInternalForces();
            
            /*--------------------------------------------------
             *    Check Convergence
             *----------------------------------------------------*/

            addouble disp_factor =   structure->dU.norm()/TotalLength;
            
            std::cout << " disp_factor = "  <<  disp_factor << std::endl;
            
            if (disp_factor <= input->Get_ConvCriteria())
            {
                converged = true;
                std::cout << " CONVERGED at iteration  "  << iIter <<   std::endl;
                totalIter += iIter;
                break;
            }
        }
        std::cout << "#####    EXITING ITERATIVE SEQUENCE   #####" << std::endl;
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

void CBeamSolver::RegisterLoads(void){
    
    register_loads = true;
    
    /*--- Initialize vector to store the gradient ---*/
    loadGradient = new passivedouble[nTotalDOF];
    
    for (int iLoad = 0; iLoad < nTotalDOF; iLoad++){
        
        /*--- Register the load vector ---*/
        AD::RegisterInput(loadVector[iLoad]);
        
        /*--- Initialize the load gradient to 0 ---*/
        loadGradient[iLoad] = 0.0;
    }
    
}

void CBeamSolver::ComputeAdjoint(void){
       
    AD::SetDerivative(objective_function, 1.0);
    
    AD::ComputeAdjoint();
    
    if (register_loads){
        for (int iLoad = 0; iLoad < nTotalDOF; iLoad++){
            loadGradient[iLoad] = AD::GetValue(AD::GetDerivative(loadVector[iLoad]));
        }
    }
}



