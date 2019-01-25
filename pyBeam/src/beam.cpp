/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rauno Cavallaro
 * 
 * Developers: Tim Albring, Ruben Sanchez (SciComp, TU Kaiserslautern)
 *             Rauno Cavallaro (Carlos III University Madrid)
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
	
	input = new CInput();

}

CBeamSolver::~CBeamSolver(void) {
	
}

void CBeamSolver::Initialize(void){

    thickness = 1.99*1e-2;

    //su2double::TapeType& globalTape = su2double::getGlobalTape();

    //globalTape.setActive();

    //globalTape.registerInput(thickness);

    input->SetParameters(thickness);

    nDOF = input->Get_nDOF();
    nTotalDOF = input->Get_nNodes() * input->Get_nDOF();

    loadVector = new su2double[nTotalDOF];
    for (int iLoad = 0; iLoad < nTotalDOF; iLoad++)
        loadVector[iLoad] = 0.0;


    //==============================================================
    //      Finite Element initialization
    //==============================================================

    std::cout << "=========  Finite Element Initialization  ====" << std::endl;
    unsigned long nFEM = input->Get_nFEM();
    element = new CElement*[nFEM];
    for (unsigned long iFEM = 0; iFEM < nFEM; iFEM++){
        element[iFEM] = new CElement(iFEM, input);
    }

	structure = NULL;

	//===============================================
	//  Initialize structural solver
	//===============================================

	structure = new CStructure(input, element);


}

void CBeamSolver::SetLoads(int iNode, int iDOF, passivedouble loadValue){

    std::cout << "Load Value" << std::endl;
    std::cout << iNode << " " << iDOF << " " << loadValue << endl;

    int index;
    index = iNode*nDOF + iDOF;

    loadVector[index] = loadValue;

	std::cout << "Reading External Forces" << std::endl;
	structure->ReadForces(input->Get_Load());
	structure->ReadForces(nTotalDOF, loadVector);

}

void CBeamSolver::Solve(void){

	//===============================================
	// LOAD STEPPING
	//===============================================
	
	std::cout << "#####    STARTING LOAD STEPPING   #####" << std::endl;

	su2double  lambda = 1.0;
	su2double dlambda =  1.0/input->Get_LoadSteps() ;
	unsigned long iIter;
	unsigned long totalIter = 0;
	unsigned long loadStep = 1;

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
			structure->UpdateExtForces(lambda,input->Get_FollowerFlag());
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

			su2double disp_factor =   structure->dU.norm()/input->Get_l();
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
  
  //su2double pos1, pos2, pos3;
  //su2double grad_t;
  
  //pos1 = structure->GetDisplacement(100, 0);
  //pos2 = structure->GetDisplacement(100, 1);
  //pos3 = structure->GetDisplacement(100, 2);
  
  //cout << pos1 << " " << pos2 << " " << pos3 << endl;
  
  //globalTape.registerOutput(pos2);
  
  //globalTape.setPassive();
  
  //pos2.setGradient(1.0);
  
  //globalTape.evaluate();
  
  //grad_t = thickness.getGradient();
   
  //std::cout << " t' is " << grad_t << endl;

}

passivedouble CBeamSolver::ExtractDisplacements(int iNode, int iDim){

  passivedouble pos;

  pos = structure->GetDisplacement(iNode, iDim).getValue();

  return pos;

}

passivedouble CBeamSolver::ExtractCoordinates(int iNode, int iDim){

  passivedouble pos;

  pos = structure->GetCoordinates(iNode, iDim).getValue();

  return pos;

}

passivedouble CBeamSolver::ExtractInitialCoordinates(int iNode, int iDim){

  passivedouble pos;

  pos = structure->GetInitialCoordinates(iNode, iDim).getValue();

  return pos;

}


