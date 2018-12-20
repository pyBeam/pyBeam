/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Ruben Sanchez, Rauno Cavallaro
 * 
 * Developers: Ruben Sanchez (SciComp, TU Kaiserslautern)
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
    
}

CBeamSolver::~CBeamSolver(void) {
	
}

void CBeamSolver::solve_beam(void){
	
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	std::cout << "A"   << input->Get_A() << std::endl;
	std::cout << "Iyy" << input->Get_Iyy() << std::endl;
	std::cout << "Izz" << input->Get_Izz() << std::endl;
	std::cout << "Jx"  << input->Get_Jx() << std::endl;

	ofstream echoes;
	echoes.open ("./output/echoes.out");         //   Some echo results than will be read directly from the post-processing routine

	// In the debug version these files are written //
#ifdef DEBG

	ofstream Xcoord;
	Xcoord.open ("./output/echo_Xcoord.out");         // X - coord

	ofstream udisp;
	udisp.open ("./output/echo_udisp.out");       // dU - disp

	ofstream myfile3;
	myfile3.open ("./output/echo_dX.out");       // dX

	ofstream myfile4;
	myfile4.open ("./output/echo_R_Re.out");       // R/Re

	ofstream echo_Fint;
	echo_Fint.open ("./output/echo_Fint.out");      //  Internal Forces

	//---------- Very experimental, to check how the "Delft's method" works
	ofstream echo_dFint;
	echo_dFint.open ("./output/echo_dFint.out");            //  Incremental Internal Forces method 1

	ofstream echo_dUel;             //  incremental elastic displacement in new coord
	echo_dUel.open ("./output/echo_dUel.out");

	ofstream echo_Fext;             //  External Forces in GLOBAL REF
	echo_Fext.open ("./output/echo_Fext.out");

	ofstream echo_Residual;             //  Residual in GLOBAL REF
	echo_Residual.open ("./output/echo_Residual.out");


#endif

	//======================     INITIAL ADJUSTEMENTS   ==================================

	//===============================================
	//
	//  Now the system is considered as a whole
	//
	//===============================================
	
	structure = new CStructure(input, element);

	structure->EchoCoord();

	std::cout << "Reading External Forces" << std::endl;
	structure->ReadForces(input->Get_Load());                                 // reads the external forces and assembly of Fnom

	#ifdef DEBG
		for (int id_fe=1; id_fe<=structure->nfem; id_fe++)
		{
			std::ofstream myfile4 ("./output/echo_R_Re.out", std::ios_base::out | std::ios_base::app);
			myfile4  << element[id_fe-1]->Rrig.block(0,0,3,3) << std::endl;
			myfile4  << element[id_fe-1]->R.block(0,0,3,3)  << std::endl;
		}

	#endif

	//===============================================
	//
	//     Writing Informations in the echo file
	//
	//===============================================

	echoes << input->Get_nNodes() << std::endl;   // Total number of Nodes
	echoes << input->Get_ConvCriteria() << std::endl;    // Convergence factor
	echoes << input->Get_Load() << std::endl;    // Load
	echoes << input->Get_Load() << std::endl;    // load steps
	echoes << input->Get_FollowerFlag() << std::endl;    // Follower Flag (0/1/2)


	std::chrono::steady_clock::time_point partial= std::chrono::steady_clock::now();

	//===============================================
	//
	// LOAD STEPPING
	//
	//===============================================
	
	std::cout << "#####    STARTING LOAD STEPPING   #####" << std::endl;

	double  lambda = 1.0;
	double dlambda =  1.0/input->Get_LoadSteps() ;
	unsigned long iIter;
	unsigned long totalIter = 0;
	int ls = 1;

	for  ( ls = 0; ls < input->Get_LoadSteps(); ls++)
	{
		std::cout << "==============================" << std::endl;
		std::cout << "==       LOAD STEP     " << ls << std::endl;

		lambda = dlambda*(ls+1);
		std::cout << "lambda = " << lambda << std::endl;

		//===============================================
		//               ITERATIVE SEQUENCE
		//===============================================
		
		bool converged = false;

		for (iIter = 0; iIter < input->Get_nIter(); iIter++)
		{
			std::cout << "   ----- ITERATION  " << iIter << "-----" << std::endl;


			/*--------------------------------------------------
			 *   Updates  Fext, Residual,
			*----------------------------------------------------*/

			// Update External Forces (if they are follower )
			structure->UpdateExtForces(lambda,input->Get_FollowerFlag());

			// Evaluate the Residual
			structure->EvalResidual();

			std::cout << " Residual = "  <<  structure->Residual.norm() << std::endl;

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

			double disp_factor =   structure->dU.norm()/input->Get_l();
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

	//  Writing the number of iterations in the echoes
	echoes <<  totalIter  << std::endl;

	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();

	std::cout << "Time difference INI-FIN  = " << chrono::duration_cast<chrono::microseconds>(end - begin).count() <<std::endl;
	std::cout << "Time difference INI-PART = " << chrono::duration_cast<chrono::microseconds>(partial - begin).count() <<std::endl;
	std::cout << "Time difference PART-FIN = " << chrono::duration_cast<chrono::microseconds>(end  - partial).count() <<std::endl;

	echoes << (chrono::duration_cast<chrono::microseconds>(partial - begin).count())/1.0e6 << std::endl;
	echoes << (chrono::duration_cast<chrono::microseconds>(end - partial).count())/1.0e6 << std::endl;

	std::cout << "LA COMMEDIA E' FINITA "  <<std::endl;
	
	
}

void CBeamSolver::Solve(void){

	//===============================================
	//  Initialize structural solver
	//===============================================
	
	structure = new CStructure(input, element);

	std::cout << "Reading External Forces" << std::endl;
	structure->ReadForces(input->Get_Load());

	//===============================================
	// LOAD STEPPING
	//===============================================
	
	std::cout << "#####    STARTING LOAD STEPPING   #####" << std::endl;

	double  lambda = 1.0;
	double dlambda =  1.0/input->Get_LoadSteps() ;
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

			double disp_factor =   structure->dU.norm()/input->Get_l();
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

