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
	
	//##################     Numerical Inputs     ###########################

	Nn = 11; 				// number of overall nodes along the wing (no collapsed)
	nFEM = Nn - 1;

	DOF = 6;                // number of rigid modes to be calculated
	
	load = 5000; 			// [N];
	follower_flag = 0;		// (0) Nonfollower (1) follower (2) approx follower
	load_steps = 1;			// Number of load steps
	n_iter = 30;			// Number of iterations  --30

	//##############    Wing Inputs  ###########################
	// Units Sys: SI

	t = 2*1e-2; 			// web & flange thickness [m]
	h = 40*1e-2;			// web height [m]
	b = 20*1e-2;			// flange width [m]
	E = 70*1e9; 			// Elastic modulus [GPa]
	Poiss = 0.3; 			// Poisson Ratio
	ro = 2.7e3;				// Beam Density [kg/m^3]
	G = E/(2*(1+Poiss) );	// Shear modulus
	l = 30; 				// Wing Length [m]
	A = t*h+2*t*b;			// cross section area
	As_z = t*h; 			// z Effective shear area
	As_y = 5/6.0*2*t*b;		// y Effective shear area

	Iyy = (pow(t,3)*h)/12.0 + 2.0* (  b*pow(t,3)/12.0 + b*t*pow( (h+t)/2.0 , 2)    );  //cross section Izz [m^4]
	Izz = ( h*pow(t,3) + 2*t*pow(b,3))/12.0;
	Jx = Iyy+Izz; 			//Polar Moment of Inertia


	Mwing = A*l*ro;		//Wing's mass [kg]
	EIy = E*Iyy;
	EIz = E*Izz;
	GJ = G*Jx;
	AE = A*E;

	Clalpha = 2*   atan(1)*4;  //  recall pi = atan(1)*4;
	Cldelta = 1;


	//#################    Elements properties    ############################

	le = l/(Nn-1);      	//element length
	m_e = Mwing/(Nn-1); 	//Element mass
	m = m_e; 				//Element's mass

	m_w = ro*le*t*h; 		//web mass
	m_f = ro*le*t*b; 		//flange mass
	
	Iz = ( m_w*(pow(t,2)+pow(le,2)) +
		 2*m_f*(pow(b,2)+pow(le,2)) )/12.0;     //[kg*m^2]
		 
	Ix = m_w/12.0*(pow(h,2)+pow(t,2))+
		 2*m_f*( (pow(t,2)+pow(b,2))/12+
		 pow((t/2+h/2),2)) ; 		//[kg*m^2]


	//################     Convergence Parameters    ###########################

	conv_disp = 1e-4;
	
	//==============================================================
	//
	//      Input initialization
	//
	//==============================================================	
	
	input = new CInput();

	//==============================================================
	//
	//      Finite Element initialization
	//
	//==============================================================

	std::cout << "=========  Finite Element Initialization  ====" << std::endl;		
	element = new CElement*[nFEM];
	for (unsigned long iFEM = 0; iFEM < nFEM; iFEM++){
		element[iFEM] = new CElement(iFEM, input);
	}
    
}

CBeamSolver::~CBeamSolver(void) {
	
}

void CBeamSolver::solve_beam(void){
	
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	std::cout << "A"   << A << std::endl;
	std::cout << "Iyy" << Iyy << std::endl;
	std::cout << "Izz" << Izz << std::endl;
	std::cout << "Jx"  << Jx << std::endl;


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


	unsigned long iFem, nFem;
	nFem = Nn - 1;

	int tempor;

	//==============================================================
	//
	//      FEM Matrices assembly
	//
	//==============================================================

	std::cout << "=========  Finite Element Initialization  ====" << std::endl;

	//===============================================
	//
	//  Now the system is considered as a whole
	//
	//===============================================

	StructSyst structsyst;

	structsyst.Initializer(nFem, DOF , follower_flag, element);    // Initialize structure

	structsyst.EchoCoord();

	std::cout << "Reading External Forces" << std::endl;
	structsyst.ReadForces(load);                                 // reads the external forces and assembly of Fnom

	#ifdef DEBG
		for (int id_fe=1; id_fe<=structsyst.nfem; id_fe++)
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

	echoes <<       Nn       << std::endl;   // Total number of Nodes
	echoes <<  conv_disp    << std::endl;    // Convergence factor
	echoes <<   load        << std::endl;    // Load
	echoes <<  load_steps   << std::endl;    // load steps
	echoes << follower_flag << std::endl;    // Follower Flag (0/1/2)


	std::chrono::steady_clock::time_point partial= std::chrono::steady_clock::now();


	//===============================================
	//
	// LOAD STEPPING
	//
	//===============================================
	
	std::cout << "#####    STARTING LOAD STEPPING   #####" << std::endl;

	double  lambda = 0;
	double dlambda =  1.0/load_steps ;
	int i = 0;
	int i_cum = 1;
	int ls = 1;

	for  ( ls=1; ls<=load_steps; ls++)
	{
		std::cout << "==============================" << ls << std::endl;
		std::cout << "==       LOAD STEP     " << ls << std::endl;

		lambda = dlambda*ls;
		std::cout << "lambda = " << lambda << std::endl;

		//===============================================
		//
		//               ITERATIVE SEQUENCE
		//
		//===============================================
		int converged = 0;
		i = 1;

		while ( i<=n_iter && converged == 0)
		{
			std::cout << "   ----- ITERATION  -----" << i << std::endl;


			/*--------------------------------------------------
			 *   Updates  Fext, Residual,
			*----------------------------------------------------*/

			// Update External Forces (if they are follower )
			structsyst.UpdateExtForces(lambda,follower_flag);
			structsyst.EchoFext();

			// Evaluate the Residual
			structsyst.EvalResidual();

			std::cout << " Residual = "  <<  structsyst.Residual.norm() << std::endl;
			structsyst.EchoRes();

			/*--------------------------------------------------
			 *   Assembly Ktang, Solve System
			*----------------------------------------------------*/

			// Reassembling Stiffness Matrix + Applying Boundary Conditions
			structsyst.AssemblyTang();

			// Solve Linear System   K*dU = Res = Fext - Fin
			structsyst.SolveLinearStaticSystem();

			/*--------------------------------------------------
			 *   Updates Coordinates, Updates Rotation Matrices
			*----------------------------------------------------*/

			structsyst.UpdateCoord();

			structsyst.EchoDisp();   structsyst.EchoCoord();

			// Now only X is updated!

			structsyst.UpdateRotationMatrix();  // based on the rotational displacements
			structsyst.UpdateLength();          // Updating length, important


			/*--------------------------------------------------
			 *   Update Internal Forces
			*----------------------------------------------------*/
            // Now, X, R, l are updated

			structsyst.UpdateInternalForces();

			/*--------------------------------------------------
			 *    Check Convergence
			*----------------------------------------------------*/

			double disp_factor =   structsyst.dU.norm()/l;
			std::cout << " disp_factor = "  <<  disp_factor << std::endl;

			if (disp_factor <= conv_disp)
			{
				converged = 1;
				std::cout << " CONVERGED at iteration  "  << i <<   std::endl;
			}
			else
				i+=1;
			    i_cum += 1;
		}
		std::cout << "#####    EXITING ITERATIVE SEQUENCE   #####" << std::endl;
	}

	//  Writing the number of iterations in the echoes
	echoes <<  i_cum  << std::endl;

	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();

	std::cout << "Time difference INI-FIN  = " << chrono::duration_cast<chrono::microseconds>(end - begin).count() <<std::endl;
	std::cout << "Time difference INI-PART = " << chrono::duration_cast<chrono::microseconds>(partial - begin).count() <<std::endl;
	std::cout << "Time difference PART-FIN = " << chrono::duration_cast<chrono::microseconds>(end  - partial).count() <<std::endl;

	echoes << (chrono::duration_cast<chrono::microseconds>(partial - begin).count())/1.0e6 << std::endl;
	echoes << (chrono::duration_cast<chrono::microseconds>(end - partial).count())/1.0e6 << std::endl;

	std::cout << "LA COMMEDIA E' FINITA "  <<std::endl;
	
	
}
