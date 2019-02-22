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


#include "../include/input.h"

using namespace std;

CInput::CInput(void) {
  
}

void CInput::SetParameters(){

	
    //##################     Numerical Inputs     ###########################

	nNodes = 21; 			// number of overall nodes along the wing (no collapsed)
	nFEM = nNodes - 1;
	
	nDOF = 6;                // To be removed
	//load = 5000; 			// [N];
	//follower_flag = 0;		// (0) Nonfollower (1) follower (2) approx follower
	//loadSteps = 1;			// Number of load steps
	//nIter = 30;			// Number of iterations  --30
	   
  
    nFEM = nNodes - 1;
	//##############    Wing Inputs  ###########################
	// Units Sys: SI

	//t = thickness;		// web & flange thickness [m]
	//h = 40*1e-2;			// web height [m]
	//b = 20*1e-2;			// flange width [m]
	//E = 70*1e9; 			// Elastic modulus [GPa]
	//Poiss = 0.3; 			// Poisson Ratio
	//ro = 2.7e3;				// Beam Density [kg/m^3]
	G = E/(2*(1+Poiss) );	// Shear modulus
	//l = 30; 				// Wing Length [m]
	A = t*h+2*t*b;			// cross section area
	//////As_z = t*h; 			// z Effective shear area
	///////As_y = 5/6.0*2*t*b;		// y Effective shear area

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

	le = l/(nNodes-1);      	//element length
	m_e = Mwing/(nNodes-1); 	//Element mass
	m = m_e; 				//Element's mass

	m_w = ro*le*t*h; 		//web mass
	m_f = ro*le*t*b; 		//flange mass
	
	/////Iz = ( m_w*(pow(t,2)+pow(le,2)) +
	/////	 2*m_f*(pow(b,2)+pow(le,2)) )/12.0;     //[kg*m^2]
		 
	/////Ix = m_w/12.0*(pow(h,2)+pow(t,2))+
	/////	 2*m_f*( (pow(t,2)+pow(b,2))/12+
	/////	 pow((t/2+h/2),2)) ; 		//[kg*m^2]

    
}

CInput::~CInput(void) {
	
}


