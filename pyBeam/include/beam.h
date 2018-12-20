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

#pragma once
#include <math.h>       /* exp */

#include "../include/FiniteElement.h"
#include "../include/StructSyst.h"

class CBeamSolver
{

private:
	
protected:

	//##################     Numerical Inputs     ###########################

	int Nn;					// number of overall nodes along the wing (no collapsed)
	int nFEM;				// number of finite elements

	int DOF; 				// number of rigid modes to be calculated

	int lin;  				// flag for linear/nonlinear solution
	double load; 			// [N];
	int follower_flag;		// (0) Nonfollower (1) follower (2) approx follower
	int load_steps;			// Number of load steps
	int n_iter;				// Number of iterations

	double end_time;		// [sec] for SS calculation
	double dt;   			// [sec] time increment for SS calculation

	//#########################    Airfoil properties    ##########################

	double c;            	// Wing Chord
	double cg;          	// Airfoil Center of Garavity
	double ac;           	// Airfoil Aerodynamic Center
	double qc;
	double ea;           	// Elastic Axis Location

	double e;            	// Distance of E.A & A.C
	double d;				// Distance of C.G & A.C
	double a;

	//##############    Wing Inputs  ###########################
	// Units Sys: SI

	double t; 				// web & flange thickness [m]
	double h;				// web height [m]
	double b;				// flange width [m]
	double E; 				// Elastic modulus [GPa]
	double Poiss; 			// Poisson Ratio
	double ro;				// Beam Density [kg/m^3]
	double G;				// Shear modulus
	double l; 				// Wing Length [m]
	double A;				// cross section area
	double As_z; 			// z Effective shear area
	double As_y;			// y Effective shear area

	double Iyy, Izz;
	double Jx; 				//Polar Moment of Inertia


	double Mwing;			//Wing's mass [kg]
	double EIy, EIz, GJ, AE;

	double Clalpha;  		//  recall pi = atan(1)*4;
	double Cldelta;
	double S=l*c; 			//Wing Surface


	//#################    Elements properties    ############################

	double le;      	//element length
	double m, m_e; 		//Element's mass
	double m_w, m_f; 	//web and flange mass
	double Ix, Iz; 		//[kg*m^2]

	//################     Convergence Parameters    ###########################

	double conv_disp;

	
public:

  CElement** element;  /*!< \brief Vector which the define the elements. */

  CBeamSolver(void);
  
  virtual ~CBeamSolver(void);
  
  double GetX1(void) { return t; }  
    
  double GetX2(void) { return h; }
  
  double GetX3(void) { return b; }
  
  void solve_beam(void);     
  
};
