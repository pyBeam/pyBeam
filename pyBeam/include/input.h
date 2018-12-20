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
#include <math.h>

class CInput
{

private:
	
protected:

	//##################     Numerical Inputs     ###########################

	unsigned long nNodes;	// Number of overall nodes along the wing (no collapsed)
	unsigned long nFEM;		// Number of finite elements

	unsigned short nDOF; 	// Number of degrees of freedom
	
	double load; 			// [N];
	int follower_flag;		// (0) Nonfollower (1) follower (2) approx follower
	int load_steps;			// Number of load steps
	int n_iter;				// Number of iterations

	double end_time;		// [sec] for SS calculation
	double dt;   			// [sec] time increment for SS calculation

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

	//#################    Elements properties    ############################

	double le;      	//element length
	double m, m_e; 		//Element's mass
	double m_w, m_f; 	//web and flange mass
	double Ix, Iz; 		//[kg*m^2]

	//################     Convergence Parameters    ###########################

	double conv_disp;

	
public:

  CInput(void);
  
  virtual ~CInput(void);
  
  unsigned long Get_nNodes(void) { return nNodes; }  
    
  unsigned long Get_nFEM(void) { return nFEM; }
  
  unsigned short Get_nDOF(void) { return nDOF; }
  
  unsigned short Get_le(void) { return le; }
  
  unsigned short Get_Jx(void) { return Jx; }
  
  unsigned short Get_m_e(void) { return m_e; } 
  
  unsigned short Get_A(void) { return A; } 
      
  unsigned short Get_EIz(void) { return EIz; } 
        
  unsigned short Get_EIy(void) { return EIy; } 
  
  unsigned short Get_GJ(void) { return GJ; } 
  
  unsigned short Get_AE(void) { return AE; } 
  
  unsigned short Get_m(void) { return m; }
  
  unsigned short Get_Iyy(void) { return Iyy; }     
          
  unsigned short Get_Izz(void) { return Izz; }  
    

  
};

