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

#pragma once
#include <math.h>

#include "../include/types.h"

class CInput
{

private:
	
protected:

	//##################     Numerical Inputs     ###########################

	unsigned long nNodes;	// Number of overall nodes along the wing (no collapsed)
	unsigned long nFEM;		// Number of finite elements

	unsigned short nDOF; 	// Number of degrees of freedom
	
	su2double load; 			// [N];
	int follower_flag;		// (0) Nonfollower (1) follower (2) approx follower
	unsigned long loadSteps;			// Number of load steps
	unsigned long nIter;				// Number of iterations

	su2double end_time;		// [sec] for SS calculation
	su2double dt;   			// [sec] time increment for SS calculation

	//##############    Wing Inputs  ###########################
	// Units Sys: SI

	su2double t; 				// web & flange thickness [m]
	su2double h;				// web height [m]
	su2double b;				// flange width [m]
	su2double E; 				// Elastic modulus [GPa]
	su2double Poiss; 			// Poisson Ratio
	su2double ro;				// Beam Density [kg/m^3]
	su2double G;				// Shear modulus
	su2double l; 				// Wing Length [m]
	su2double A;				// cross section area
	su2double As_z; 			// z Effective shear area
	su2double As_y;			// y Effective shear area

	su2double Iyy, Izz;
	su2double Jx; 				//Polar Moment of Inertia


	su2double Mwing;			//Wing's mass [kg]
	su2double EIy, EIz, GJ, AE;

	su2double Clalpha;  		//  recall pi = atan(1)*4;
	su2double Cldelta;

	//#################    Elements properties    ############################

	su2double le;      	//element length
	su2double m, m_e; 		//Element's mass
	su2double m_w, m_f; 	//web and flange mass
	su2double Ix, Iz; 		//[kg*m^2]

	//################     Convergence Parameters    ###########################

	su2double convCriteria;

	
public:

  CInput(void);
  
  virtual ~CInput(void);
  
  void SetParameters(su2double thickness);
  
  unsigned long Get_nNodes(void) { return nNodes; }  
    
  unsigned long Get_nFEM(void) { return nFEM; }
  
  unsigned short Get_nDOF(void) { return nDOF; }
  
  unsigned short Get_FollowerFlag(void) { return follower_flag; }  
  
  unsigned long Get_LoadSteps(void) { return loadSteps;}
  
  unsigned long Get_nIter(void) { return nIter;}  
  
  su2double Get_Load(void) { return load;}  
  
  su2double Get_l(void) { return l; }  
  
  su2double Get_le(void) { return le; }
  
  su2double Get_Jx(void) { return Jx; }
  
  su2double Get_m_e(void) { return m_e; } 
  
  su2double Get_A(void) { return A; } 
      
  su2double Get_EIz(void) { return EIz; } 
        
  su2double Get_EIy(void) { return EIy; } 
  
  su2double Get_GJ(void) { return GJ; } 
  
  su2double Get_AE(void) { return AE; } 
  
  su2double Get_m(void) { return m; }
  
  su2double Get_Iyy(void) { return Iyy; }     
          
  su2double Get_Izz(void) { return Izz; }  
    
  su2double Get_ConvCriteria(void) { return convCriteria; }  
    
  
};

