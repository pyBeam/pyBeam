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
	
	addouble load; 			// [N];
	int follower_flag;		// (0) Nonfollower (1) follower (2) approx follower
	unsigned long loadSteps;			// Number of load steps
	unsigned long nIter;				// Number of iterations

	addouble end_time;		// [sec] for SS calculation
	addouble dt;   			// [sec] time increment for SS calculation

	//##############    Wing Inputs  ###########################
	// Units Sys: SI

	addouble t; 				// web & flange thickness [m]
	addouble h;				// web height [m]
	addouble b;				// flange width [m]
	addouble E; 				// Elastic modulus [GPa]
	addouble Poiss; 			// Poisson Ratio
	addouble ro;				// Beam Density [kg/m^3]
	addouble G;				// Shear modulus
	addouble l; 				// Wing Length [m]
	addouble A;				// cross section area
	addouble As_z; 			// z Effective shear area
	addouble As_y;			// y Effective shear area

	addouble Iyy, Izz;
	addouble Jx; 				//Polar Moment of Inertia


	addouble Mwing;			//Wing's mass [kg]
	addouble EIy, EIz, GJ, AE;

	addouble Clalpha;  		//  recall pi = atan(1)*4;
	addouble Cldelta;

	//#################    Elements properties    ############################

	addouble le;      	//element length
	addouble m, m_e; 		//Element's mass
	addouble m_w, m_f; 	//web and flange mass
	addouble Ix, Iz; 		//[kg*m^2]

	//################     Convergence Parameters    ###########################

	addouble convCriteria;

	
public:

  CInput(void);
  
  virtual ~CInput(void);
  
  void SetParameters(addouble thickness);
  
  unsigned long Get_nNodes(void) { return nNodes; }  
    
  unsigned long Get_nFEM(void) { return nFEM; }
  
  unsigned short Get_nDOF(void) { return nDOF; }
  
  unsigned short Get_FollowerFlag(void) { return follower_flag; }  
  
  unsigned long Get_LoadSteps(void) { return loadSteps;}
  
  unsigned long Get_nIter(void) { return nIter;}  
  
  addouble Get_Load(void) { return load;}  
  
  addouble Get_l(void) { return l; }  
  
  addouble Get_le(void) { return le; }
  
  addouble Get_Jx(void) { return Jx; }
  
  addouble Get_m_e(void) { return m_e; } 
  
  addouble Get_A(void) { return A; } 
      
  addouble Get_EIz(void) { return EIz; } 
        
  addouble Get_EIy(void) { return EIy; } 
  
  addouble Get_GJ(void) { return GJ; } 
  
  addouble Get_AE(void) { return AE; } 
  
  addouble Get_m(void) { return m; }
  
  addouble Get_Iyy(void) { return Iyy; }     
          
  addouble Get_Izz(void) { return Izz; }  
    
  addouble Get_ConvCriteria(void) { return convCriteria; }  
    
  
};

