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

#pragma once
#include <math.h>
#include <map>
#include <stdio.h>      /* printf, fopen */
/*We can think of using #include "./mpi_structure.hpp" and the command SU2_MPI::Error
 * once fully integrated with SU2 core
 */ 
#include <stdlib.h>     /* exit, EXIT_FAILURE */ 

#include "../include/types.h"





class CInput
{

private:
	
protected:

	//##################     Numerical Inputs     ###########################

     int nConstr;       // Total number of constraints
    MatrixXdDiff  Constr_matrix;  // COnstraint matrix [ NODE_ID DOF_ID ]
    
	unsigned long nNodes;	// Number of overall nodes along the wing (no collapsed)   // To be removed
	unsigned long nFEM;		// Number of finite elements     // To be removed

	unsigned short nDOF; 	// Number of degrees of freedom     // To be removed
	
	addouble load; 			// [N];              // To be removed
	int follower_flag;		// (0) Nonfollower (1) follower (2) approx follower
	unsigned long loadSteps;			// Number of load steps
	unsigned long nIter;				// Number of iterations

	addouble end_time;		// [sec] for SS calculation                    // To be removed
	addouble dt;   			// [sec] time increment for SS calculation     // To be removed

	//##############    Wing Inputs  ###########################
	// Units Sys: SI

	addouble t; 				// web & flange thickness [m]    // To be removed
	addouble h;				// web height [m]               // To be removed
	addouble b;				// flange width [m]             // To be removed
	addouble E; 				// Elastic modulus [GPa]  
	addouble Poiss; 			// Poisson Ratio
	addouble ro;				// Beam Density [kg/m^3]
	addouble G;				// Shear modulus
	addouble l; 				// Wing Length [m]                    // To be removed
	addouble A;				// cross section area               // To be removed
	addouble As_z; 			// z Effective shear area                 // To be removed
	addouble As_y;			// y Effective shear area                  // To be removed

	addouble Iyy, Izz;                                                          // To be removed
	addouble Jx; 				//Polar Moment of Inertia            // To be removed


	addouble Mwing;			//Wing's mass [kg]                             // To be removed
	addouble EIy, EIz, GJ, AE;                                                        // To be removed

	addouble Clalpha;  		//  recall pi = atan(1)*4;
	addouble Cldelta;

	//#################    Elements properties    ############################

	addouble le;      	//element length                                                // To be removed
	addouble m, m_e; 		//Element's mass                                  // To be removed
	addouble m_w, m_f; 	//web and flange mass                                     // To be removed
	addouble Ix, Iz; 		//[kg*m^2]                                            // To be removed

	//################     Convergence Parameters    ###########################

	addouble convCriteria;

	
public:

  CInput(int py_nPoint, int py_nElem);
  
  virtual ~CInput(void);
  
  void SetParameters();
  
  void SetWebThickness(passivedouble thickness) {t = thickness; }  //To Be Removed

  void SetWebHeight(passivedouble height) {h = height; }  //To Be Removed
  
  void SetFlangeWidth(passivedouble FlangeWidth) {b = FlangeWidth; }   //To Be Removed
 
  void SetYoungModulus(passivedouble YoungModulus) {E = YoungModulus; } 

  void SetPoisson(passivedouble Poisson) {Poiss = Poisson; }  
  
  void SetDensity(passivedouble Density) {ro = Density; }
    
  void SetBeamLength(passivedouble BeamLength) {l = BeamLength; }   //To Be Removed
  
  void SetLoad(passivedouble Load) {load = Load; }  //To Be Removed
  
  void SetFollowerFlag(int FollowerFlag) {follower_flag = FollowerFlag; cout << "Warning! Follower loads are not implemented yet!! FOLLOWER_FLAG =0" << endl; }  
  
  void SetLoadSteps(unsigned long LoadSteps) { loadSteps = LoadSteps; }   
  
  void SetNStructIter(unsigned long NStructIter) {nIter = NStructIter; }   
  
  void SetConvCriterium(passivedouble ConvCriterium) {convCriteria = ConvCriterium; }   

  void SetnConstr(int nconstr) {nConstr = nconstr; Constr_matrix = MatrixXdDiff::Zero(nconstr,2);};
  
  void SetSingleConstr(int iConstr, int node_id, int DOF_id) {Constr_matrix(iConstr,1 -1) = node_id; Constr_matrix(iConstr,2 -1) = DOF_id;};
  
  MatrixXdDiff  GetConstrMatrix() {return Constr_matrix;};
  
  passivedouble GetYoungModulus() {return AD::GetValue(E); }

  passivedouble GetPoisson() {return AD::GetValue(Poiss); }
  
  passivedouble GetShear() {return AD::GetValue(G); }
  
  passivedouble GetDensity() {return AD::GetValue(ro); }
  
  unsigned long Get_nNodes(void) { return nNodes; }  //To Be Removed
    
  unsigned long Get_nFEM(void) { return nFEM; }  //To Be Removed
  
  unsigned short Get_nDOF(void) { return nDOF; }  //To Be Removed
  
  unsigned short Get_FollowerFlag(void) { return follower_flag; }  
  
  unsigned long Get_LoadSteps(void) { return loadSteps;}
  
  unsigned long Get_nIter(void) { return nIter;}  
  
  addouble Get_Thickness(void) { return t;}  //To Be Removed
  
  addouble Get_Load(void) { return load;}  //To Be Removed
  
  addouble Get_l(void) { return l; }  //To Be Removed
  
  addouble Get_le(void) { return le; }  //To Be Removed
  
  addouble Get_Jx(void) { return Jx; }  //To Be Removed
  
  addouble Get_m_e(void) { return m_e; }  //To Be Removed
  
  addouble Get_A(void) { return A; }  //To Be Removed
      
  addouble Get_EIz(void) { return EIz; } //To Be Removed
        
  addouble Get_EIy(void) { return EIy; } //To Be Removed
  
  addouble Get_GJ(void) { return GJ; } //To Be Removed
  
  addouble Get_AE(void) { return AE; } //To Be Removed
  
  addouble Get_m(void) { return m; }  //To Be Removed
  
  addouble Get_Iyy(void) { return Iyy; }     //To Be Removed
          
  addouble Get_Izz(void) { return Izz; }  //To Be Removed
    

  addouble Get_ConvCriteria(void) { return convCriteria; }  

 
};


