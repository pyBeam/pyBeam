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
    
	unsigned long nNodes;	// Number of overall nodes along the wing (no collapsed)   
	unsigned long nFEM;		// Number of finite elements   
	unsigned long nRBE2;		// Number of RBE2 rigid elements          

	unsigned short nDOF; 	// Number of degrees of freedom     
        
        unsigned short rigid = 1; 	// [0] RBE2 [1] penalty
        addouble penalty = pow(10,11);   // penalty coefficient
	
	addouble load; 			// [N];              // To be removed
	int follower_flag;		// (0) Nonfollower (1) follower (2) approx follower
	unsigned long loadSteps;			// Number of load steps
	unsigned long nIter;				// Number of iterations

	addouble end_time;		// [sec] for SS calculation                    // To be removed
	addouble dt;   			// [sec] time increment for SS calculation     // To be removed

	//##############    Material inputs (only ONE homogeneous material is allowed by now)  ###########################
	// Units Sys: SI
	addouble E; 				// Elastic modulus [GPa]  
	addouble Poiss; 			// Poisson Ratio
	addouble ro;				// Beam Density [kg/m^3]
	addouble G;				// Shear modulus                                        // To be removed

	//################     Convergence Parameters    ###########################

	addouble convCriteria;

	
public:

    
  // Double constructor to make it compatible for old cases in which RBE2 were not defined  
  CInput(int py_nPoint, int py_nElem);  
    
  CInput(int py_nPoint, int py_nElem, int py_nRBE2);
  
  virtual ~CInput(void);
  
  void SetParameters();
  
  void SetYoungModulus(passivedouble YoungModulus) {E = YoungModulus; } 

  void SetPoisson(passivedouble Poisson) {Poiss = Poisson; }  
  
  void SetDensity(passivedouble Density) {ro = Density; }
    
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
  
  unsigned long Get_nNodes(void) { return nNodes; } 
    
  unsigned long Get_nFEM(void) { return nFEM; } 
  
  unsigned long Get_nRBE2(void) { return nRBE2; }   
  
  unsigned short Get_nDOF(void) { return nDOF; } 
  
  unsigned short Get_FollowerFlag(void) { return follower_flag; }  
  
  unsigned long Get_LoadSteps(void) { return loadSteps;}
  
  unsigned long Get_nIter(void) { return nIter;}  

  addouble Get_ConvCriteria(void) { return convCriteria; }  

  unsigned short Get_RigidCriteria(void) { return rigid; }  
  
  addouble GetPenalty(void) { return penalty; }   
  
 
};


