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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include "../include/FiniteElement.h"
#include "../include/Rotations.h"
#include "../include/input.h"

#include <iostream>

#ifdef DEBG
#include <fstream>
#endif

class CStructure
{
	
private:

public:

    int nfem;                // number of finite elements

	int DOF;                  // In space, 6
	int FollFlag;             // Flag for Follower forces (1)

	CElement **fem;      // Pointer to the first finite element

	Eigen::MatrixXd M;      // Recall in Eigen X stays for dynamic, d for double:  (nfem+1)*6  X   (nfem+1)*6
	Eigen::MatrixXd Ksys;


	Eigen::VectorXd dU;           // Displacement array (iterative)
	Eigen::VectorXd X;            // Position of the fem nodes in global coordinate system


	Eigen::VectorXd Fint;        // Array of internal forces
	Eigen::VectorXd Fext;        // Array of External Forces
	Eigen::VectorXd Residual;    // Array of Unbalanced Forces

	Eigen::Vector3d Ftip;      // (vector read from the input file) - needed as basis to update the Fext
	Eigen::VectorXd Fnom;        // Array of nominal forces

	CStructure(CInput *input, CElement **element);

	~CStructure();

	/*##############################################################
	 *
	 *         External Forces, Residual, Intenral Forces
	 *
	 *###############################################################*/

	void ReadForces(double forces);

	void UpdateExtForces(double , int );

	void EvalResidual();

	//===================================================
	//      Assembly System Stiffness Matrix
	//===================================================

	void AssemblyTang();

	void EvalSensRot();  // Evaluate the sensitivity of Rotation Matrix - need for Jacobian

	//===================================================
	//      Solve linear static system
	//===================================================
	// Assembles LHS and RHS and solves the linear static problem

	void SolveLinearStaticSystem();

	//===================================================
	//      Update Coordinates
	//===================================================
	/* This member function upates the coordinates (expressed in global reference system) of the finite element nodes.
	 * This is necessary for "booking" the position, as the compatiblity and the equations are based on the dispalcements.
	 */

	void UpdateCoord();

	void InitialCoord();

	void UpdateLength();

	void UpdateRotationMatrix();

	//===================================================
	//     TOOLS: WRITING COORDINATES
	//===================================================

	void EchoCoord();

	void EchoDisp();

	void EchoRes();

	void EchoFext();

	void EchoMatrixK();

	//===================================================
	//      INTERNAL FORCES
	//===================================================

	void UpdateInternalForces();




};
