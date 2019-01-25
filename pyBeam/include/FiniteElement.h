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
#include <Eigen/Dense>

#include "../include/types.h"

#include "../include/Rotations.h"
#include "../include/input.h"

class CElement
{
private:
public:
	int elemdofs;

	su2double le;
	su2double Jx;
	su2double m_e;
	su2double A;

	su2double EIz;
	su2double EIy;
	su2double GJ;
	su2double AE;

	su2double m;
	su2double Iyy;
	su2double Izz;


	MatrixXdDiff  Rrig;         // Rodriguez Rotation Matrix (local reference system in GCS)
	MatrixXdDiff  R;             // Local reference system in GCS
	MatrixXdDiff  Rprev;         // m_R previous value

	// Added for Battini's CS
	su2double l_ini;                    // Initial Length
	su2double l_act;                    // Actual Length
	su2double l_prev;                   // Previous length


public:
	MatrixXdDiff Mfem;
	VectorXdDiff fint;  // This is the elemen't Internal Forces

	VectorXdDiff eps ;  // Elastic Deformational Status
	VectorXdDiff phi ;  // Elastic Cumulative Tension
    MatrixXdDiff Kprim;



private:


public:
	//	FiniteElement(){}
	
	CElement(void);	
	
	CElement(unsigned long iElement, CInput *input);
	
	virtual ~CElement(void);

	//	// Default constructors with also parameter definition
	//	Segment(int valore_init , su2double valore_kinit );
	// Default constructors with also parameter definition

	void Initializer(su2double val_le, su2double val_Jx , su2double val_m_e ,su2double val_A,
			su2double val_EIz , su2double val_EIy , su2double val_GJ, su2double val_AE ,
			su2double val_m ,  su2double val_Iyy, su2double val_Izz,
			int val_elemdofs=12)
	{
		le  = val_le;
		Jx  = val_Jx;
		m_e = val_m_e;
		A   = val_A;
		EIz = val_EIz;
		EIy = val_EIy;
		GJ  = val_GJ;
		AE  = val_AE;
		m   = val_m ;
		Iyy = val_Iyy;
		Izz = val_Izz;


		elemdofs = val_elemdofs;

		Mfem  = MatrixXdDiff::Zero(elemdofs,elemdofs);

		fint = VectorXdDiff::Zero(elemdofs);


		Rrig = MatrixXdDiff::Zero(6,6);         // Rotation Matrix
		R     = MatrixXdDiff::Identity(6,6);    // Initial Rotation Matrix
		Rprev = MatrixXdDiff::Identity(6,6);    // Initial Rotation Matrix

		l_act  = le;
		l_ini  = le;
		l_prev = le;

		eps  = VectorXdDiff::Zero(6);   // Elastic Cumulative deformation
        phi  = VectorXdDiff::Zero(6);   // Elastic Cumulative tension

        // INITIALIZATION of KPRIM

        VectorXdDiff diagonale = VectorXdDiff::Zero(6);
        Kprim = MatrixXdDiff::Zero(6,6);

        diagonale << AE/l_ini  , GJ/l_ini  ,  4*EIy/l_ini  ,   4*EIz/l_ini , 4*EIy/l_ini , 4*EIz/l_ini ;

        // Writing the diagonal
        for (int iii=1; iii<=6; iii++)
        {
            Kprim(iii-1,iii-1) = diagonale(iii-1);
        }

        Kprim(3-1,5-1) = 2*EIy/l_ini;  Kprim(5-1,3-1) = Kprim(3-1,5-1);
        Kprim(4-1,6-1) = 2*EIz/l_ini;  Kprim(6-1,4-1) = Kprim(4-1,6-1);

	}

	// Evaluates FEM element matrix
	void ElementMass_Rao();

	void EvalNaNb(MatrixXdDiff &Na,  MatrixXdDiff  &Nb);

	// Evaluates FEM element matrix
	void ElementElastic_Rao(MatrixXdDiff &Kel);

	// Evaluates FEM element matrix
	void ElementTang_Rao(MatrixXdDiff &Ktang);

	void EvalRotMat(VectorXdDiff &dU_AB , VectorXdDiff &X_AB );

	void EvalRotMatDEBUG(VectorXdDiff &dU_AB , VectorXdDiff &X_AB , MatrixXdDiff &Rtest);


};

