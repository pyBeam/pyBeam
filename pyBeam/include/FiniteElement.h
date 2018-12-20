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

#include "../include/Rotations.h"

class CElement
{
private:
public:
	int elemdofs;

	double le;
	double Jx;
	double m_e;
	double A;

	double EIz;
	double EIy;
	double GJ;
	double AE;

	double m;
	double Iyy;
	double Izz;


	Eigen::MatrixXd  Rrig;         // Rodriguez Rotation Matrix (local reference system in GCS)
	Eigen::MatrixXd  R;             // Local reference system in GCS
	Eigen::MatrixXd  Rprev;         // m_R previous value

	// Added for Battini's CS
	double l_ini;                    // Initial Length
	double l_act;                    // Actual Length
	double l_prev;                   // Previous length


public:
	Eigen::MatrixXd Mfem;
	Eigen::VectorXd fint;  // This is the elemen't Internal Forces

	Eigen::VectorXd eps ;  // Elastic Deformational Status
	Eigen::VectorXd phi ;  // Elastic Cumulative Tension
    Eigen::MatrixXd Kprim;



private:


public:
	//	FiniteElement(){}
	
	CElement(void);
	
	virtual ~CElement(void);

	//	// Default constructors with also parameter definition
	//	Segment(int valore_init , double valore_kinit );
	// Default constructors with also parameter definition

	void Initializer(double val_le, double val_Jx , double val_m_e ,double val_A,
			double val_EIz , double val_EIy , double val_GJ, double val_AE ,
			double val_m ,  double val_Iyy, double val_Izz,
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

		Mfem  = Eigen::MatrixXd::Zero(elemdofs,elemdofs);

		fint = Eigen::VectorXd::Zero(elemdofs);


		Rrig = Eigen::MatrixXd::Zero(6,6);         // Rotation Matrix
		R     = Eigen::MatrixXd::Identity(6,6);    // Initial Rotation Matrix
		Rprev = Eigen::MatrixXd::Identity(6,6);    // Initial Rotation Matrix

		l_act  = le;
		l_ini  = le;
		l_prev = le;

		eps  = Eigen::VectorXd::Zero(6);   // Elastic Cumulative deformation
        phi  = Eigen::VectorXd::Zero(6);   // Elastic Cumulative tension

        // INITIALIZATION of KPRIM

        Eigen::VectorXd diagonale = Eigen::VectorXd::Zero(6);
        Kprim = Eigen::MatrixXd::Zero(6,6);

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

	void EvalNaNb(Eigen::MatrixXd &Na,  Eigen::MatrixXd  &Nb);

	// Evaluates FEM element matrix
	void ElementElastic_Rao(Eigen::MatrixXd &Kel);

	// Evaluates FEM element matrix
	void ElementTang_Rao(Eigen::MatrixXd &Ktang);

	void EvalRotMat(Eigen::VectorXd &dU_AB , Eigen::VectorXd &X_AB );

	void EvalRotMatDEBUG(Eigen::VectorXd &dU_AB , Eigen::VectorXd &X_AB , Eigen::MatrixXd &Rtest);


};

