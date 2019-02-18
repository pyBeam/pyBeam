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
#include <Eigen/Dense>

#include "../include/types.h"

#include "../include/Rotations.h"
#include "../include/input.h"

class CElement
{
private:
public:
	int elemdofs = 12;

	addouble le; // beam length
	addouble Jx;
	addouble m_e; // Element mass (has to be avaluated here)
	addouble A;

	addouble EIz; // Has to be calc
	addouble EIy; // Has to be calc
	addouble GJ; // Has to be calc
	addouble AE; // Has to be calc

	addouble m;  // same as m_e to be removed
	addouble Iyy; // taken from property
	addouble Izz; // taken from property


	MatrixXdDiff  Rrig;         // Rodriguez Rotation Matrix (local reference system in GCS)
	MatrixXdDiff  R;             // Local reference system in GCS
	MatrixXdDiff  Rprev;         // m_R previous value

	// Added for Battini's CS
	addouble l_ini;                    // Initial Length
	addouble l_act;                    // Actual Length
	addouble l_prev;                   // Previous length


	MatrixXdDiff Mfem;
	VectorXdDiff fint;  // This is the elemen't Internal Forces

	VectorXdDiff eps ;  // Elastic Deformational Status
	VectorXdDiff phi ;  // Elastic Cumulative Tension
    MatrixXdDiff Kprim;






