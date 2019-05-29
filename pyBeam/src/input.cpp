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

#include <iostream>
#include <fstream>
#include <chrono>


#include "../include/input.h"

using namespace std;

CInput::CInput(int py_nPoint, int py_nElem) {

  nNodes = py_nPoint;
  nFEM = py_nElem;
  nRBE2 = 0;
  
}

CInput::CInput(int py_nPoint, int py_nElem, int py_nRBE2) {

  nNodes = py_nPoint;
  nFEM = py_nElem;
  nRBE2 = py_nRBE2;
  
}

void CInput::SetParameters(){	
    //##################     Numerical Inputs     ###########################
	
	nDOF = 6;                // To be removed

	//##############    Material inputs (only ONE homogeneous material is allowed by now)  ###########################
	// Units Sys: SI


	G = E/(2*(1+Poiss) );	// Shear modulus

}

CInput::~CInput(void) {
	
}


