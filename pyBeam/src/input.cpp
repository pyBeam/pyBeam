/*
 * pyBeam, an open-source Beam Solver
 *
 * Copyright (C) 2019 by the authors
 *
 * File developers: Rocco Bombardieri (Carlos III University Madrid)
 *                  Rauno Cavallaro (Carlos III University Madrid)
 *                  Ruben Sanchez (SciComp, TU Kaiserslautern)
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
    nFEM   = py_nElem;
    nRBE2  = py_nRBE2;

}

CInput::CInput(int py_nPoint, int py_nElem, int py_nRBE2, int py_nDV) {

    nNodes = py_nPoint;
    nFEM   = py_nElem;
    nRBE2  = py_nRBE2;
    nDV    = py_nDV;

}

CInput::CInput(int py_nPoint, int py_nElem, int py_nRBE2, int py_nDV,int py_nProp) {

    nNodes = py_nPoint;
    nFEM   = py_nElem;
    nRBE2  = py_nRBE2;
    nDV    = py_nDV;
    nProp  = py_nProp;
}



void CInput::SetParameters(){

    nDOF = 6;                // To be removed
    E = E_dimensional/E_dimensional;
    G = E/(2*(1+Poiss) );	// Shear modulus

}

CInput::~CInput(void) {

}
