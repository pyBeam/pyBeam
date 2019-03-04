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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "../include/property.h"

using namespace std;
using namespace Eigen;

CProperty::CProperty(int ID) {PropertyID = ID; }

void CProperty::SetSectionProperties(passivedouble A_in, passivedouble Iyy_in, passivedouble Izz_in, passivedouble Jt_in) { 
    
    
    A = A_in;  Iyy = Iyy_in; Izz = Izz_in;  J0 = Iyy + Izz;  Jt = Jt_in;
    //cout << A << endl;
    //cout << Iyy << endl;
    //cout << Izz << endl;
    //cout << J0 << endl;
} 


CProperty::~CProperty(void) {};