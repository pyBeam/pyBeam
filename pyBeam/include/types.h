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
#include <cmath> 
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <iostream>
#include <vector>
#include <string>

#include "../externals/CoDiPack/include/codi.hpp"


using namespace std;


#ifdef CODI_REVERSE_TYPE
#include "./datatypes/ad_reverse.hpp"
#elif CODI_FORWARD_TYPE
#include "./datatypes/ad_forward.hpp"
#else
#include "./datatypes/ad_passive.hpp"
#endif

typedef double passivedouble;

// Redefine Eigen types depending on compilation
typedef Eigen::Matrix<addouble, Eigen::Dynamic, Eigen::Dynamic> MatrixXdDiff; // MatrixXd
typedef Eigen::Matrix<addouble, Eigen::Dynamic, 1> VectorXdDiff;              // VectorXd
typedef Eigen::Matrix<addouble, 3, 3> Matrix3dDiff;                           // Matrix3d
typedef Eigen::Matrix<addouble, 3, 1> Vector3dDiff;                           // Vector3d

