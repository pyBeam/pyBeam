/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rocco Bombardieri, Rauno Cavallaro 
 * 
 * File developers: Rauno Cavallaro (Carlos III University Madrid)
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

#include "../include/types.h"

#include <Eigen/Dense>

#include <iostream>

void RotToPseudo(Vector3dDiff& pseudo , Matrix3dDiff R);

void PseudoToRot(Vector3dDiff pseudo , Matrix3dDiff& R,  int print=0);

//void PseudoToRotDer(Vector3dDiff pseudo , Matrix3dDiff& dR_1, Matrix3dDiff& dR_2, Matrix3dDiff& dR_3)


