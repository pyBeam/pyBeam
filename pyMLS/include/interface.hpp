/*
 * pyMLS, an open-source Moving Least Squares library
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


#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iostream> 
#include <vector>
#include <fstream>

// External open-source dependencies
#include "../externals/ann/include/ANN/ANN.h"
#include "../externals/libigl/include/igl/slice.h"

void mls_interface (std::vector<double> &interpolation_matrix_std,
                    std::vector<double> &norm_err_std,
                    int str_nodenumb,
                    int aero_nodenumb,
                    std::vector<double> str_data_std,
                    std::vector<double> aero_data_std,
                    int poly, int weight, long int points,
                    double rmax,double delta,  double toll);
