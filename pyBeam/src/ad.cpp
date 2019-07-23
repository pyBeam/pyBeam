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

#ifdef CODI_REVERSE_TYPE
#include "../include/datatypes/ad_reverse.hpp"
#include "../include/types.h"
#elif CODI_FORWARD_TYPE
#include "../include/datatypes/ad_forward.hpp"
#else
#include "../include/datatypes/ad_passive.hpp"
#endif


#ifdef CODI_REVERSE_TYPE

namespace AD {

  /*--- Initialization of the external function helper ---*/
  ExtFuncHelper* FuncHelper;

  /*--- Initialization of the tape ---*/
  addouble::TapeType& globalTape = addouble::getGlobalTape();

}

void SolveAdjSys::SolveSys(const codi::RealReverse::Real* x, codi::RealReverse::Real* x_b, size_t m,
                           const codi::RealReverse::Real* y, const codi::RealReverse::Real* y_b, size_t n,
                           codi::DataStore* d) {

    int nNode_b;
    d->getData(nNode_b);

    unsigned short kind_linSol_b;
    d->getData(kind_linSol_b);

    /*--- Initialize vectors ---*/
    VectorXdDiff dU_bar;
    dU_bar = VectorXdDiff::Zero(nNode_b*6);

    VectorXdDiff Residual_bar;
    Residual_bar = VectorXdDiff::Zero(nNode_b*6);

    MatrixXdDiff Ksys_b;
    d->getData(Ksys_b);
    addouble a;

    /*--- Initialize the right-hand side with the gradient of the solution of the primal linear system ---*/

    for (unsigned long i = 0; i < n; i ++) {
      Residual_bar(i) = y_b[i];
    }

    switch(kind_linSol_b){
    case PartialPivLu:
      dU_bar = Ksys_b.transpose().partialPivLu().solve(Residual_bar); break;
    case FullPivLu:
      dU_bar = Ksys_b.transpose().fullPivLu().solve(Residual_bar); break;
    case HouseholderQr:
      dU_bar = Ksys_b.transpose().householderQr().solve(Residual_bar); break;
    case ColPivHouseholderQr:
      dU_bar = Ksys_b.transpose().colPivHouseholderQr().solve(Residual_bar); break;
    case FullPivHouseholderQr:
      dU_bar = Ksys_b.transpose().fullPivHouseholderQr().solve(Residual_bar); break;
    case LLT:
      dU_bar = Ksys_b.transpose().llt().solve(Residual_bar); break;
    case LDLT:
      dU_bar = Ksys_b.transpose().ldlt().solve(Residual_bar); break;
    default:
      dU_bar = Ksys_b.transpose().fullPivHouseholderQr().solve(Residual_bar); break;
    }

    for (unsigned long i = 0; i < n; i ++) {
      x_b[i] = AD::GetValue(dU_bar(i));
    }

}

#endif
