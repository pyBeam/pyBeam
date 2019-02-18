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

#include "../../externals/CoDiPack/include/codi.hpp"

typedef codi::RealForward addouble;

namespace AD{
  
  /*--- Set the value of a variable ---*/
  inline void SetValue(addouble& data, const double &val) {data.setValue(val);}
  inline double GetValue(const addouble& data) { return data.getValue();}
  
  /*--- Set the derivative of a variable ---*/
  inline void SetDerivative(addouble& data, const double &val) {data.setGradient(val);}
  inline double GetDerivative(const addouble& data) { return data.getGradient();}
  
  /*--- Overloaded functions ---*/
  inline void RegisterInput(addouble &data) { }
  inline void RegisterOutput(addouble& data) { }
  inline void StartRecording() { }
  inline void StopRecording() { }
  inline void ClearAdjoints() { }
  inline void ComputeAdjoint() { }
  inline void Reset() { }
}
