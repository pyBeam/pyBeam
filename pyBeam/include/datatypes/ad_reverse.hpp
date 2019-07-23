/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rocco Bombardieri, Rauno Cavallaro 
 * 
 * File developers: Ruben Sanchez (SciComp, TU Kaiserslautern)
 *                  Tim Albring (SciComp, TU Kaiserslautern)
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

typedef codi::RealReverse addouble;

namespace AD{
	
  /*--- Reference to the tape ---*/
  extern addouble::TapeType& globalTape;
  typedef codi::ExternalFunctionHelper<addouble> ExtFuncHelper;
  extern ExtFuncHelper* FuncHelper;

  /*--- Register variables as input/output ---*/
  inline void RegisterInput(addouble &data) {AD::globalTape.registerInput(data);}
  inline void RegisterOutput(addouble& data) {AD::globalTape.registerOutput(data);}
  
  /*--- Activate/deactivate the tape ---*/
  inline void StartRecording() {AD::globalTape.setActive();}
  inline void StopRecording() {AD::globalTape.setPassive();}

  /*--- Clear the tape ---*/
  inline void ClearAdjoints() {AD::globalTape.clearAdjoints(); }
  
  /*--- Evaluate the tape ---*/
  inline void ComputeAdjoint() {AD::globalTape.evaluate(); }

  /*--- Set the value of a variable ---*/
  inline void SetValue(addouble& data, const double &val) {data.setValue(val);}
  inline double GetValue(const addouble& data) { return data.getValue();}
  
  /*--- Set the derivative of a variable ---*/
  inline void SetDerivative(addouble& data, const double &val) {data.setGradient(val);}    
  inline double GetDerivative(const addouble& data) { return data.getGradient();} 

  /*--- Reset the tape ---*/
  inline void Reset() { globalTape.reset(); }

  inline void SetExtFuncIn(addouble &data) {
    FuncHelper->addInput(data);
  }

  inline void SetExtFuncOut(addouble &data) {
    if (AD::globalTape.isActive()) {
      FuncHelper->addOutput(data);
    }
  }

  inline void StartExtFunc(bool storePrimalInput, bool storePrimalOutput){
    FuncHelper = new ExtFuncHelper(true);
    if (!storePrimalInput){
      FuncHelper->disableInputPrimalStore();
    }
    if (!storePrimalOutput){
      FuncHelper->disableOutputPrimalStore();
    }
  }
	
}
