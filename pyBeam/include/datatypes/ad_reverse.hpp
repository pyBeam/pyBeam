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
typedef codi::ExternalFunctionHelper<addouble> ExtFuncHelper;

namespace AD{
	
  /*--- Reference to the tape ---*/

  extern ExtFuncHelper* FuncHelper;

  extern addouble::TapeType& globalTape;

  extern std::vector<addouble::GradientData> inputValues;
  extern int adjointVectorPosition;

  /*--- Register variables as input/output ---*/
  inline void RegisterInput(addouble &data) {AD::globalTape.registerInput(data);
                                             inputValues.push_back(data.getGradientData());}
  inline void RegisterOutput(addouble& data) {AD::globalTape.registerOutput(data);}
  
  /*--- Activate/deactivate the tape ---*/
  inline void StartRecording() {AD::globalTape.setActive();}
  inline void StopRecording() {AD::globalTape.setPassive();}

  /*--- Clear the tape ---*/
  inline void ClearAdjoints() {AD::globalTape.clearAdjoints(); }
  
  /*--- Evaluate the tape ---*/
  inline void ComputeAdjoint() {AD::globalTape.evaluate(); adjointVectorPosition = 0;}

  /*--- Set the value of a variable ---*/
  inline void SetValue(addouble& data, const double &val) {data.setValue(val);}
  inline double GetValue(const addouble& data) { return data.getValue();}
  
  /*--- Set the derivative of a variable ---*/
  inline void SetDerivative(addouble& data, const double &val) {data.setGradient(val);}    
  inline double GetDerivative(const addouble& data) { return AD::globalTape.getGradient(AD::inputValues[AD::adjointVectorPosition++]);}

  /*--- Reset the tape ---*/
  inline void Reset() {
    if (inputValues.size() != 0) {
      globalTape.reset();
      adjointVectorPosition = 0;
      inputValues.clear();
    }
    else{
      globalTape.reset();
    }
  }

  inline void SetExtFuncIn(const addouble &data) {AD::FuncHelper->addInput(data);}

  inline void SetExtFuncOut(addouble &data){
    if (AD::globalTape.isActive()) {
      AD::FuncHelper->addOutput(data);
    }
  }

  inline void StartExtFunc(bool storePrimalInput, bool storePrimalOutput){
    FuncHelper = new ExtFuncHelper(true);
    if (!storePrimalInput){
      AD::FuncHelper->disableInputPrimalStore();
    }
    if (!storePrimalOutput){
      AD::FuncHelper->disableOutputPrimalStore();
    }
  }

  inline void EndExtFunc() {delete AD::FuncHelper;}
	
}

class SolveAdjSys{

public:
  static void SolveSys(const codi::RealReverse::Real* x, codi::RealReverse::Real* x_b, size_t m,
                       const codi::RealReverse::Real* y, const codi::RealReverse::Real* y_b, size_t n,
                       codi::DataStore* d);
};
