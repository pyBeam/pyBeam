/*
 * pyBeam, an open-source Beam Solver
 *
 * Copyright (C) 2019 by the authors
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

typedef double addouble;

namespace AD{

    /*--- Set the value of a variable ---*/
    inline void SetValue(addouble& data, const addouble &val) {data = val;}
    inline double GetValue(const addouble& data) { return data;}

    /*--- Set the derivative of a variable ---*/
    inline void SetDerivative(addouble &data, const addouble &val) {}
    inline double GetDerivative(const addouble& data) { return 0.0;}

    /*--- Overloaded functions ---*/
    inline void RegisterInput(addouble &data) { }
    inline void RegisterOutput(addouble& data) { }
    inline void StartRecording() { }
    inline void StopRecording() { }
    inline void ClearAdjoints() { }
    inline void ComputeAdjoint() { }
    inline void Reset() { }

    inline void SetExtFuncIn(addouble &data) { }
    inline void SetExtFuncOut(addouble &data) { }

    inline void StartExtFunc(bool storePrimalInput, bool storePrimalOutput) { }

    inline void EndExtFunc(){ }
}
