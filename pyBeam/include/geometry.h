/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rocco Bombardieri, Rauno Cavallaro 
 * 
 * File developers: Rocco Bombardieri (Carlos III University Madrid)
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

#pragma once
#include <iostream>
#include <fstream>
#include <chrono>

#include <math.h>

#include "../include/types.h"

class CNode 
{
private:

protected:

  unsigned long ID;  
    
  Vector3dDiff coord = VectorXdDiff::Zero(3);
  
  Vector3dDiff coordOld = VectorXdDiff::Zero(3);  
  
  Vector3dDiff coord0 = VectorXdDiff::Zero(3); 
  
  Vector3dDiff Vel = VectorXdDiff::Zero(3);

  Vector3dDiff Force = VectorXdDiff::Zero(3);

public:

  CNode( int id);

  ~CNode(void);

  inline void SetCoordinate(int iDim, passivedouble val_coor) {AD::SetValue(coord(iDim), val_coor);}

  inline void SetCoordinateOld(int iDim, passivedouble val_coor) {AD::SetValue(coordOld(iDim), val_coor);}
  
  inline void SetCoordinate0(int iDim, passivedouble val_coor) {AD::SetValue(coord0(iDim), val_coor);}
  
  inline void SetVel(int iDim, passivedouble val_vel) {AD::SetValue(Vel(iDim), val_vel);}
  
  inline void SetForce(int iDim, passivedouble val_force) {AD::SetValue(Force(iDim), val_force);}

#ifdef CODI_REVERSE_TYPE
  inline void SetCoordinate(int iDim, addouble val_coor) {coord(iDim) = val_coor;}

  inline void SetCoordinateOld(int iDim, addouble val_coor) {coordOld(iDim) = val_coor;}

  inline void SetCoordinate0(int iDim, addouble val_coor) {coord0(iDim) = val_coor;}
#endif
  
  inline addouble GetCoordinate(int iDim) {return coord(iDim);}

  inline addouble GetCoordinateOld(int iDim) {return coordOld(iDim);}  
  
  inline addouble GetCoordinate0(int iDim) {return coord0(iDim);}  
  
  inline addouble GetVel(int iDim) {return Vel(iDim);} 
  
  inline addouble GetForce(int iDim) {return Force(iDim);}   
  
  inline int GeID() {return ID;}    

};
