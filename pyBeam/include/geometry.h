/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rauno Cavallaro
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
  
  Vector3dDiff coord0 = VectorXdDiff::Zero(3); 
  
  Vector3dDiff Vel = VectorXdDiff::Zero(3);

  Vector3dDiff Force = VectorXdDiff::Zero(3);

public:

  CNode( int id);

  ~CNode(void);

  inline void SetCoordinate(int iDim, passivedouble val_coor) {coord(iDim) = val_coor;}

  inline void SetCoordinate0(int iDim, passivedouble val_coor) {coord0(iDim) = val_coor;}
  
  inline void SetVel(int iDim, passivedouble val_vel) {Vel(iDim) = val_vel;}  
  
  inline void SetForce(int iDim, passivedouble val_force) {Vel(iDim) = val_force;}   
  
  inline addouble GetCoordinate(int iDim) {return coord(iDim);}
  
  inline addouble GetCoordinate0(int iDim) {return coord0(iDim);}  
  
  inline addouble GetVel(int iDim) {return Vel(iDim);} 
  
  inline addouble GetForce(int iDim) {return Force(iDim);}   
  
  inline int GeID() {return ID;}    

};
/*

class CConnectivity {
private:

protected:

public:

  unsigned long nodeA;
  
  unsigned long nodeB;
  
  unsigned long property;
  
  VectorXdDiff aux_vector;

  CConnectivity(void);

  virtual ~CConnectivity(void);

  inline void SetNode_i(unsigned long node_i) {nodeA = node_i;};

  void SetNode_j(unsigned long node_j) {nodeB = node_j;};
  
  void SetProperty(unsigned long Prop) {property= Prop;} ; 

};
*/
