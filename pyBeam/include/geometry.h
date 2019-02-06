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
#include <math.h>

#include "../include/types.h"

class CNode {
private:

protected:

  addouble *coord;

public:

  CNode(void);

  virtual ~CNode(void);

  inline void SetCoordinate(int iDim, addouble val_coor) {coord[iDim] = val_coor;}

  inline addouble GetCoordinate(int iDim) {return coord[iDim];}

};


class CConnectivity {
private:

protected:

public:

  CNode **node;

  CConnectivity(void);

  virtual ~CConnectivity(void);

  inline void SetNode_i(CNode *node_i) {node[0] = node_i;}

  void SetNode_j(CNode *node_j) {node[1] = node_j;}

};
