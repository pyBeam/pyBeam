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
#include <math.h>       /* exp */

#include "../include/types.h"

#include "../include/element.h"
#include "../include/structure.h"
#include "../include/geometry.h"
#include "../include/input.h"

class CBeamSolver
{

private:

  addouble objective_function;
  bool register_loads;
  passivedouble *loadGradient;

  CNode **node;                     /*!< \brief Vector which stores the node initial coordinates. */
  //CConnectivity **connectivity;      /*!< \brief Vector which stores the connectivity. */

  CInput* input;

  CElement** element;  	  /*!< \brief Vector which the define the elements. */

  CStructure* structure;  /*!< \brief Pointer which the defines the structure. */

  int nDOF, nTotalDOF, nDim;
  unsigned long nFEM;
  addouble *loadVector;
  addouble thickness;

protected:

public:

  CBeamSolver(void);
  
  virtual ~CBeamSolver(void);
  
  void InitializeInput(CInput *py_input);

  void InitializeNode(CNode *py_node, unsigned long iNode);

  void InitializeElement(CElement *py_element, unsigned long iFEM);

  void InitializeStructure(void);

  void RegisterLoads(void);

  void Solve(void);

  passivedouble OF_NodeDisplacement(int iNode);

  passivedouble ComputeAdjoint(void);

  // Inlined functions

  //inline void SetThickness(passivedouble val_thickness) {thickness = val_thickness;}

  inline void SetLoads(int iNode, int iDOF, passivedouble loadValue) { loadVector[iNode*nDOF + iDOF] = loadValue; }

  inline passivedouble ExtractDisplacements(int iNode, int iDim) {return AD::GetValue(structure->GetDisplacement(iNode, iDim));}

  inline passivedouble ExtractCoordinates(int iNode, int iDim) {return AD::GetValue(structure->GetCoordinates(iNode, iDim));}

  inline passivedouble ExtractInitialCoordinates(int iNode, int iDim) {return AD::GetValue(structure->GetInitialCoordinates(iNode, iDim));}

  inline void StartRecording(void) { AD::StartRecording(); AD::RegisterInput(thickness);}

  inline void StopRecording(void) { AD::RegisterOutput(objective_function); AD::StopRecording(); }

  inline passivedouble ExtractLoadGradient(int iNode, int iDOF) {return loadGradient[iNode*nDOF + iDOF];}

  inline unsigned long Get_nNodes(void) {return input->Get_nNodes();}

};
