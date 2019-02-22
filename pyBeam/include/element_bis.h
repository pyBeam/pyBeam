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
#include <Eigen/Dense>

#include "../include/types.h"

#include "../include/rotations.h"
#include "../include/input.h"
#include "../include/property.h"
#include "../include/geometry.h"

class CElement
{
private:
public:
    
    
    int nodeIndexA;
    int nodeIndexB;
    CNode* nodeA;
    CNode* nodeB;
    CProperty* property;
    CInput* input;
    
    int elemdofs = 12;   // Beam element DOFs
    
    addouble le;
    addouble J0;
    addouble m_e;
    addouble A;
    
    addouble EIz;
    addouble EIy;
    addouble GJ;
    addouble AE;
    
    addouble m;
    addouble Iyy;
    addouble Izz;
    
    
    MatrixXdDiff  Rrig;         // Rodriguez Rotation Matrix (local reference system in GCS)
    MatrixXdDiff  R;             // Local reference system in GCS
    MatrixXdDiff  Rprev;         // m_R previous value
    
    // Added for Battini's CS
    addouble l_ini;                    // Initial Length
    addouble l_act;                    // Actual Length
    addouble l_prev;                   // Previous length
    VectorXdDiff aux_vector;                  // Versor directed from nodeA to nodeB
    VectorXdDiff  =  VectorXdDiff::Zero(6);
    
public:
    MatrixXdDiff Mfem;
    VectorXdDiff fint;  // This is the elemen't Internal Forces
    
    VectorXdDiff eps ;  // Elastic Deformational Status
    VectorXdDiff phi ;  // Elastic Cumulative Tension
    MatrixXdDiff Kprim;
    
    
    
private:
    
    
public:
    //	FiniteElement(){}
    
    //CElement(void);	
    
    // In the constructor we assign to the element its nodes and properties
    CElement(unsigned long iElement, CInput *Input, CNode* Node1, CNode* Node2, CProperty* Property  ) {
        nodeA = Node1; nodeB = Node2; property = Property; input = Input;
    };
    
    virtual ~CElement(void);
    
    //	// Default constructors with also parameter definition
    //	Segment(int valore_init , addouble valore_kinit );
    // Default constructors with also parameter definition
    
    void setGlobalDOFs();    
    
    void setLength();
    
    void setElementMass( ro);        
    
    void Initializer();
    
    // Evaluates FEM element matrix
    void ElementMass_Rao();
    
    void EvalNaNb(MatrixXdDiff &Na,  MatrixXdDiff  &Nb);
    
    // Evaluates FEM element matrix
    void ElementElastic_Rao(MatrixXdDiff &Kel);
    
    // Evaluates FEM element matrix
    void ElementTang_Rao(MatrixXdDiff &Ktang);
    
    void EvalRotMat(VectorXdDiff &dU_AB , VectorXdDiff &X_AB );
    
    void EvalRotMatDEBUG(VectorXdDiff &dU_AB , VectorXdDiff &X_AB , MatrixXdDiff &Rtest);
    
    
};

