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
#include "../include/types.h"

#include "../include/rotations.h"
#include "../include/input.h"
#include "../include/property.h"
#include "../include/geometry.h"

class CElement
{
private:
public:
    
    int iElement;
    CNode* nodeA;
    CNode* nodeB;
    CProperty* property;
    CInput* input;
    
    int elemdofs = 12;   // Beam element DOFs
    VectorXdDiff GlobalDOFs  =  VectorXdDiff::Zero(12); // Global DOFs 
    
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
    Vector3dDiff aux_vector  =  Vector3dDiff::Zero(3);                // Versor directed from nodeA to nodeB
    
    
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
    CElement(int element_ID) ; // { iElement = element_ID; };
    
    ~CElement(void);
    
    //	// Default constructors with also parameter definition
    //	Segment(int valore_init , addouble valore_kinit );
    // Default constructors with also parameter definition
    
    inline void SetNode_1( CNode* Node1) { nodeA = Node1; };
    
    inline void SetNode_2( CNode* Node2) { nodeB = Node2;};    
    
    inline void SetProperty(CProperty* Property) {property = Property;};
    
    inline void SetInput(CInput* Input) {input = Input;};
    
    inline void SetAuxVector(addouble x, addouble y, addouble z) {aux_vector(0) = x; aux_vector(1) = y; aux_vector(2) = z;}
    
    void setGlobalDOFs();    
    
    void setLength();
    
    //addouble getLength() {return l_ini;}
    
    void setElementMass();        
    
    void Initializer(CNode* Node1, CNode* Node2, CProperty* Property, CInput* Input, passivedouble AuxVector_x, passivedouble AuxVector_y, passivedouble AuxVector_z);
    
    
    // Evaluates FEM element matrix
    void ElementMass_Rao();
    
    void EvalNaNb(MatrixXdDiff &Na,  MatrixXdDiff  &Nb);
    
    // Evaluates FEM element matrix
    void ElementElastic_Rao(MatrixXdDiff &Kel);
    
    // Evaluates FEM element matrix
    void ElementTang_Rao(int iIter, MatrixXdDiff &Ktang);
    
    void EvalRotMat(VectorXdDiff &dU_AB , VectorXdDiff &X_AB );
    
    void EvalRotMatDEBUG(VectorXdDiff &dU_AB , VectorXdDiff &X_AB , MatrixXdDiff &Rtest);
    
    void InitializeRotMats();
};

