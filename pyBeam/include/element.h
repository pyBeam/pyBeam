/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rocco Bombardieri, Rauno Cavallaro 
 * 
 * File developers: Rocco Bombardieri (Carlos III University Madrid)
 *                  Rauno Cavallaro (Carlos III University Madrid)
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

    int iElement;

    // Added for Battini's CS
    addouble l_ini;              // Initial Length
    addouble l_curr;             // Current Length
    addouble l_prev;             // Previous length

    addouble m_e;                // Element mass
    addouble A;                  // Element area

    addouble Iyy;                // Beam inertia over y axis
    addouble Izz;                // Beam inertia over y axis

    addouble GJ;                 // Torsional stiffness
    addouble J0;                 // Polar inertia

    addouble AE;
    addouble EIz;
    addouble EIy;

public:

    CNode* nodeA;
    CNode* nodeB;
    CProperty* elprop;
    CInput* input;

    VectorXdDiff GlobalDOFs; // Global DOFs

    MatrixXdDiff Rrig;          // Rodriguez Rotation Matrix (local reference system in GCS)
    MatrixXdDiff R;             // Local reference system in GCS
    MatrixXdDiff Rprev;         // m_R previous value

    Vector3dDiff aux_vector;   // Versor directed from nodeA to nodeB

    MatrixXdDiff Mfem;
    VectorXdDiff fint;  // This is the elemen't Internal Forces
    
    VectorXdDiff eps ;  // Elastic Deformational Status
    VectorXdDiff phi ;  // Elastic Cumulative Tension
    MatrixXdDiff Kprim;
    
private:
    
    
public:
    
    // In the constructor we assign to the element its nodes and properties
    CElement(int element_ID) ; 
    
    ~CElement(void);

    // Methods to access initial element length
    inline addouble GetInitial_Length(void) { return l_ini; }

    // Methods to access current element length
    inline void SetCurrent_Length( addouble val_length) { l_curr = val_length; }
    inline addouble GetCurrent_Length(void) { return l_curr; }

    // Methods to access/set the previous element length
    inline void SetPrevious_Length(void) { l_prev = l_curr; }
    inline addouble GetPrevious_Length(void) { return l_prev; }

    // Methods to set initial properties to the element
    inline void SetNode_1( CNode* Node1) { nodeA = Node1; }
    inline void SetNode_2( CNode* Node2) { nodeB = Node2;}
    inline void SetProperty(CProperty* Property) {elprop = Property;}
    inline void SetInput(CInput* Input) {input = Input;}
    
    inline void SetAuxVector(addouble x, addouble y, addouble z) {aux_vector(0) = x; aux_vector(1) = y; aux_vector(2) = z;}

    inline void setElementMass() {m_e = elprop->GetA()*l_ini* input->GetDensity();}
    
    void setGlobalDOFs();    
    
    void setLength();
    
    void Initializer(CNode* Node1, CNode* Node2, CProperty* Property, CInput* Input, passivedouble AuxVector_x, passivedouble AuxVector_y, passivedouble AuxVector_z);
    
    // Evaluates FEM element mass matrix
    void ElementMass_Rao();

    // Evaluates FEM kinematic matrix
    void EvalNaNb(MatrixXdDiff &Na,  MatrixXdDiff  &Nb);
    
    // Evaluates FEM element elastic matrix
    void ElementElastic_Rao(MatrixXdDiff &Kel);
    
    // Evaluates FEM element stiffness matrix
    void ElementTang_Rao(int iIter, MatrixXdDiff &Ktang);

    // Evaluates FEM element rotation matrix
    void EvalRotMat(VectorXdDiff &dU_AB , VectorXdDiff &X_AB );

    // Initially rotates the elements
    void InitializeRotMats();    

    // Evaluates the rotation matrix for a value of a small dU_AB_eps around a given position of the reference system evaluated with current positions X_AB and rotations dU_AB X_AB
    //void EvalRotMatFiniteDifferences(VectorXdDiff dU_AB_eps, VectorXdDiff  X_AB, Matrix3dDiff &R_eps);
    
    // Set the element dependencies
    void SetDependencies(void);

};

