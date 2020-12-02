/*
 * pyBeam, an open-source Beam Solver
 *
 * Copyright (C) 2019 by the authors
 *
 * File developers: Rocco Bombardieri (Carlos III University Madrid)
 *                  Rauno Cavallaro (Carlos III University Madrid)
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
#include <math.h>       /* exp */

#include "../include/types.h"

#include "../include/element.h"
#include "../include/rigid_element.h"
#include "../include/structure.h"
#include "../include/geometry.h"
#include "../include/input.h"

class CBeamSolver
{

private:

    addouble objective_function;
    addouble resp_weight;
    addouble resp_KS_stress ;
    addouble  resp_KS_buckl ;    //DBG
    addouble resp_EA;
    addouble resp_Izz_b;    
    addouble resp_Nint;
    addouble resp_sigmaboom;
    //END BDG 
    bool register_loads;
    
    passivedouble *loadGradient;
    passivedouble E_grad, Nu_grad;
    passivedouble *propGradient;
    
    ///DEBUG
    passivedouble A_grad;  
    

    bool verbose = true;

    CNode **node;                         /*!< \brief Vector which stores the node initial coordinates. */
    //CConnectivity **connectivity;       /*!< \brief Vector which stores the connectivity. */

    CInput* input;

    CElement** element;                   /*!< \brief Vector which the define the elements. */

    CRBE2** RBE2;                         /*!< \brief Vector which the define the elements. */

    CStructure* structure;                /*!< \brief Pointer which the defines the structure. */
    
    CProperty** Prop;                     /*!< \brief Pointer to the Properties. */

    
    int nDOF, nTotalDOF, nRBE2, nDim;
    int iRigid = 0;     // iRigid = 0: Penalty method iRigid = 1: Lagrange Multiplier method
    unsigned long nFEM;
    unsigned long nProp;                    ///< Number of Properties
    int nPropDVs; 
    
    unsigned long totalIter;
    addouble initDispNorm;
    addouble initResNorm;
    addouble *loadVector;
    addouble thickness;
    addouble *propDVsVector;   ///< property DVs gathered in a vector
    
//protected:

public:

    CBeamSolver(void);

    virtual ~CBeamSolver(void);

    void InitializeInput(CInput *py_input);

    inline void InitializeNode(CNode *py_node, unsigned long iNode) {node[iNode] = py_node;}

    inline void InitializeElement(CElement *py_element, unsigned long iFEM) {element[iFEM] = py_element;}

    inline void InitializeRBE2(CRBE2* py_RBE2,unsigned long iRBE2) {RBE2[iRBE2] = py_RBE2;}
    
    inline void InitializeProp(CProperty* py_prop,unsigned long iProp) {Prop[iProp] = py_prop;}

    void InitializeStructure(void);

    void Solve(int FSIIter);
    
    void SolveLin(int FSIIter);    

    passivedouble OF_NodeDisplacement(int iNode);
    
    passivedouble EvalWeight();
       
    passivedouble EvalKSStress();  
    
    passivedouble EvalKSBuckling();
    
    passivedouble EvalEA();   ///DEBUG 
    
    passivedouble EvalIzz_b();   ///DEBUG     
    
    passivedouble EvalNint();   /// FOR DEBUG PURPOSES
    
    passivedouble EvalSigmaBoom();   ///DEBUG    

    void ComputeAdjoint(void);
    
    void ComputeAdjointNint(void);   /// DEBUG    
    void ComputeAdjointSigmaBoom(void);   /// DEBUG  
    void ComputeAdjointEA(void);       ///DEBUG 
    void ComputeAdjointIzz_b(void);       ///DEBUG 
    
    void ComputeAdjointWeight(void);  
    
    void ComputeAdjointKSstresses(void); 
    void ComputeAdjointKSbuckling(void); 
    
    
   
    

    inline void SetLoads(int iNode, int iDOF, passivedouble loadValue) { loadVector[iNode*nDOF + iDOF] = loadValue; }
    
    // In FSI if Solve is called again with a new set of loads, old loadset needs to be reset
    inline void ResetLoads() {
        for (int iLoad = 0; iLoad < nTotalDOF; iLoad++){
        loadVector[iLoad] = 0.0;}
        structure->ResetForces();}

    inline passivedouble ExtractDisplacements(int iNode, int iDim) {
        return AD::GetValue(structure->GetDisplacement(iNode, iDim));
    }

    inline passivedouble ExtractCoordinate(int iNode, int iDim) {
        return AD::GetValue(node[iNode]->GetCoordinate(iDim));
    }

    inline passivedouble ExtractCoordinate0(int iNode, int iDim) {
        return AD::GetValue(node[iNode]->GetCoordinate0(iDim));
    }

    inline passivedouble ExtractCoordinateOld(int iNode, int iDim) {
        return AD::GetValue(node[iNode]->GetCoordinateOld(iDim));
    }

    inline passivedouble ExtractInitialCoordinates(int iNode, int iDim) {
        return AD::GetValue(structure->GetInitialCoordinates(iNode, iDim));
    }

    inline passivedouble GetInitialCoordinates(int iNode, int iDim) {
        return AD::GetValue(structure->node[iNode]->GetCoordinate(iDim));
    }

    inline void StartRecording(void) { AD::Reset(); AD::StartRecording();}

    inline void RegisterThickness(void) { AD::RegisterInput(thickness);}

    void SetDependencies(void);

    void StopRecording(void);
    
    void StopRecordingWeight(void);
    
    void StopRecordingKSstresses(void);   
    void StopRecordingKSbuckling(void);   
    
    void StopRecordingEA(void); /// DEBUG
    void StopRecordingIzz_b(void); /// DEBUG    
    void StopRecordingNint(void);  /// DEBUG    
    void StopRecordingSigmaBoom(void);  /// DEBUG

    inline passivedouble ExtractLoadGradient(int iNode, int iDOF) {return loadGradient[iNode*nDOF + iDOF];}

    inline passivedouble ExtractGradient_E(void) {return E_grad;}

    inline passivedouble ExtractGradient_Nu(void) {return Nu_grad;}
    
    inline passivedouble ExtractPropGradient(int iPDV){return propGradient[iPDV];}
    
    inline passivedouble ExtractGradient_A(void) {return A_grad;}

    

    inline unsigned long Get_nNodes(void) {return input->Get_nNodes();}

    void StoreDisplacementAdjoint(int iNode, int iDim, passivedouble val_adj);

    void RunRestart(int FSIIter);
    
    void RunRestartLin(int FSIIter);    

    void WriteRestart();

    void ReadRestart();

    void UExtract(std::string line ,int &nNode, double &Ux, double &Uy,
                  double &Uz, double &Urx, double &Ury, double &Urz);

    inline void SetLowVerbosity(void) { verbose = false; structure->SetLowVerbosity();}
    inline void SetHighVerbosity(void) { verbose = true; structure->SetLowVerbosity();}
};