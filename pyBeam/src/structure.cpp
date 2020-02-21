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


#include "../include/structure.h"

CStructure::CStructure(CInput *input, CElement **container_element, CNode **container_node)
{

    // Initialize the pointers
    RBE2    = nullptr;
    element = nullptr;
    node    = nullptr;

    tol_LinSol = input->GetTolerance_LinSol();
    kind_linSol = input->GetKind_LinSol();

    DOF = input->Get_nDOF();

    // Links to the finite-element object
    element = container_element;
    node = container_node;

    nNode = input->Get_nNodes();
    nfem = input->Get_nFEM();
    
    // Get the constrain matrix [NODE_ID DOF_ID]
    Constr_matrix = input->GetConstrMatrix();

    YoungModulus = input->GetYoungModulus_dimensional();
    
    // Resizes and zeros the K matrices
    Ksys.resize(nNode*6,nNode*6);
    Ksys = MatrixXdDiff::Zero(nNode*6,nNode*6);


    U   = VectorXdDiff::Zero(nNode*6);         // Whole system displacements (Cumulative)
    dU  = VectorXdDiff::Zero(nNode*6);         // Incremental system displacements

    U_adj = VectorXdDiff::Zero(nNode*6);       // Whole system displacement adjoint
    
    X  = VectorXdDiff::Zero(nNode*3);          // Current coordinates of the system
    X0 = VectorXdDiff::Zero(nNode*3);          // Initial coordinates of the system


    cross_term = VectorXdDiff::Zero(nNode*6);    // Displacement adjoint cross-term storage

    // Forces nodal Vector
    Fnom     =  VectorXdDiff::Zero(nNode*6);
    Fext     =  VectorXdDiff::Zero(nNode*6);
    Fint     =  VectorXdDiff::Zero(nNode*6);
    Residual =  VectorXdDiff::Zero(nNode*6);

}

CStructure::~CStructure(void)
{
    if (RBE2 != nullptr) delete [] RBE2;  // Pointer to the first RBE2 element
    if (element != nullptr){
        for (int iElem = 0; iElem < nfem; iElem++)
            if (element[iElem] != nullptr) delete [] element[iElem];
        delete [] element;
    }
    if (node != nullptr){
        for (int iNode = 0; iNode < nNode; iNode++)
            if (node[iNode] != nullptr) delete [] node[iNode];
        delete [] node;
    }

}

/**
 * @todo Finalize implementation of RBE2 elements
 * @body It's necessary to finalize and verify the implementation of rigid elements. Before opening a PR, cleanup
 *       debug comments.
 */
/*
//===================================================
//      Assembly RBE2 rigid constraint matrix
//===================================================
void CStructure::AssemblyRigidConstr() 
{

    int i; int j; int iRBE2;
    // Identification of the master DOFS and SLAVE DOFS
    std::vector<int> dofs_all(6*nNode) ;
    for(int i = 0; i < 6*nNode; i++)  { dofs_all[i] = i+1; }   // here dofs go form 1 to 6 differently than for beams (just a convention)
    
    
    std::vector<int> master;
    std::vector<int> master_all(6*nNode);   // full initialization (this vector is going to be for sure smaller)
    master.resize(RBE2[0]->MasterDOFs.size());
    std::vector<int> slave;
    std::vector<int> slave_all;
    slave.resize(RBE2[0]->SlaveDOFs.size());
    //
    std::vector<int>::iterator it;
    
    //finding master dofs and slave dofs
    // EYE here: only in this case DOFS start from 1 instead than from 0. Explained below
    for(iRBE2 = 0; iRBE2 < nRBE2; iRBE2++)
    {
        //std::vector<int> master;
        //master.resize(RBE2[iRBE2]->MasterDOFs.size());
        //VectorXi::Map(&master[0], RBE2[iRBE2]->MasterDOFs.size()) = RBE2[iRBE2]->MasterDOFs;
        //vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
        //master_all.insert( master_all.end(), master.begin(), master.end() );

        VectorXi::Map(&slave[0], RBE2[iRBE2]->SlaveDOFs.size()) = RBE2[iRBE2]->SlaveDOFs;
        slave_all.insert( slave_all.end(), slave.begin(), slave.end() );
        
        
    }
    // all unique slave Dofs for the system (from 1 to 6 for each node )
    //sort( master_all.begin(), master_all.end() );     //master_all.erase( unique( master_all.begin(), master_all.end() ), master_all.end() );
    sort( slave_all.begin(), slave_all.end() ); slave_all.erase( unique( slave_all.begin(), slave_all.end() ), slave_all.end() );
    //Here I do the difference between all the DOFs and the slave DOFS to get only the master Dofs (in the sense all the not slave Dofs)
    it= std::set_difference (dofs_all.begin(), dofs_all.end(), slave_all.begin(), slave_all.end(), master_all.begin());

    // all the master DOFs (from 1 to 6 for each node)
    master_all.resize(it-master_all.begin());
    
    
    // transform the std vector into Eigen matrices occupying the same memory location
    Eigen::Map< VectorXi > master_all_eig(&master_all[0],master_all.size());
    Eigen::Map< VectorXi > slave_all_eig(&slave_all[0],slave_all.size());

    // Evaluation of the full_to_red and red_to_full
    // that's the reason i need dofs from 1 to 6 as 0 represents the slave dofs position
    VectorXi red_to_full = master_all_eig;
    VectorXi full_to_red = VectorXi::Zero(6*nNode);
    // Evaluation of the master_all_eig_red so: the DOF of the master element ordered in the reduced coordinate vector
    VectorXi master_all_eig_red = VectorXi::Zero(master_all.size());
    for (i = 0; i < master_all.size(); i++)
    {
        full_to_red(master_all_eig(i) -1) = i+1 ;
        master_all_eig_red(i) = i+1 ;
    }

    // Initialization of the KRBE matrix and KRBE_ext
    KRBE = MatrixXdDiff::Zero(6*nNode,master_all.size());
    KRBE_ext = MatrixXdDiff::Zero(master_all.size(),master_all.size());
    

    // KRBE assembly
    for (i = 0; i < master_all_eig.size(); i++)
    {
        for (j = 0; j < master_all_eig_red.size(); j++)
        {
            KRBE(master_all_eig(i) -1,master_all_eig_red(i) -1) = 1;
        }
    }

    for(int iRBE2 = 0; iRBE2 < nRBE2; iRBE2++)
    {
        for (i = 0; i < DOF ; i++)
        {
            for (j = 0; j < DOF ; j++)
            {
                KRBE((RBE2[iRBE2]->node_slave-> GeID() -1)*6 +i ,full_to_red((RBE2[iRBE2]->node_master-> GeID()-1)*6+j ) -1 ) = RBE2[iRBE2]->Kinem_matrix(i,j);

            }
        }
        
        /*KRBE_ext = [-Vz*z-Vy*y,       Vy*x,             Vz*x
         *               Vx*y,      -Vz*z - Vx*x,         Vz*y
         *              Vx*z,           Vy*z,          -Vy*y -Vx*x]
         
        //cout << "Axis vector = " << RBE2[iRBE2]->axis_vector.transpose() << endl;
        KRBE_ext.block(full_to_red((RBE2[iRBE2]->node_master-> GeID()-1)*6+3 ) -1,full_to_red((RBE2[iRBE2]->node_master-> GeID()-1)*6+3 ) -1,3,3) <<
                                                                                                                                                     - Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+3 -1)*RBE2[iRBE2]->axis_vector(3 -1) - Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+2 -1)*RBE2[iRBE2]->axis_vector(2 -1) ,
                Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+2 -1)*RBE2[iRBE2]->axis_vector(1 -1) ,
                Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+3 -1)*RBE2[iRBE2]->axis_vector(1 -1) ,
                Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+1 -1)*RBE2[iRBE2]->axis_vector(2 -1) ,
                - Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+3 -1)*RBE2[iRBE2]->axis_vector(3 -1) - Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+1 -1)*RBE2[iRBE2]->axis_vector(1 -1),
                Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+3 -1)*RBE2[iRBE2]->axis_vector(2 -1),
                Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+1 -1)*RBE2[iRBE2]->axis_vector(3 -1),
                Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+2 -1)*RBE2[iRBE2]->axis_vector(3 -1),
                - Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+2 -1)*RBE2[iRBE2]->axis_vector(2 -1) - Fext((RBE2[iRBE2]->node_slave-> GeID()-1)*6+1 -1)*RBE2[iRBE2]->axis_vector(1 -1);
        
        //cout << "KRBE_ext = \n" << KRBE_ext.block(full_to_red((RBE2[iRBE2]->node_master-> GeID()-1)*6+3 ) -1,full_to_red((RBE2[iRBE2]->node_master-> GeID()-1)*6+3 ) -1,3,3) << endl;


    }
    // cout << "KRBE = \n" <<KRBE << endl;
}
*/
//===================================================
//      Add rigid penalty contribution to residual
//===================================================
void CStructure::RigidResidual_FD()
{
    VectorXdDiff Um = VectorXdDiff::Zero(6);
    VectorXdDiff Us = VectorXdDiff::Zero(6);
    VectorXdDiff residual_rigid = VectorXdDiff::Zero(12);
    VectorXdDiff gplus ;
    VectorXdDiff gminus ;
    VectorXdDiff step;
    VectorXdDiff g;
    addouble delta = pow(10,-3);
    int constr_dof_id = 0;
    for (int iRBE2 = 0; iRBE2 < nRBE2; iRBE2++) {
        RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
        //std::cout << "g  = " << RBE2[iRBE2]->g(0)  << " " << RBE2[iRBE2]->g(1)  << " " << RBE2[iRBE2]->g(2)  << " " << RBE2[iRBE2]->g(3)  << " " << RBE2[iRBE2]->g(4)  << " " << RBE2[iRBE2]->g(5)  << " " << std::endl;
        
        for  (int i = 0; i < 12; i++) {
            //std::cout << "i  = " << i << "\n"<<std::endl; 
            // Plus  
            Um = VectorXdDiff::Zero(6);
            Us = VectorXdDiff::Zero(6);
            step = VectorXdDiff::Zero(12);
            step(i) = delta;
            //Master and slave cumulative displacements
            Um = U.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1, 6);
            Us = U.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1, 6);
            Um = Um + step.segment(0,6);
            Us = Us + step.segment(7-1,6);
            //std::cout << "Um  = " << Um(0)  << " " << Um(1)  << " " << Um(2)  << " " << Um(3)  << " " << Um(4)  << " " << Um(5)  << " " << std::endl;
            //std::cout << "Us  = " << Us(0)  << " " << Us(1)  << " " << Us(2)  << " " << Us(3)  << " " << Us(4)  << " " << Us(5)  << " " << std::endl;             
            RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
            gplus = RBE2[iRBE2]->g; 
            //std::cout << "gplus  = " << RBE2[iRBE2]->g(0)  << " " << RBE2[iRBE2]->g(1)  << " " << RBE2[iRBE2]->g(2)  << " " << RBE2[iRBE2]->g(3)  << " " << RBE2[iRBE2]->g(4)  << " " << RBE2[iRBE2]->g(5)  << " " << std::endl;

            
            // minus  
            Um = VectorXdDiff::Zero(6);
            Us = VectorXdDiff::Zero(6);
            step = VectorXdDiff::Zero(12);
            step(i) = -delta;
            //Master and slave cumulative displacements
            Um = U.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1, 6);
            Us = U.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1, 6);
            Um = Um + step.segment(0,6);
            Us = Us + step.segment(7-1,6);
            //std::cout << "Um  = " << Um(0)  << " " << Um(1)  << " " << Um(2)  << " " << Um(3)  << " " << Um(4)  << " " << Um(5)  << " " << std::endl;
            //std::cout << "Us  = " << Us(0)  << " " << Us(1)  << " " << Us(2)  << " " << Us(3)  << " " << Us(4)  << " " << Us(5)  << " " << std::endl;             
            
            RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
            gminus = RBE2[iRBE2]->g;             
            //std::cout << "gminus  = " << RBE2[iRBE2]->g(0)  << " " << RBE2[iRBE2]->g(1)  << " " << RBE2[iRBE2]->g(2)  << " " << RBE2[iRBE2]->g(3)  << " " << RBE2[iRBE2]->g(4)  << " " << RBE2[iRBE2]->g(5)  << " " << std::endl;            
            
            residual_rigid(i) = penalty*( gplus.dot(gplus)*0.5 - gminus.dot(gminus)*0.5)/ (delta*2);
            
        }
            
        Residual.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1,6) = Residual.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1,6) -residual_rigid.segment(1 -1,6);
        Residual.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1,6) = Residual.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1,6) -residual_rigid.segment(7 -1,6);

    }    
    
    // We need to impose again the BC on the residual
    // BC on the residuals
    for (int iii =1; iii<= Constr_matrix.rows(); iii++) {
        constr_dof_id = round(AD::GetValue((Constr_matrix(iii-1,1-1) -1) *6 + Constr_matrix(iii-1,2-1)));
        Residual(constr_dof_id-1) = 0.0;
    }   
    
    std::ofstream myfile;
    myfile.open ("residual_rigid.pyBeam");
    myfile << "residual_rigid= "; myfile <<setprecision(15)<<  residual_rigid/penalty ; myfile << "\n ";
    myfile.close();    
    //debug
    
    //Um = U.segment(RBE2[0]->MasterDOFs(1 - 1) - 1, 6);
    //Us = U.segment(RBE2[0]->SlaveDOFs(1 - 1) - 1, 6);
    //RBE2[0]->EvalConstraintEquation( Um,  Us);
    /*
    std::cout << "g  = " << RBE2[0]->g(0)  << " " << RBE2[0]->g(1)  << " " << RBE2[0]->g(2)  << " " << RBE2[0]->g(3)  << " " << RBE2[0]->g(4)  << " " << RBE2[0]->g(5)  << " " << std::endl;            
    RBE2[0]->EvalJacobian( Um);
    cout << "G = \n" <<RBE2[0]->G.transpose() << endl;  
    residual_rigid = RBE2[0]->G.transpose()*RBE2[0]->g;
    cout << "residual_rigid/penalty = \n" <<residual_rigid << endl;    
    */
    
    //myfile.open ("g.pyBeam", fstream::in | fstream::out | fstream::app);
    //; myfile << RBE2[0]->g ; myfile << "\n "; myfile << "\n ";
    //myfile.close();
}



//===================================================
//      Add rigid penalty contribution to residual
//===================================================
void CStructure::RigidResidual()
{
    
    VectorXdDiff Um = VectorXdDiff::Zero(6); // REVISE
    VectorXdDiff Us = VectorXdDiff::Zero(6);
    VectorXdDiff residual_rigid = VectorXdDiff::Zero(12);
    int constr_dof_id = 0;
    for (int iRBE2 = 0; iRBE2 < nRBE2; iRBE2++) {
        
        // Evaluating constraint equation
        Um = U.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1, 6);
        Us = U.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1, 6);
        
        RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
        RBE2[iRBE2]->EvalConstraintEquation_2( Um,  Us);
        RBE2[iRBE2]->EvalJacobian( Um);
        RBE2[iRBE2]->EvalJacobian_2( Um);
        
        //Evaluating the rigid component of the residual 
        residual_rigid = penalty*RBE2[iRBE2]->G.transpose()*RBE2[iRBE2]->g;
        //cout << "residual_rigid = \n" <<setprecision(20)<<residual_rigid << endl;
        //cout << "Residual = \n" <<setprecision(20)<<Residual << endl;
        Residual.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1,6) = Residual.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1,6) +residual_rigid.segment(1 -1,6);
        Residual.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1,6) = Residual.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1,6) +residual_rigid.segment(7 -1,6);
        
    }
    
    // We need to impose again the BC on the residual
    // BC on the residuals
    for (int iii =1; iii<= Constr_matrix.rows(); iii++) {
        constr_dof_id = round(AD::GetValue((Constr_matrix(iii-1,1-1) -1) *6 + Constr_matrix(iii-1,2-1)));
        Residual(constr_dof_id-1) = 0.0;
    }
    //cout << "Residual+RBE = \n" <<setprecision(20)<<Residual << endl;
}

//===================================================
//      Assembly penalty matrix and vector for rigid constraints (using centered Finite Differences )
//===================================================
void CStructure::AssemblyRigidPenalty_FD()
{
    int n_eq = 6;
    int n_RBEdofs =12;
    int constr_dof_id = 0;
    VectorXdDiff Um;
    VectorXdDiff Us ;    
    VectorXdDiff gplusplus ;
    VectorXdDiff gminusminus ;
    VectorXdDiff gplusminus ;
    VectorXdDiff gminusplus ;
    VectorXdDiff step;
    VectorXdDiff g;
    addouble delta = pow(10,-3);
    
    // Setting to Zero the SYSTEM penalty matrix and residual vector
    K_penal = MatrixXdDiff::Zero(nNode*6,nNode*6);
    MatrixXdDiff Krbe =  MatrixXdDiff::Zero(12,12);
    MatrixXdDiff Krbe2 =  MatrixXdDiff::Zero(12,12);
    //std::cout << "U " <<  " = \n" << setprecision(15)<< U  <<std::endl;
    for (int iRBE2 = 0; iRBE2 < nRBE2; iRBE2++) {   
        Um = U.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1, 6);
        Us = U.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1, 6);
        RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
        g = RBE2[iRBE2]->g;
        for  (int i = 0; i < 12; i++) {
            for  (int j = 0; j < 12; j++) { 
                //std::cout << "i , j = " << i  << " " << j  <<"\n" << std::endl;
                if (i == j)
                {
                    // Plus  
                    Um = VectorXdDiff::Zero(6);
                    Us = VectorXdDiff::Zero(6);
                    step = VectorXdDiff::Zero(12);
                    step(i) = delta;
                    step(j) = delta;
                    //Master and slave cumulative displacements
                    Um = U.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1, 6);
                    Us = U.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1, 6);
                    Um = Um + step.segment(0,6);
                    Us = Us + step.segment(7-1,6);
                    //std::cout << "Um  = " << Um(0)  << " " << Um(1)  << " " << Um(2)  << " " << Um(3)  << " " << Um(4)  << " " << Um(5)  << " " << std::endl;
                    //std::cout << "Us  = " << Us(0)  << " " << Us(1)  << " " << Us(2)  << " " << Us(3)  << " " << Us(4)  << " " << Us(5)  << " " << std::endl;             
                    RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
                    gplusplus = RBE2[iRBE2]->g;
                    //std::cout << "gPlus  = " << RBE2[iRBE2]->g(0)  << " " << RBE2[iRBE2]->g(1)  << " " << RBE2[iRBE2]->g(2)  << " " << RBE2[iRBE2]->g(3)  << " " << RBE2[iRBE2]->g(4)  << " " << RBE2[iRBE2]->g(5)  << " " << std::endl;
                    
                    // Minus  
                    Um = VectorXdDiff::Zero(6);
                    Us = VectorXdDiff::Zero(6);
                    step = VectorXdDiff::Zero(12);
                    step(i) = -delta;
                    step(j) = -delta;
                    //Master and slave cumulative displacements
                    Um = U.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1, 6);
                    Us = U.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1, 6);
                    Um = Um + step.segment(0,6);
                    Us = Us + step.segment(7-1,6);
                    //std::cout << "Um  = " << Um(0)  << " " << Um(1)  << " " << Um(2)  << " " << Um(3)  << " " << Um(4)  << " " << Um(5)  << " " << std::endl;
                    //std::cout << "Us  = " << Us(0)  << " " << Us(1)  << " " << Us(2)  << " " << Us(3)  << " " << Us(4)  << " " << Us(5)  << " " << std::endl;             
                    
                    RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
                    gminusminus = RBE2[iRBE2]->g;    
                    //std::cout << "gMinus  = " << RBE2[iRBE2]->g(0)  << " " << RBE2[iRBE2]->g(1)  << " " << RBE2[iRBE2]->g(2)  << " " << RBE2[iRBE2]->g(3)  << " " << RBE2[iRBE2]->g(4)  << " " << RBE2[iRBE2]->g(5)  << " " << std::endl;                    
                    
                    Krbe(i,j) =  penalty*0.5*( gplusplus.dot(gplusplus) + gminusminus.dot(gminusminus) - 2*g.dot(g))/ (pow(delta,2)); 
                    
                    
                }
                else
                {
                    // PlusPlus  
                    Um = VectorXdDiff::Zero(6);
                    Us = VectorXdDiff::Zero(6);
                    step = VectorXdDiff::Zero(12);
                    step(i) = delta;
                    step(j) = delta;
                    //Master and slave cumulative displacements
                    Um = U.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1, 6);
                    Us = U.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1, 6);
                    Um = Um + step.segment(0,6);
                    Us = Us + step.segment(7-1,6);
                    //std::cout << "Um  = " << Um(0)  << " " << Um(1)  << " " << Um(2)  << " " << Um(3)  << " " << Um(4)  << " " << Um(5)  << " " << std::endl;
                    //std::cout << "Us  = " << Us(0)  << " " << Us(1)  << " " << Us(2)  << " " << Us(3)  << " " << Us(4)  << " " << Us(5)  << " " << std::endl;             
                    
                    RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
                    gplusplus = RBE2[iRBE2]->g;                 
                    //std::cout << "gPlusPlus  = " << RBE2[iRBE2]->g(0)  << " " << RBE2[iRBE2]->g(1)  << " " << RBE2[iRBE2]->g(2)  << " " << RBE2[iRBE2]->g(3)  << " " << RBE2[iRBE2]->g(4)  << " " << RBE2[iRBE2]->g(5)  << " " << std::endl;
                    
                    // Minus  
                    Um = VectorXdDiff::Zero(6);
                    Us = VectorXdDiff::Zero(6);
                    step = VectorXdDiff::Zero(12);
                    step(i) = -delta;
                    step(j) = -delta;
                    //Master and slave cumulative displacements
                    Um = U.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1, 6);
                    Us = U.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1, 6);
                    Um = Um + step.segment(0,6);
                    Us = Us + step.segment(7-1,6);
                    //std::cout << "Um  = " << Um(0)  << " " << Um(1)  << " " << Um(2)  << " " << Um(3)  << " " << Um(4)  << " " << Um(5)  << " " << std::endl;
                    //std::cout << "Us  = " << Us(0)  << " " << Us(1)  << " " << Us(2)  << " " << Us(3)  << " " << Us(4)  << " " << Us(5)  << " " << std::endl;             
                    
                    RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
                    gminusminus = RBE2[iRBE2]->g;     
                    //std::cout << "gMinusMinus  = " << RBE2[iRBE2]->g(0)  << " " << RBE2[iRBE2]->g(1)  << " " << RBE2[iRBE2]->g(2)  << " " << RBE2[iRBE2]->g(3)  << " " << RBE2[iRBE2]->g(4)  << " " << RBE2[iRBE2]->g(5)  << " " << std::endl;
                    
                    // Plus minus  
                    Um = VectorXdDiff::Zero(6);
                    Us = VectorXdDiff::Zero(6);
                    step = VectorXdDiff::Zero(12);
                    step(i) = delta;
                    step(j) = -delta;
                    //Master and slave cumulative displacements
                    Um = U.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1, 6);
                    Us = U.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1, 6);
                    Um = Um + step.segment(0,6);
                    Us = Us + step.segment(7-1,6);
                    //std::cout << "Um  = " << Um(0)  << " " << Um(1)  << " " << Um(2)  << " " << Um(3)  << " " << Um(4)  << " " << Um(5)  << " " << std::endl;
                    //std::cout << "Us  = " << Us(0)  << " " << Us(1)  << " " << Us(2)  << " " << Us(3)  << " " << Us(4)  << " " << Us(5)  << " " << std::endl;             
                    
                    RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
                    gplusminus = RBE2[iRBE2]->g; 
                    //std::cout << "gPlusminus  = " << RBE2[iRBE2]->g(0)  << " " << RBE2[iRBE2]->g(1)  << " " << RBE2[iRBE2]->g(2)  << " " << RBE2[iRBE2]->g(3)  << " " << RBE2[iRBE2]->g(4)  << " " << RBE2[iRBE2]->g(5)  << " " << std::endl;
                    
                    // minus plus  
                    Um = VectorXdDiff::Zero(6);
                    Us = VectorXdDiff::Zero(6);
                    step = VectorXdDiff::Zero(12);
                    step(i) = -delta;
                    step(j) = delta;
                    //Master and slave cumulative displacements
                    Um = U.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1, 6);
                    Us = U.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1, 6);
                    Um = Um + step.segment(0,6);
                    Us = Us + step.segment(7-1,6);
                    //std::cout << "Um  = " << Um(0)  << " " << Um(1)  << " " << Um(2)  << " " << Um(3)  << " " << Um(4)  << " " << Um(5)  << " " << std::endl;
                    //std::cout << "Us  = " << Us(0)  << " " << Us(1)  << " " << Us(2)  << " " << Us(3)  << " " << Us(4)  << " " << Us(5)  << " " << std::endl;             
                    
                    RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
                    gminusplus = RBE2[iRBE2]->g;                    
                    //std::cout << "gminusplus  = " << RBE2[iRBE2]->g(0)  << " " << RBE2[iRBE2]->g(1)  << " " << RBE2[iRBE2]->g(2)  << " " << RBE2[iRBE2]->g(3)  << " " << RBE2[iRBE2]->g(4)  << " " << RBE2[iRBE2]->g(5)  << " " << std::endl;
                    
                    Krbe(i,j) =  penalty*( gplusplus.dot(gplusplus) + gminusminus.dot(gminusminus) - gplusminus.dot(gplusminus) - gminusplus.dot(gminusplus) )*0.5/ (4*pow(delta,2));
                }
                
                
                
                
            }
        }
        
        // Expansion in to the system's tangent matrix
        K_penal.block(RBE2[iRBE2]->MasterDOFs(1 -1) -1,RBE2[iRBE2]->MasterDOFs(1 -1) -1, 6 , 6) = Krbe.block(1 -1, 1 -1, 6, 6) ;
        K_penal.block(RBE2[iRBE2]->MasterDOFs(1 -1) -1,RBE2[iRBE2]->SlaveDOFs(1 -1) -1, 6 , 6) = Krbe.block(1 -1, 7 -1, 6, 6) ;
        K_penal.block(RBE2[iRBE2]->SlaveDOFs(1 -1) -1,RBE2[iRBE2]->MasterDOFs(1 -1) -1, 6 , 6) = Krbe.block(7 -1, 1 -1, 6, 6) ;
        K_penal.block(RBE2[iRBE2]->SlaveDOFs(1 -1) -1,RBE2[iRBE2]->SlaveDOFs(1 -1) -1, 6 , 6) = Krbe.block(7 -1, 7 -1, 6, 6) ;
        
        Ksys  +=    K_penal;
        
    } 
    
    // We need to impose again the boundary conditions
    // Imposing BC
    for (int iii = 1; iii <= Constr_matrix.rows(); iii++) {
        constr_dof_id = round(AD::GetValue((Constr_matrix(iii - 1, 1 - 1) - 1) *6 + Constr_matrix(iii - 1, 2 - 1)));
        Ksys.row(constr_dof_id - 1) = VectorXdDiff::Zero(nNode * 6);
        Ksys.col(constr_dof_id - 1) = VectorXdDiff::Zero(nNode * 6);
        Ksys(constr_dof_id - 1, constr_dof_id - 1) = 1;
    }
    
    //cout << "K_penal = \n" <<K_penal << endl; 
    
    std::ofstream myfile;
    myfile.open ("Krbe.pyBeam");
    myfile <<setprecision(15)<<  Krbe/penalty ; myfile << "\n ";
    myfile.close();
    
    //debug
    
    //Um = U.segment(RBE2[0]->MasterDOFs(1 - 1) - 1, 6);
    //Us = U.segment(RBE2[0]->SlaveDOFs(1 - 1) - 1, 6);
    //RBE2[0]->EvalConstraintEquation( Um,  Us);
    //RBE2[0]->EvalJacobian( Um);
    //Krbe2 =  RBE2[0]->g(0)*RBE2[0]->H_0 + RBE2[0]->g(1)*RBE2[0]->H_1 + RBE2[0]->g(2)*RBE2[0]->H_2 + RBE2[0]->g(3)*RBE2[0]->H_3 + RBE2[0]->g(4)*RBE2[0]->H_4 + RBE2[0]->g(5)*RBE2[0]->H_5;
    //cout << "Krbe2 = \n" <<Krbe2 << endl;
    
      }

//===================================================
//      Assembly penalty matrix and vector for rigid constraints 
//===================================================
void CStructure::AssemblyRigidPenalty()
{
    
    int n_eq = 6;
    int n_RBEdofs =12;
     int constr_dof_id = 0;
    VectorXdDiff Um;
    VectorXdDiff Us ;    
    // Setting to Zero the SYSTEM penalty matrix and residual vector
    K_penal = MatrixXdDiff::Zero(nNode*6,nNode*6);
    
  
    MatrixXdDiff Krbe1 =  MatrixXdDiff::Zero(12,12); // first contribution to the tangent matrix: G^T*G  (single element)
    MatrixXdDiff Krbe2 =  MatrixXdDiff::Zero(12,12); // second contribution to the tangent matrix: sum_i g_i*H_i (single element)



    for (int iRBE2 = 0; iRBE2 < nRBE2; iRBE2++) {
        // EYE here: only in this case DOFS start from 1 instead than from 0. Explained in function AssemblyRigidConstraint

        //Master and slave cumulative displacements
        Um = U.segment(RBE2[iRBE2]->MasterDOFs(1 - 1) - 1, 6);
        Us = U.segment(RBE2[iRBE2]->SlaveDOFs(1 - 1) - 1, 6);

        RBE2[iRBE2]->EvalConstraintEquation( Um,  Us);
        RBE2[iRBE2]->EvalJacobian( Um);
        RBE2[iRBE2]->EvalHessian( Um);
        //// Penalty matrix has 2 contributions:
        
        Krbe1 =  MatrixXdDiff::Zero(12,12);
        Krbe2 =  MatrixXdDiff::Zero(12,12);
        // G^T*G (the Jacobian of the constraint set of equations)
        Krbe1 = RBE2[iRBE2]->G.transpose()*RBE2[iRBE2]->G;
        // sum_i g_i*H_i (from the Hessian of the constraint set of equations)
        //Krbe2 =  RBE2[iRBE2]->g(0)*RBE2[iRBE2]->H_0 + RBE2[iRBE2]->g(1)*RBE2[iRBE2]->H_1 + RBE2[iRBE2]->g(2)*RBE2[iRBE2]->H_2;// + RBE2[iRBE2]->g(3)*RBE2[iRBE2]->H_3 + RBE2[iRBE2]->g(4)*RBE2[iRBE2]->H_4 + RBE2[iRBE2]->g(5)*RBE2[iRBE2]->H_5;
        
        // Expansion in to the system's tangent matrix
        K_penal.block(RBE2[iRBE2]->MasterDOFs(1 -1) -1,RBE2[iRBE2]->MasterDOFs(1 -1) -1, 6 , 6) = Krbe1.block(1 -1, 1 -1, 6, 6) + Krbe2.block(1 -1, 1 -1, 6, 6);
        K_penal.block(RBE2[iRBE2]->MasterDOFs(1 -1) -1,RBE2[iRBE2]->SlaveDOFs(1 -1) -1, 6 , 6) =  Krbe1.block(1 -1, 7 -1, 6, 6) + Krbe2.block(1 -1, 7 -1, 6, 6);
        K_penal.block(RBE2[iRBE2]->SlaveDOFs(1 -1) -1,RBE2[iRBE2]->MasterDOFs(1 -1) -1, 6 , 6) =  Krbe1.block(7 -1, 1 -1, 6, 6) + Krbe2.block(7 -1, 1 -1, 6, 6);
        K_penal.block(RBE2[iRBE2]->SlaveDOFs(1 -1) -1,RBE2[iRBE2]->SlaveDOFs(1 -1) -1, 6 , 6) =   Krbe1.block(7 -1, 7 -1, 6, 6) + Krbe2.block(7 -1, 7 -1, 6, 6);
        
                //cout << "Krbe1 = \n" <<Krbe1 << endl;
                //cout << "Krbe2 = \n" <<Krbe2 << endl;
                //cout << "g = \n" <<RBE2[iRBE2]->g << endl;

 
    }
    
    // Penalty application
    Ksys  +=  - K_penal*penalty;
           
    
    // We need to impose again the boundary conditions
    // Imposing BC
    for (int iii = 1; iii <= Constr_matrix.rows(); iii++) {
        constr_dof_id = round(AD::GetValue((Constr_matrix(iii - 1, 1 - 1) - 1) *6 + Constr_matrix(iii - 1, 2 - 1)));
        Ksys.row(constr_dof_id - 1) = VectorXdDiff::Zero(nNode * 6);
        Ksys.col(constr_dof_id - 1) = VectorXdDiff::Zero(nNode * 6);
        Ksys(constr_dof_id - 1, constr_dof_id - 1) = penalty;
    }
    
    //cout << "RBE2[iRBE2]->G.transpose()*RBE2[iRBE2]->G = \n" <<RBE2[iRBE2]->G.transpose()*RBE2[iRBE2]->G << endl;  
    //cout << "penalty = \n" <<penalty << endl;
    
            std::ofstream myfile;
    myfile.open ("Krbe.pyBeam", fstream::in | fstream::out | fstream::app);
    myfile << "Krbe= "; myfile <<setprecision(15)<<  Krbe1+Krbe2 ; myfile << "\n ";
    myfile.close();
    
}


//===================================================
//      Assembly System Tangent Matrix
//===================================================
/*
 * Requirements:
 *   (0) current length
 *   (1) current Rotation Matrices
 *   (2) current internal stress/deform
 *
 */

void CStructure::AssemblyTang(int iIter)
{

    //std::cout  << " Assembly Tangent Matrix"  << std::endl;
    
    int iii = 0; int dof = 0;   int dof_jjj = 0;   int dof_kkk = 0;
    int constr_dof_id;
    int nodeA_id = 0; int nodeB_id = 0;

    // Intermediate rotation matrix
    MatrixXdDiff Krotated = MatrixXdDiff::Zero(6,6);
    
    // Setting to Zero the stiffness matrix
    Ksys = MatrixXdDiff::Zero(nNode*6,nNode*6);
    
    // Element's contribution to Ktang
    MatrixXdDiff Ktang(12,12);

    /*------------------------------------
     * Elastic contribution to the stiffness matrix
     * Linear elastic + stretch part
     * Computed in ElementTang_Rao
     *------------------------------------*/
    
    for (int id_el=1; id_el<= nfem; id_el++) {
        
        nodeA_id = element[id_el-1]->nodeA->GeID();
        nodeB_id = element[id_el-1]->nodeB->GeID();
        
        //  To evaluate the Tangent the Elastic Matrix needs to be updated
        element[id_el-1]->ElementTang_Rao(iIter, Ktang);

        // Reorganize the element tangent matrix into the global matrix according to the element DOFs
        for (int jjj=1; jjj<= 2; jjj++){
            
            // Determine the global element node id "j"
            if (jjj==1) {dof_jjj = (nodeA_id-1)*6 +1;} else {dof_jjj = (nodeB_id-1)*6 +1;}
            
            for (int kkk=1; kkk<= 2; kkk++){
                
                // Determine the global element node id "k"
                if (kkk==1) {dof_kkk  =  (nodeA_id-1)*6 +1;} else {dof_kkk  =  (nodeB_id-1)*6 +1;}
                
                // Rotates the element's SUBMATRIX tangent
                Krotated = (element[id_el-1]->R * Ktang.block((jjj-1)*6+1 -1,(kkk-1)*6+1 -1,6,6)  ) * element[id_el-1]->R.transpose() ;
                
                // Contribution to the appropriate location of the global matrix
                Ksys.block(dof_jjj-1,dof_kkk-1,6,6) += Krotated;
                
            }
        }
        
    }
    
    /*--------------------------------------------------------
     *    Rigid rotation contribution to the Stiffness Matrix
     * -------------------------------------------------------*/
    
    //EvalSensRot(iIter);
    if (iIter !=0)
    {
        //EvalSensRotFiniteDifferences();
        EvalSensRot();
    }
    /*------------------------------------
     *    Imposing  B.C.
     *------------------------------------*/
    
    // Imposing BC
    for (iii =1; iii<= Constr_matrix.rows(); iii++) {
        constr_dof_id = round(AD::GetValue((Constr_matrix(iii-1,1-1) -1) *6 + Constr_matrix(iii-1,2-1)));
        Ksys.row(constr_dof_id-1) = VectorXdDiff::Zero(nNode*6);
        Ksys.col(constr_dof_id-1) = VectorXdDiff::Zero(nNode*6);
        Ksys(constr_dof_id-1,constr_dof_id-1) = 1.0;
    }
    
    
}

/*===================================================
 *        Evaluate the Sensitivity of Rotation Matrix
 *===================================================*/
/* Given the element's internal forces and the R, evaluates the
 * contribution to the tangent matrix dF/dU = dR/dU*f
 * */

void CStructure::EvalSensRot(){
    
    VectorXdDiff dl_dU =  VectorXdDiff::Zero(12);
    Vector3dDiff e1    = Vector3dDiff::Zero();
    Vector3dDiff e3    = Vector3dDiff::Zero();
    Vector3dDiff e2    = Vector3dDiff::Zero();
    MatrixXdDiff de1 = MatrixXdDiff::Zero(3,12);
    MatrixXdDiff de2 = MatrixXdDiff::Zero(3,12);
    MatrixXdDiff de3 = MatrixXdDiff::Zero(3,12);

    MatrixXdDiff Krot = MatrixXdDiff::Zero(12,12);
    VectorXdDiff fint =  VectorXdDiff::Zero(12);
    
    MatrixXdDiff I = MatrixXdDiff::Identity(3,3);
    
    addouble onetol = 0.0;
    
    MatrixXdDiff de1_part1 = MatrixXdDiff::Zero(3,12);
    
    de1_part1.block(1-1,1-1,3,3) = - I;
    de1_part1.block(1-1,7-1,3,3) =   I;


    VectorXdDiff dU_AB = VectorXdDiff::Zero(12);
    VectorXdDiff  X_AB = VectorXdDiff::Zero(6);

    MatrixXdDiff alpha = MatrixXdDiff::Zero(3,12);
    MatrixXdDiff beta = MatrixXdDiff::Zero(3,12);

    Vector3dDiff p= Vector3dDiff::Zero();

    
    Matrix3dDiff e1_star = Matrix3dDiff::Zero(); // rearrangement of e1 in matrix form
    Matrix3dDiff p_star = Matrix3dDiff::Zero();  // rearrangement of p in matrix form
    
    MatrixXdDiff dp_1 = MatrixXdDiff::Zero(3,12);
    Matrix3dDiff gamma = Matrix3dDiff::Zero();
    Matrix3dDiff dp_block = Matrix3dDiff::Zero();
    Matrix3dDiff e3_star = Matrix3dDiff::Zero(); // rearrangement of e3 in matrix form
    Matrix3dDiff E2 = Matrix3dDiff::Zero();
    Matrix3dDiff E3 = Matrix3dDiff::Zero();
    
    //-------------------------------
    
    VectorXdDiff XbmXa = Vector3dDiff::Zero(3);
    
    int nodeA_id = 0; int nodeB_id = 0;
    
    int dof_jjj = 0; int dof_kkk = 0;

    for (int id_el=1; id_el<= nfem; id_el++) {
        
        // Setting to zero the various coefficients
        //alpha = Vector3dDiff::Zero();
        alpha = MatrixXdDiff::Zero(3,12);
        de1 = MatrixXdDiff::Zero(3,12);
        de3 = MatrixXdDiff::Zero(3,12);
        dp_1 = MatrixXdDiff::Zero(3,12);
        dp_block = Matrix3dDiff::Zero();
        p = Vector3dDiff::Zero();
        p_star = Matrix3dDiff::Zero();
        e1_star = Matrix3dDiff::Zero();
        E2 = Matrix3dDiff::Zero();
        E3 = Matrix3dDiff::Zero();
        beta = MatrixXdDiff::Zero(3,12);
        gamma = Matrix3dDiff::Zero();
        e3_star = Matrix3dDiff::Zero();
        
        nodeA_id = element[id_el-1]->nodeA->GeID();
        nodeB_id = element[id_el-1]->nodeB->GeID();
        
        // Position of the A and B (initial and final) nodes of the element
        // They are already in the global CS (but not the updated final one).
        X_AB.head(3) = X.segment((nodeA_id-1)*3+1 -1,3);
        X_AB.tail(3) = X.segment((nodeB_id-1)*3+1 -1,3);
        
        // Displacements of the A and B (initial and final) nodes of the element
        // They are already in the global CS (but not the updated final one).
        dU_AB.head(6) = dU.segment((nodeA_id-1)*6+1 -1,6);
        dU_AB.tail(6) = dU.segment((nodeB_id-1)*6+1 -1,6);
        
        fint = element[id_el-1]->fint;
        
        e1 = element[id_el-1]->R.block(1-1,1-1,3,1);
        e2 = element[id_el-1]->R.block(1-1,2-1,3,1);
        e3 = element[id_el-1]->R.block(1-1,3-1,3,1);
        
        //===   de1_du ===
        XbmXa = X.segment((nodeB_id-1)*3+1 -1,3) - X.segment((nodeA_id-1)*3+1 -1,3);
        onetol = 1.0/element[id_el-1]->GetCurrent_Length(); //  1/l
        
        dl_dU.head(3)        = -element[id_el-1]->R.block(1-1,1-1,3,1);
        dl_dU.segment(7-1,3) =  element[id_el-1]->R.block(1-1,1-1,3,1);
        
        de1 = (-onetol*onetol)*( XbmXa * dl_dU.transpose());    //
        de1 += onetol*de1_part1;
        
        //===   de3_du === 0
        p = e2;

        
        /* e1_star = [ 0      -e_1(3)  e_1(2)
         *             e_1(3)    0    -e_1(1)
         *            -e_1(2)  e_1(1)   0     ]*/
        e1_star(2-1,1-1) =  e1(3-1);      e1_star(1-1,2-1) =  -e1(3-1);    e1_star(1-1,3-1) =  e1(2-1);
        e1_star(3-1,1-1) =  -e1(2-1);   e1_star(3-1,2-1) =  e1(1-1);      e1_star(2-1,3-1) =  -e1(1-1);

        /* p_star  = [ 0      -p(3)   p(2)
         *             p(3)     0    -p(1)
         *            -p(2)    p(1)    0     ]    alias: e2_star*/

        p_star(2-1,1-1) =  p(3-1);
        p_star(3-1,1-1) =  -p(2-1); //
        p_star(1-1,2-1) =  -p(3-1);
        p_star(3-1,2-1) =  p(1-1); //
        p_star(1-1,3-1) =  p(2-1);
        p_star(2-1,3-1) =  -p(1-1);

        /* e3_star = [ 0      -e_3(3)  e_3(2)
         *             e_3(3)    0    -e_3(1)
         *            -e_3(2)  e_3(1)   0     ]*/
        // rearrangement of e3 in matrix form
        e3_star(2-1,1-1) =  e3(3-1);      e3_star(1-1,2-1) =  -e3(3-1);    e3_star(1-1,3-1) =  e3(2-1);
        e3_star(3-1,1-1) =  -e3(2-1);   e3_star(3-1,2-1) =  e3(1-1);      e3_star(2-1,3-1) =  -e3(1-1);
        
        E2 = e2*e2.transpose();
        
        E3 = e3*e3.transpose();
        
        /* dp_block  = [ 0      e2(3)   -e2(2)
         *             -e2(3)     0    -e2(1)
         *             e2(2)    -e2(1)    0     ]*/
        dp_block(2-1,1-1) =  -e2(3-1);      dp_block(1-1,2-1) =  e2(3-1);    dp_block(1-1,3-1) =  -e2(2-1);
        dp_block(3-1,1-1) =  e2(2-1);   dp_block(3-1,2-1) =  -e2(1-1);      dp_block(2-1,3-1) =  e2(1-1);
        
        dp_1.block(1-1,4-1,3,3) = dp_block;
        dp_1.block(1-1,10-1,3,3) =   dp_block;
        
        dp_1 = dp_1*0.5;
        
        alpha = e1_star*dp_1 - p_star*de1 - E3*(e1_star*dp_1 - p_star*de1);
        
        beta = e3_star*de1 - E2*e3_star*de1;
        
        gamma = I - (e1_star- E3*e1_star)*(- e1_star + E2*e1_star);
        
        
        de3 = alpha;
        
        de2 = (- e1_star + E2*e1_star)* de3 + beta;
        
        // ====== Krot
        
        Krot.block(1-1,1-1,3,12)  =  de1*fint(1-1)  + de2*fint(2-1)  + de3*fint(3-1) ;
        Krot.block(4-1,1-1,3,12)  =  de1*fint(4-1)  + de2*fint(5-1)  + de3*fint(6-1) ;
        Krot.block(7-1,1-1,3,12)  =  de1*fint(7-1)  + de2*fint(8-1)  + de3*fint(9-1) ;
        Krot.block(10-1,1-1,3,12) =  de1*fint(10-1) + de2*fint(11-1) + de3*fint(12-1) ;

        // ================= > insert in the right position
        
        
        for (int jjj=1; jjj<= 2; jjj++) {
            
            if (jjj==1) {dof_jjj  =  (nodeA_id-1)*6 +1;} else {dof_jjj  =  (nodeB_id-1)*6 +1;}
            
            for (int kkk=1; kkk<= 2; kkk++) {
                if (kkk==1) {dof_kkk  =  (nodeA_id-1)*6 +1;} else {dof_kkk  =  (nodeB_id-1)*6 +1;}
                
                Ksys.block(dof_jjj-1,dof_kkk-1,6,6) += Krot.block((jjj-1)*6+1 -1,(kkk-1)*6+1 -1,6,6);
                
            }
            
        }
        
    }

}


/*===================================================
 //      Evaluate the residual
 //===================================================
 This subroutine:
 (a) evaluates the residual
 (b) applies the B.C.
 WARNING: (Fext and Fint need to be updated before)    */

void CStructure::EvalResidual() {
    Residual =  VectorXdDiff::Zero(nNode*6);
    Residual = Fext - Fint;

    int iii = 0; int constr_dof_id = 0;
    
    // BC on the residuals
    for (iii =1; iii<= Constr_matrix.rows(); iii++) {
        constr_dof_id = round(AD::GetValue((Constr_matrix(iii-1,1-1) -1) *6 + Constr_matrix(iii-1,2-1)));
        Residual(constr_dof_id-1) = 0.0;
    }
    
}



/*===================================================
 /      Solve linear static system
 /===================================================
 Solves the linear static problem.
 It needs Ksys and Residual updated*/

void CStructure::SolveLinearStaticSystem(int iIter, std::ofstream &history, int print) {

    bool TapeActive = false;

#ifdef CODI_REVERSE_TYPE

    TapeActive = AD::globalTape.isActive();

    AD::StartExtFunc(false, false);
    unsigned long index = 0;

    for (unsigned long iRes = 0; iRes < Residual.size(); iRes++){
        AD::SetExtFuncIn(Residual(iRes));
    }

    /*--- Stop the recording for the linear solver ---*/

    AD::StopRecording();
#endif
    
    switch(kind_linSol){
        case PartialPivLu:
            dU = Ksys.partialPivLu().solve(Residual); break;
        case FullPivLu:
            dU = Ksys.fullPivLu().solve(Residual); break;
        case HouseholderQr:
            dU = Ksys.householderQr().solve(Residual); break;
        case ColPivHouseholderQr:
            dU = Ksys.colPivHouseholderQr().solve(Residual); break;
        case FullPivHouseholderQr:
            dU = Ksys.fullPivHouseholderQr().solve(Residual); break;
        case LLT:
            dU = Ksys.llt().solve(Residual); break;
        case LDLT:
            dU = Ksys.ldlt().solve(Residual); break;
        default:
            dU = Ksys.fullPivLu().solve(Residual); break;
    }

    if(TapeActive) {

        /*--- Start recording if it was stopped for the linear solver ---*/

        AD::StartRecording();

        for (unsigned long iRes = 0; iRes < Residual.size(); iRes++)
            AD::SetExtFuncOut(dU(iRes));

#ifdef CODI_REVERSE_TYPE
        AD::FuncHelper->addUserData(nNode);
        AD::FuncHelper->addUserData(kind_linSol);
        //     AD::FuncHelper->addUserData(Residual);
        AD::FuncHelper->addUserData(Ksys);

        AD::FuncHelper->addToTape(SolveAdjSys::SolveSys);
#endif

        AD::EndExtFunc();
    }

    addouble relative_error = (Ksys*dU -Residual).norm() / Residual.norm(); // norm() is L2 norm

    if (verbose){
        std::cout.width(17); std::cout << log10(relative_error);
        if (print=1) {history.width(17); history << log10(relative_error); };
    }

    if (relative_error > tol_LinSol)
    {
        std::cout << std:: endl << "Solution of Linear System not precise enough!" << std:: endl;
        std::cout << "Use Keyword TOLERANCE_LINSOL" << std:: endl;
        throw std::exception();
    }
    
    //	Decomposition  	                   Method     Requirements          Speed   Accuracy
    //	PartialPivLU 	             partialPivLu()  Invertible             ++      +
    //	FullPivLU 	                    fullPivLu() 	None                -       +++
    //	HouseholderQR 	             householderQr() 	None                ++      +
    //	ColPivHouseholderQR 	colPivHouseholderQr() 	None                + 	++
    //	FullPivHouseholderQR 	fullPivHouseholderQr() 	None                - 	+++
    //	LLT 	                          llt() 	Positive definite       +++ 	+
    //	LDLT 	                         ldlt() Positive or negative semidefinite 	+++ 	++
    std::ofstream myfile;
    myfile.open ("Residual.pyBeam", fstream::in | fstream::out | fstream::app);
    myfile << "Residual= "; myfile <<setprecision(15)<<  Residual ; myfile << "\n ";
    myfile.close();
}

/*===================================================
 *            Update Coordinates
 * ==================================================
 This member function upates the coordinates XYZ
 (expressed in global reference system) of the finite element nodes.
 
 */

void CStructure::UpdateCoord() {

    //std::cout << "-->  Update Global Coordinates "  << std::endl;
    
    /* We have the X array, we need to add the displacements referred to the pre-last displacement local reference system.
     Thus, this operation need to eb done before the rottion matrix is updated. */
    
    // to correctly update the rotational Dofs of the nodes
    Vector3dDiff U_rot = Vector3dDiff::Zero();
    Vector3dDiff U_rot_new = Vector3dDiff::Zero();
    Vector3dDiff dU_rot = Vector3dDiff::Zero();
    Matrix3dDiff R_U = Matrix3dDiff::Zero();
    Matrix3dDiff R_U_new = Matrix3dDiff::Zero();
    Matrix3dDiff R_dU = Matrix3dDiff::Zero();

    int posX = 1;    // current  position in the X array
    int posU = 1;    // current position in the U array

    VectorXdDiff DX;
    DX = VectorXdDiff::Zero(nNode*3);

 
    // Browsing all the nodes of the current segment
    for (int id_node=1-1; id_node<=nNode-1 ; id_node++) {
        
        DX.segment(posX-1,3) = dU.segment(posU-1,3);  //  *
        X.segment(posX-1,3) += DX.segment(posX-1,3);

        // Update displacements
        U.segment(posU-1,3) += dU.segment(posU-1,3);

        // Update the rotations TODO: check
        
        //The rotation matrix is extracted from the rotational degrees of freedom of each node
        U_rot = U.segment(posU+3 -1,3);
        R_U = Matrix3dDiff::Zero();
        PseudoToRot(U_rot, R_U);

        //same thing is done for the new rotation
        dU_rot = dU.segment(posU+3 -1,3);
        R_dU = Matrix3dDiff::Zero();
        PseudoToRot(dU_rot, R_dU);

        //Rotation is updated
        R_U_new = Matrix3dDiff::Zero();
        R_U_new = R_dU*R_U;

        //and into the vector
        U_rot_new = Vector3dDiff::Zero();
        RotToPseudo(U_rot_new , R_U_new);
        U.segment(posU+3 -1,3) = U_rot_new;
                  
        // Updating the node's coordinates
        for (int iDim=0; iDim < 3; iDim++) {
            node[id_node]->SetCoordinate(iDim, X(posX+iDim-1)) ;
        }

        posX += 3;
        posU += 6;
    }

}

/*===================================================
 *            Restart Coordinates
 * ==================================================
 This member function upates the coordinates XYZ
 (expressed in global reference system) of the finite element nodes.
 
 */

void CStructure::RestartCoord() {

    //std::cout << "-->  Update Global Coordinates "  << std::endl;
    
    /* We have the X array, we need to add the displacements referred to the pre-last displacement local reference system.
     Thus, this operation need to be done before the rotation matrix is updated. */
    

    int posX = 1;    // current  position in the X array
    int posU = 1;    // current position in the U array

    // Browsing all the nodes of the current segment
    for (int id_node=1-1; id_node<=nNode-1 ; id_node++) {

        X.segment(posX-1,3) +=  U.segment(posU-1,3);
        // Updating the node's coordinates
        for (int iDim=0; iDim < 3; iDim++) {
            node[id_node]->SetCoordinate(iDim, X(posX+iDim-1)) ;
        }

        posX += 3;
        posU += 6;
    }
    
}

//===================================================
//     Initialize  Coordinates
//===================================================
/* This member function initialize the coordinates of the fem nodes. It should be changed to a more general function reading the grid ecc ecc
 *
 *
 */
void CStructure::SetCoord0()
{
    std::cout << "-->  Setting the Initial Coordinates "  << std::endl;
    // Here we should refer directly to the node objects
    X0 = VectorXdDiff::Zero(nNode*3);

    int posX = 1;    // current  position in the X array
    int count = 0;   // number of fem upstream the node

    //Browse the nodes    (again this is not related to the number of fem elements)
    for (int id_node=1-1; id_node<= nNode -1; id_node++) {

        for (int iDim=0; iDim < 3; iDim++) {

            // I need the old position of the nodes before the iterative procedure starts
            node[id_node]->SetCoordinateOld(iDim, node[id_node]->GetCoordinate(iDim));

            X0(posX+iDim-1) = node[id_node]->GetCoordinate0(iDim);
        }

        posX += 3;
        count += 1;
    }


}
/*
//===================================================
//     Initialize  Coordinates
//===================================================
/*  This member function initialize the coordinates of the fem nodes.
 *  It should be changed to a more general function reading the grid ecc ecc
 */
void CStructure::InitialCoord()
{
    //std::cout << "---------  Resetting Initial Coordinate Values "  << std::endl;
    // Here we should refer directly to the node objects
    X  = VectorXdDiff::Zero(nNode*3);
    
    int posX = 1;    // current  position in the X array
    int count = 0;   // number of fe upstream the node
    
    //Browse the nodes    (again this is not related to the number of fem elements)
    for (int id_node=1-1; id_node<= nNode -1; id_node++) {

        for (int iDim=0; iDim < 3; iDim++) {
            X(posX+iDim-1) = node[id_node]->GetCoordinate0(iDim);
        }
        
        posX += 3;
        count += 1;
    }
    
    
}

//
//===================================================
//    Update Length of the segment
//===================================================
/*  This member function updates the length of the segment thought as connecting
 * the first and last nodes, and being straight.
 *
 */
void CStructure::UpdateLength()
{
    //std::cout << "-->  Updating Length "  << std::endl;

    int nodeA_id = 0;
    int nodeB_id = 0;
    
    Vector3dDiff Xa = Vector3dDiff::Zero();
    Vector3dDiff Xb = Vector3dDiff::Zero();
    Vector3dDiff temp = Vector3dDiff::Zero();
    
    for (int id_fe=1; id_fe<=nfem; id_fe++)
    {
        nodeA_id = element[id_fe-1]->nodeA->GeID();
        nodeB_id = element[id_fe-1]->nodeB->GeID();
        
        Xa.head(3) = X.segment((nodeA_id-1)*3+1  -1,3);
        Xb.head(3) = X.segment((nodeB_id-1)*3+1  -1,3);
        
        temp = Xb - Xa;
        element[id_fe-1]->SetPrevious_Length();
        element[id_fe-1]->SetCurrent_Length(temp.norm());
    }
    
}

/*===================================================
 *  n Update Rotation Matrix
 *===================================================
 Given the incremental displacement and the previous (cumulative) Rotation Matrix,
 
 (a)  the (cumulative) rotation matrix
 (b)  incremental rotation matrix is updated
 */

void CStructure::UpdateRotationMatrix() {

    //=============   Updating Rotation Matrix   ======================
    
    //std::cout << "-->  Updating Rotation Matrix "  << std::endl;

    VectorXdDiff dU_AB = VectorXdDiff::Zero(12);
    VectorXdDiff  X_AB = VectorXdDiff::Zero(6);
    
    int nodeA_id = 0;
    int nodeB_id = 0;

    int i_fe;
    
    for ( i_fe =1; i_fe<=nfem; i_fe++) {
        
        nodeA_id = element[i_fe-1]->nodeA->GeID();
        nodeB_id = element[i_fe-1]->nodeB->GeID();

        // Position of the A and B (initial and final) nodes of the element
        // They are already in the global CS (but not the updated final one).
        X_AB.head(3) = X.segment((nodeA_id-1)*3+1 -1,3);
        X_AB.tail(3) = X.segment((nodeB_id-1)*3+1 -1,3);

        // Displacements of the A and B (initial and final) nodes of the element
        // They are already in the global CS (but not the updated final one).
        dU_AB.head(6) = dU.segment((nodeA_id-1)*6+1 -1,6);
        dU_AB.tail(6) = dU.segment((nodeB_id-1)*6+1 -1,6);

        // Calling the coordinate update routine
        element[i_fe-1]->EvalRotMat(dU_AB,X_AB);
    }

    i_fe = i_fe -1;
    
}

/*===================================================
 *  n Update Rotation Matrix
 *===================================================
 Given the incremental displacement and the previous (cumulative) Rotation Matrix,
 
 (a)  the (cumulative) rotation matrix
 (b)  incremental rotation matrix is updated
 */

void CStructure::UpdateRotationMatrix_FP() {

    //=============   Updating Rotation Matrix   ======================
    
    //std::cout << "-->  Updating Rotation Matrix "  << std::endl;
    
    VectorXdDiff dU_AB = VectorXdDiff::Zero(12);
    VectorXdDiff  X_AB = VectorXdDiff::Zero(6);
    
    int nodeA_id = 0;
    int nodeB_id = 0;
    int i_fe;
    for (i_fe=1; i_fe<=nfem; i_fe++) {
        
        nodeA_id = element[i_fe-1]->nodeA->GeID();
        nodeB_id = element[i_fe-1]->nodeB->GeID();
        
        // Position of the A and B (initial and final) nodes of the element
        // They are already in the global CS (but not the updated final one).
        X_AB.head(3) = X.segment((nodeA_id-1)*3+1 -1,3);
        X_AB.tail(3) = X.segment((nodeB_id-1)*3+1 -1,3);
        
        // Displacements of the A and B (initial and final) nodes of the element
        // They are already in the global CS (but not the updated final one). 
        dU_AB.head(6) = U.segment((nodeA_id-1)*6+1 -1,6);  // REVISE change name dU -> U
        dU_AB.tail(6) = U.segment((nodeB_id-1)*6+1 -1,6);
        
        // Calling the coordinate update routine
        element[i_fe-1]->EvalRotMat_FP(dU_AB,X_AB);
        
    }
    i_fe = i_fe -1;

}

/*===================================================
 *           INTERNAL FORCES
 *===================================================
 Given the incremental displacement dU,  updated X,l,R
 
 (a) UPDATES the element's elastic matrix
 (b) evaluates the increment of ELEMENT'S INTERNAL ELASTIC DISP and Adds to the CUMULATIVE ELASTIC DISPL.
 (c) evaluates the NODAL INTERNAL FORCE ARRAY
 
 */

void CStructure::UpdateInternalForces()
{
    
    //std::cout << "-->  Updating Internal Forces "   << std::endl;
    
    // dU is the incremental displacement
    // Need to evaluate the displacements in the new reference system.
    // Re is the matrix which rotates from one to the other one.
    
    int nodeA_id = 0;
    int nodeB_id = 0;
    
    // Nodal full Rotation PseudoVectors and Matrices
    Vector3dDiff pseudo_A = Vector3dDiff::Zero();
    Vector3dDiff pseudo_B = Vector3dDiff::Zero();
    
    Vector3dDiff pseudo_A_el = Vector3dDiff::Zero();
    Vector3dDiff pseudo_B_el = Vector3dDiff::Zero();
    
    Matrix3dDiff Rnode_A = Matrix3dDiff::Zero();
    Matrix3dDiff Rnode_B = Matrix3dDiff::Zero();
    
    // Nodal ELASTIC Rotation Matrix
    Matrix3dDiff Rel_A = Matrix3dDiff::Zero();
    Matrix3dDiff Rel_B = Matrix3dDiff::Zero();
    
    Matrix3dDiff Rreduc = Matrix3dDiff::Zero();;
    Matrix3dDiff Rtransp = Matrix3dDiff::Zero();;
    Matrix3dDiff R_rigtransp;

    // Auxiliary matrices for the elastic component
    MatrixXdDiff Na = MatrixXdDiff::Zero(6,6);
    MatrixXdDiff Nb = MatrixXdDiff::Zero(6,6);
    
    int tot_dofs= nNode*6;
    
    // Element's level incremental forces/elastic displ
    VectorXdDiff du_el = VectorXdDiff::Zero(12);
    
    // Nodal vector of internal forces
    // VERY IMPORTANT to reset it to 0 every time
    Fint = VectorXdDiff::Zero(nNode*6);
    
    /*-------------------------------
     //     LOOPING FINITE ELEMENTS
     * -------------------------------*/
    
    int id_fe;
    
    for (id_fe=1;     id_fe <= nfem ; id_fe++) {

        du_el = VectorXdDiff::Zero(12);
        
        nodeA_id = element[id_fe-1]->nodeA->GeID();
        nodeB_id = element[id_fe-1]->nodeB->GeID();

        Vector3dDiff node1_disp = dU.segment((nodeA_id-1)*6+1 -1,3);
        Vector3dDiff node2_disp = dU.segment((nodeB_id-1)*6+1 -1,3);
        
        /*----------------------------
         //      TRANSLATIONAL PART
         * ---------------------------*/

        // Relative displacement of the second node is only along the new axis direction
        du_el(7-1) = element[id_fe-1]->GetCurrent_Length() - element[id_fe-1]->GetPrevious_Length();
        
        /*----------------------------
         *       ROTATIONAL PART
         * ----------------------------*/
        // (a) incremental pseudo-vector is in global CS
        
        // Extracting incremental rotation
        pseudo_A = dU.segment((nodeA_id-1)*6+4 -1,3);
        pseudo_B = dU.segment((nodeB_id-1)*6+4 -1,3);
        
        // (b) transforming nodal pseudo-vector in Rotation/Transformation Matrix
        /* CAREFUL, this rotation does not directly lead from global to nodal triad. But is is
         * an INCREMENTAL rotation given in global coordinates. Which means that, it should be augmented with
         * the rotation from global to old_local.
         * Rnodal = Rnode_A * Rprev  */
        
        PseudoToRot(pseudo_A , Rnode_A);
        PseudoToRot(pseudo_B , Rnode_B);
        
        // (C) Using identity Rnode*Rprev = R*Relastic
        /* Relastic = R'*Rnode_A*Rprev  */
        //
        Rreduc = element[id_fe-1]->R.block(0,0,3,3);
        Rtransp = Rreduc.transpose();
        Rel_A = Rtransp  *  Rnode_A  * element[id_fe-1]->Rprev.block(0,0,3,3);
        Rel_B = Rtransp  *  Rnode_B  * element[id_fe-1]->Rprev.block(0,0,3,3);

        
        // (c) Transforming in pseudo-vector, since infinitesimal (elastic), the components are independent
        RotToPseudo(pseudo_A_el , Rel_A);
        RotToPseudo(pseudo_B_el , Rel_B);

        
        du_el.segment(4 -1,3)  = pseudo_A_el;
        du_el.segment(10 -1,3) = pseudo_B_el;
        
        // Incrementing the cumulative elastic displacements
        // phi is the deformational state vector
        // eps = {    DL,    DTheta ,  Theta_y_el_B ,  Theta_z_el_B , Theta_y_el_A,  Theta_z_el_A)
        element[id_fe-1]->eps(1-1) += du_el( 7-1);
        element[id_fe-1]->eps(2-1) += du_el( 10-1) - du_el( 4-1);
        element[id_fe-1]->eps(3-1) += du_el( 11-1);
        element[id_fe-1]->eps(4-1) += du_el( 12-1);
        element[id_fe-1]->eps(5-1) += du_el( 5-1);
        element[id_fe-1]->eps(6-1) += du_el( 6-1);

        // Constitutive relation between deformational and tensional state
        // phi = tensional state = { N ,  Mt , MBy , MBz , MAy , M_Az }
        
        element[id_fe-1]->phi =  element[id_fe-1]->Kprim*element[id_fe-1]->eps;
        
        Na = MatrixXdDiff::Zero(6,6);
        Nb = MatrixXdDiff::Zero(6,6);
        
        element[id_fe-1]->EvalNaNb(Na , Nb);
        
        // Updating cumulative internal forces
        element[id_fe-1]->fint.segment(1-1,6) =  Na.transpose()*element[id_fe-1]->phi;
        element[id_fe-1]->fint.segment(7-1,6) =  Nb.transpose()*element[id_fe-1]->phi;

        // Contribution to the NODAL Internal Forces ARRAY
        Fint.segment((nodeA_id-1)*6+1 -1,6) +=  element[id_fe-1]->R * element[id_fe-1]->fint.segment(1-1,6);
        Fint.segment((nodeB_id-1)*6+1 -1,6) +=  element[id_fe-1]->R * element[id_fe-1]->fint.segment(7-1,6);

    }

}

void CStructure::InitializeInternalForces()
{
    
    //std::cout << "-->  Updating Internal Forces "   << std::endl;
    
    // dU is the incremental displacement
    // Need to evaluate the displacements in the new reference system.
    // Re is the matrix which rotates from one to the other one.
    
    int nodeA_id = 0;
    int nodeB_id = 0;
    
    // Auxiliary matrices for the elastic component
    MatrixXdDiff Na = MatrixXdDiff::Zero(6,6);
    MatrixXdDiff Nb = MatrixXdDiff::Zero(6,6);

    // Nodal vector of internal forces
    // VERY IMPORTANT to reset it to 0 every time
    Fint = VectorXdDiff::Zero(nNode*6);
    
    /*-------------------------------
     //     LOOPING FINITE ELEMENTS
     * -------------------------------*/
    
    for (int id_fe=1;     id_fe <= nfem ; id_fe++) {

        nodeA_id = element[id_fe-1]->nodeA->GeID();
        nodeB_id = element[id_fe-1]->nodeB->GeID();

        
        // Constitutive relation between deformational and tensional state
        // phi = tensional state = { N ,  Mt , MBy , MBz , MAy , M_Az }
        
        
        element[id_fe-1]->phi =  element[id_fe-1]->Kprim*element[id_fe-1]->eps;
        
        Na = MatrixXdDiff::Zero(6,6);
        Nb = MatrixXdDiff::Zero(6,6);
        
        element[id_fe-1]->EvalNaNb(Na , Nb);
        
        // Updating cumulative internal forces
        element[id_fe-1]->fint.segment(1-1,6) =  Na.transpose()*element[id_fe-1]->phi;
        element[id_fe-1]->fint.segment(7-1,6) =  Nb.transpose()*element[id_fe-1]->phi;

        // Contribution to the NODAL Internal Forces ARRAY
        Fint.segment((nodeA_id-1)*6+1 -1,6) +=  element[id_fe-1]->R * element[id_fe-1]->fint.segment(1-1,6);
        Fint.segment((nodeB_id-1)*6+1 -1,6) +=  element[id_fe-1]->R * element[id_fe-1]->fint.segment(7-1,6);

    }

}


/*===================================================
 *           INTERNAL FORCES FP
 *===================================================
 Given the incremental cumulative displacement U,  updated X,l,R
 
 (a)
 (b) evaluates CUMULATIVE ELASTIC DISPL.
 (c) evaluates the NODAL INTERNAL FORCE ARRAY
 
 */

void CStructure::UpdateInternalForces_FP()
{
    
    //std::cout << "-->  Updating Internal Forces "   << std::endl;
    
    // U is the displacement
    // Need to evaluate the displacements in the new reference system.
    // Re is the matrix which rotates from one to the other one.
    
    int nodeA_id = 0;
    int nodeB_id = 0;
    
    // Nodal full Rotation PseudoVectors and Matrices
    Vector3dDiff pseudo_A = Vector3dDiff::Zero();
    Vector3dDiff pseudo_B = Vector3dDiff::Zero();
    
    Vector3dDiff pseudo_A_el = Vector3dDiff::Zero();
    Vector3dDiff pseudo_B_el = Vector3dDiff::Zero();
    
    Matrix3dDiff Rnode_A = Matrix3dDiff::Zero();
    Matrix3dDiff Rnode_B = Matrix3dDiff::Zero();
    
    // Nodal ELASTIC Rotation Matrix
    Matrix3dDiff Rel_A = Matrix3dDiff::Zero();
    Matrix3dDiff Rel_B = Matrix3dDiff::Zero();
    
    Matrix3dDiff Rreduc = Matrix3dDiff::Zero();;
    Matrix3dDiff Rtransp = Matrix3dDiff::Zero();;

    // Auxiliary matrices for the elastic component
    MatrixXdDiff Na = MatrixXdDiff::Zero(6,6);
    MatrixXdDiff Nb = MatrixXdDiff::Zero(6,6);
    
    int tot_dofs= nNode*6;
    
    // Element's level  forces/elastic displ
    VectorXdDiff u_el = VectorXdDiff::Zero(12);
    
    // Nodal vector of internal forces
    // VERY IMPORTANT to reset it to 0 every time
    Fint = VectorXdDiff::Zero(nNode*6);
    
    VectorXdDiff eps = VectorXdDiff::Zero(6);
    VectorXdDiff phi = VectorXdDiff::Zero(6);
    VectorXdDiff fint = VectorXdDiff::Zero(12);
    
    /*-------------------------------
     //     LOOPING FINITE ELEMENTS
     * -------------------------------*/
    
    int id_fe;
    
    for (id_fe=1;     id_fe <= nfem ; id_fe++) {

        u_el = VectorXdDiff::Zero(12);
        pseudo_A_el = Vector3dDiff::Zero();
        pseudo_B_el = Vector3dDiff::Zero();
        Rel_A = Matrix3dDiff::Zero();
        Rel_B = Matrix3dDiff::Zero();

        
        nodeA_id = element[id_fe-1]->nodeA->GeID();
        nodeB_id = element[id_fe-1]->nodeB->GeID();
        
        /*----------------------------
         //      TRANSLATIONAL PART
         * ---------------------------*/

        // Relative displacement of the second node is only along the new axis direction
        u_el(7-1) = element[id_fe-1]->GetCurrent_Length() - element[id_fe-1]->GetInitial_Length();

        
        /*----------------------------
         *       ROTATIONAL PART
         * ----------------------------*/
        // (a) pseudo-vector is in global CS
        
        // Extracting  rotation
        pseudo_A = U.segment((nodeA_id-1)*6+4 -1,3);
        pseudo_B = U.segment((nodeB_id-1)*6+4 -1,3);
        
        // (b) transforming nodal pseudo-vector in Rotation/Transformation Matrix
        /* CAREFUL, this rotation does not directly lead from global to nodal triad. But is is
         * an  rotation given in global coordinates. Which means that, it should be augmented with
         * the rotation from global to old_local.
         * Rnodal = Rnode_A * Rprev  */
        
        PseudoToRot(pseudo_A , Rnode_A);
        PseudoToRot(pseudo_B , Rnode_B);

        // (C) Using identity Rnode*Rprev = R*Relastic
        /* Relastic = R'*Rnode_A*Rprev  */
        //
        Rreduc = element[id_fe-1]->R.block(0,0,3,3);
        Rtransp = Rreduc.transpose();
        Rel_A = Rtransp  *  Rnode_A * element[id_fe-1]->R0.block(0,0,3,3);
        Rel_B = Rtransp  *  Rnode_B * element[id_fe-1]->R0.block(0,0,3,3) ;
        
        // (c) Transforming in pseudo-vector, since infinitesimal (elastic), the components are independent
        RotToPseudo(pseudo_A_el , Rel_A);
        RotToPseudo(pseudo_B_el , Rel_B);

        u_el.segment(4 -1,3)  = pseudo_A_el;
        u_el.segment(10 -1,3) = pseudo_B_el;
        
        
        // Incrementing the cumulative elastic displacements
        // phi is the deformational state vector
        // eps = {    DL,    DTheta ,  Theta_y_el_B ,  Theta_z_el_B , Theta_y_el_A,  Theta_z_el_A)
        element[id_fe-1]->eps(1-1) = u_el( 7-1);
        element[id_fe-1]->eps(2-1) = u_el( 10-1) - u_el( 4-1);
        element[id_fe-1]->eps(3-1) = u_el( 11-1);
        element[id_fe-1]->eps(4-1) = u_el( 12-1);
        element[id_fe-1]->eps(5-1) = u_el( 5-1);
        element[id_fe-1]->eps(6-1) = u_el( 6-1);

        // Constitutive relation between deformational and tensional state
        // phi = tensional state = { N ,  Mt , MBy , MBz , MAy , M_Az }
        
        element[id_fe-1]->phi =  element[id_fe-1]->Kprim*element[id_fe-1]->eps;
        
        Na = MatrixXdDiff::Zero(6,6);
        Nb = MatrixXdDiff::Zero(6,6);
        
        element[id_fe-1]->EvalNaNb(Na , Nb);
        
        // Updating cumulative internal forces
        element[id_fe-1]->fint.segment(1-1,6) =  Na.transpose()*element[id_fe-1]->phi;
        element[id_fe-1]->fint.segment(7-1,6) =  Nb.transpose()*element[id_fe-1]->phi;

        // Contribution to the NODAL Internal Forces ARRAY
        Fint.segment((nodeA_id-1)*6+1 -1,6) +=  element[id_fe-1]->R * element[id_fe-1]->fint.segment(1-1,6);
        Fint.segment((nodeB_id-1)*6+1 -1,6) +=  element[id_fe-1]->R * element[id_fe-1]->fint.segment(7-1,6);
    }

}

/*
/*===================================================
 *        Evaluate the Sensitivty of Rotation Matrix with finite differences
 *===================================================*/
/* Given the element's internal forces and the R, evaluates the
 * contribution to the tangent matrix dF/dU = dR/dU*f
 * */
/*
void CStructure::EvalSensRotFiniteDifferences(){

    VectorXdDiff dU_AB = VectorXdDiff::Zero(12);
    VectorXdDiff  X_AB = VectorXdDiff::Zero(6);
    Matrix3dDiff R_eps = Matrix3dDiff::Zero();
    
    MatrixXdDiff de1 = MatrixXdDiff::Zero(3,12);
    MatrixXdDiff de2 = MatrixXdDiff::Zero(3,12);
    MatrixXdDiff de3 = MatrixXdDiff::Zero(3,12);
    
    MatrixXdDiff Krot = MatrixXdDiff::Zero(12,12);
    VectorXdDiff fint =  VectorXdDiff::Zero(12);
    
    
    MatrixXdDiff de1_part1 = MatrixXdDiff::Zero(3,12);
    
    de1_part1.block(1-1,1-1,3,3) = - MatrixXdDiff::Identity(3,3);
    de1_part1.block(1-1,7-1,3,3) =   MatrixXdDiff::Identity(3,3);
    
    
    //-------------------------------
    
    
    int nodeA_id = 0; int nodeB_id = 0;
    
    int dof_jjj = 0; int dof_kkk = 0;
    
    // Finite difference for translation DOFs
    addouble fd_t = 1.0e-10;
    // Finite difference for rotational DOFs
    addouble fd_r = 1.0e-10;
    addouble fd;
    
    int ii;
    
    std::ofstream file("./FD_de1_de2_de3.txt");
    std::cout.precision(17);
    
    for (int id_el=1; id_el<= nfem; id_el++) {

        nodeA_id = element[id_el-1]->nodeA->GeID();
        nodeB_id = element[id_el-1]->nodeB->GeID();
        
        fint = element[id_el-1]->fint;
        
        // Position of the A and B (initial and final) nodes of the element
        // They are already in the global CS (but not the updated final one).
        X_AB.head(3) = X.segment((nodeA_id-1)*3+1 -1,3);
        X_AB.tail(3) = X.segment((nodeB_id-1)*3+1 -1,3);
        
        // Displacements of the A and B (initial and final) nodes of the element
        // They are already in the global CS (but not the updated final one).
        dU_AB.head(6) = dU.segment((nodeA_id-1)*6+1 -1,6);
        dU_AB.tail(6) = dU.segment((nodeB_id-1)*6+1 -1,6);
        
        for (ii=1; ii <=12; ii++)
        {
            VectorXdDiff dU_AB_eps = VectorXdDiff::Zero(12);
            if ( ii <4 or ( ii>6 and ii<10) ) {
                fd = fd_t;
            }
            else {
                fd = fd_r;
            }
            
            dU_AB_eps(ii -1) = fd;
            
            element[id_el-1]->EvalRotMatFiniteDifferences( dU_AB_eps, X_AB, R_eps);
            
            de1.block(1-1,ii-1,3,1) =  ( R_eps.block(1-1,1-1,3,1) - element[id_el-1]->R.block(1-1,1-1,3,1) ) / fd;
            de2.block(1-1,ii-1,3,1) =  ( R_eps.block(1-1,2-1,3,1) - element[id_el-1]->R.block(1-1,2-1,3,1) ) / fd;
            de3.block(1-1,ii-1,3,1) =  ( R_eps.block(1-1,3-1,3,1) - element[id_el-1]->R.block(1-1,3-1,3,1) ) / fd;
            
        }
        
        // ====== Krot
        
        Krot.block(1-1,1-1,3,12)  =  de1*fint(1-1)  + de2*fint(2-1)  + de3*fint(3-1) ;
        Krot.block(4-1,1-1,3,12)  =  de1*fint(4-1)  + de2*fint(5-1)  + de3*fint(6-1) ;
        Krot.block(7-1,1-1,3,12)  =  de1*fint(7-1)  + de2*fint(8-1)  + de3*fint(9-1) ;
        Krot.block(10-1,1-1,3,12) =  de1*fint(10-1) + de2*fint(11-1) + de3*fint(12-1) ;
        
        // ================= > insert in the right position
        
        
        for (int jjj=1; jjj<= 2; jjj++) {

            if (jjj==1) {dof_jjj  =  (nodeA_id-1)*6 +1;} else {dof_jjj  =  (nodeB_id-1)*6 +1;}
            
            for (int kkk=1; kkk<= 2; kkk++) {
                if (kkk==1) {dof_kkk  =  (nodeA_id-1)*6 +1;} else {dof_kkk  =  (nodeB_id-1)*6 +1;}
                
                Ksys.block(dof_jjj-1,dof_kkk-1,6,6) += Krot.block((jjj-1)*6+1 -1,(kkk-1)*6+1 -1,6,6);
                
            }
        }
        file  <<  "Element: " << id_el << '\n';
        file  << '\n';
        file  <<  "de1 = \n "<< de1 << '\n';
        file  << '\n';
        file  <<   "de2 = \n "<< de2 << '\n';
        file  << '\n';
        file  <<   "de3 = \n "<< de3 << '\n';
        file  << '\n';
    }
    file.close();
}
*/


