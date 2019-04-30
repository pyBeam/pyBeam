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

#include "../include/structure.h"

CStructure::CStructure(CInput *input, CElement **element, CNode **container_node)
{
    
    FollFlag = input->Get_FollowerFlag();
    
    DOF = input->Get_nDOF();               
    
    // Links to the finite-element object
    fem = element;
    
    node = container_node;
    
    // Resizes and zeros the M matrices
    nNode = input->Get_nNodes(); 
    nfem = input->Get_nFEM();            
    
    // Get the constrain matrix [NODE_ID DOF_ID]
    Constr_matrix = input->GetConstrMatrix();
    
    
    // Resizes and zeros the K matrices
    Ksys.resize(nNode*6,nNode*6);
    Ksys = MatrixXdDiff::Zero(nNode*6,nNode*6);
    
    M.resize(nNode*6,nNode*6);
    M = MatrixXdDiff::Zero(nNode*6,nNode*6);
    
    
    U  = VectorXdDiff::Zero(nNode*6);         // Whole system displacements (Cumulative)
    dU  = VectorXdDiff::Zero(nNode*6);         // Whole system displacements
    
    X  = VectorXdDiff::Zero(nNode*3);
    
    // Forces nodal Vector
    Ftip     =  Vector3dDiff::Zero();
    Fnom     =  VectorXdDiff::Zero(nNode*6);
    Fext     =  VectorXdDiff::Zero(nNode*6);
    Fint     =  VectorXdDiff::Zero(nNode*6);
    Residual =  VectorXdDiff::Zero(nNode*6);
    
}

CStructure::~CStructure(void)
{
    
}
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
         */
        cout << "Axis vector = " << RBE2[iRBE2]->axis_vector.transpose() << endl;
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

//===================================================
//      Assembly penalty matrix and vector for rigid constraints 
//===================================================
void CStructure::AssemblyRigidPenalty(addouble penalty)
{
    
    /* System matrix penalty
    
    |    0     0 |  {DU_m } = {0}
    |   D_ms  -I |  {DU_s }   {0}
      
    D_ms_rbe =>  DU_s = D_ms * DU_m 
                
       D_ms_rbe =     | I   Kinem_matrix_rbe |
                      | 0         I          |
         
    
    */
    // Setting to Zero the SYSTEM penalty matrix
    K_penal = MatrixXdDiff::Zero(nNode*6,nNode*6);
    
    for(int iRBE2 = 0; iRBE2 < nRBE2; iRBE2++)
    {
        // EYE here: only in this case DOFS start from 1 instead than from 0. Explained in function AssemblyRigidConstraint
        //K_penal.block(RBE2[iRBE2]->SlaveDOFs(1 -1) -1,RBE2[iRBE2]->MasterDOFs(1 -1) -1, 3 , 3) = MatrixXdDiff::Identity(3,3);
        //K_penal.block(RBE2[iRBE2]->SlaveDOFs(4 -1) -1,RBE2[iRBE2]->MasterDOFs(1 -1) -1, 3 , 3) = MatrixXdDiff::Identity(3,3);
        
        //K_penal.block(RBE2[iRBE2]->SlaveDOFs(1 -1) -1,RBE2[iRBE2]->MasterDOFs(4 -1) -1, 3 , 3) = RBE2[iRBE2]->Kinem_matrix.block(1 -1, 4 -1, 3, 3);
        //cout << "DEBUG \n" << endl;
        
        K_penal.block(RBE2[iRBE2]->MasterDOFs(1 -1) -1,RBE2[iRBE2]->MasterDOFs(1 -1) -1, 6 , 6) = RBE2[iRBE2]->Kinem_matrix * RBE2[iRBE2]->Kinem_matrix;
        K_penal.block(RBE2[iRBE2]->MasterDOFs(1 -1) -1,RBE2[iRBE2]->SlaveDOFs(1 -1) -1, 6 , 6) = - RBE2[iRBE2]->Kinem_matrix;
        K_penal.block(RBE2[iRBE2]->SlaveDOFs(1 -1) -1,RBE2[iRBE2]->MasterDOFs(1 -1) -1, 6 , 6) = - RBE2[iRBE2]->Kinem_matrix;
        K_penal.block(RBE2[iRBE2]->SlaveDOFs(1 -1) -1,RBE2[iRBE2]->SlaveDOFs(1 -1) -1, 6 , 6) =  MatrixXdDiff::Identity(6,6);
        
    }
    
    // Penalty application
    K_penal = K_penal*penalty;
    // Penalty vector for residual
    V_penal = VectorXdDiff::Zero(nNode*6);
    V_penal = K_penal*U;
    
    cout << "K_penal = \n" <<K_penal << endl;
    cout << "V_penal = \n" <<V_penal << endl; 
    
}


/***************************************************************
 *
 *         ReadForces
 *
 ***************************************************************
 This subroutine should read the external forces. At the moment,
 it just applies the external force at the TIP
 ***************************************************************/

void CStructure::ReadForces(int nTotalDOF, addouble *loadVector)
{
    
    int iLoad, iDim, index;
    
    for (iLoad = 0; iLoad < nTotalDOF; iLoad++){
        Fnom(iLoad) = loadVector[iLoad];
    }
    
}



/***************************************************************
 *
 *         Update EXternal Forces
 *
 ****************************************************************/
// This subroutine updates the external force.

void CStructure::UpdateExtForces(addouble lambda)
{
    if (FollFlag == 0)
        Fext = lambda* Fnom;
    else if (FollFlag == 1)
    {
        // NOT READY YET!!!
    }
}

//===================================================
//      Assembly System Tangent Matrix
//===================================================
/* WARNING 1: valid only for ordinate,subsequent numeration of beams in 3D
 *
 * WARNING 2: clamped BC are applied
 *
 * WARNING 3: what is needed is:
 *            (0) actual length
 *            (1) actual Rotation Mattrices
 *            (2) actual internal stress/deform
 *
 */
//--------------------------     Evaluates Segment Stiffness FEM
void CStructure::AssemblyTang(int iIter)
{
    //
    
    std::cout  << " Assembly Tangent Matrix"  << std::endl;
    
    int iii = 0; int dof = 0;   int dof_jjj = 0;   int dof_kkk = 0;
    int constr_dof_id;
    int nodeA_id = 0; int nodeB_id = 0;
    MatrixXdDiff Krotated = MatrixXdDiff::Zero(6,6);   // Matrice di appoggio
    
    // Setting to Zero the SYSTEM Stiffness
    Ksys = MatrixXdDiff::Zero(nNode*6,nNode*6);
    
    // Element's contribution to Ktang
    MatrixXdDiff Ktang(12,12);  // 12 is the dimension of the Ktang element level
    /*------------------------------------
     *    Cycle on the finite elements
     *------------------------------------*/
    
    for (int id_el=1; id_el<= nfem; id_el++)
    {
        
        nodeA_id = fem[id_el-1]->nodeA->GeID();
        nodeB_id = fem[id_el-1]->nodeB->GeID();
        
        //  To evaluate the Tangent the Updated Elastic Matrix needs top be Updated
        
        fem[id_el-1]->ElementTang_Rao(iIter, Ktang);       //--> writes the fem[id_el].Ktang
        //std::cout << "Ktang = \n" << Ktang <<std::endl;        
        // For a general approach the element level matrix has to be reorganized into 
        // the global matrix according to the element DOFs
        
        
        //for (int jjj=1; jjj<= 12; jjj+=6)   
        for (int jjj=1; jjj<= 2; jjj++)
        {
            
            if (jjj==1) {dof_jjj  =  (nodeA_id-1)*6 +1;} else {dof_jjj  =  (nodeB_id-1)*6 +1;}           
            for (int kkk=1; kkk<= 2; kkk++)
            {
                if (kkk==1) {dof_kkk  =  (nodeA_id-1)*6 +1;} else {dof_kkk  =  (nodeB_id-1)*6 +1;}          
                
                // Rotates the element's SUBMATRIX tangent
                Krotated = (   fem[id_el-1]->R * Ktang.block((jjj-1)*6+1 -1,(kkk-1)*6+1 -1,6,6)  ) * fem[id_el-1]->R.transpose() ;
                //std::cout << "Krotated = \n" << Krotated <<std::endl; 
                // Contribution in the appropriate SPOT of the SYSTEM TANGENT
                Ksys.block(dof_jjj-1,dof_kkk-1,6,6) += Krotated;
                
            }
        }
        
    }

    /*--------------------------------------------------------
     *    Rigid rotation contribution to the Stiffness Matrix
     * -------------------------------------------------------*/
    
    
    EvalSensRot();

    
    /*------------------------------------
     *    Imposing  B.C.
     *------------------------------------*/
    
    // Imposing BC
    for (iii =1; iii<= Constr_matrix.rows(); iii++)
    {
        constr_dof_id = round(AD::GetValue((Constr_matrix(iii-1,1-1) -1) *6 + Constr_matrix(iii-1,2-1)));
        Ksys.row(constr_dof_id-1) = VectorXdDiff::Zero(nNode*6);
        Ksys.col(constr_dof_id-1) = VectorXdDiff::Zero(nNode*6);
        Ksys(constr_dof_id-1,constr_dof_id-1) = 1.0;          
    }
    
    
}



/*===================================================
 *        Evaluate the Sensitivty of Rotation Matrix
 *===================================================*/
/* Given the element's internal forces and the R, evaluates the
 * contribution to the tangent matrix  dF = dR*f
 * */

void CStructure::EvalSensRot()
{
    VectorXdDiff dl_dU =  VectorXdDiff::Zero(12);
    MatrixXdDiff de1 = MatrixXdDiff::Zero(3,12);
    MatrixXdDiff de2 = MatrixXdDiff::Zero(3,12);
    MatrixXdDiff de3 = MatrixXdDiff::Zero(3,12);
    
    MatrixXdDiff Krot = MatrixXdDiff::Zero(12,12);
    VectorXdDiff fint =  VectorXdDiff::Zero(12);
    
    addouble onetol = 0.0;
    
    MatrixXdDiff de1_part1 = MatrixXdDiff::Zero(3,12);
    
    de1_part1.block(1-1,1-1,3,3) = - MatrixXdDiff::Identity(3,3);
    de1_part1.block(1-1,7-1,3,3) =   MatrixXdDiff::Identity(3,3);
    
    Vector3dDiff de1_i = Vector3dDiff::Zero();
    Vector3dDiff e3    = Vector3dDiff::Zero();
    
    //-------------------------------
    
    
    VectorXdDiff XbmXa = Vector3dDiff::Zero(3);
    
    int nodeA_id = 0; int nodeB_id = 0;    
    
    int dof_jjj = 0; int dof_kkk = 0;
    
    
    for (int id_el=1; id_el<= nfem; id_el++)
    {
        
        nodeA_id = fem[id_el-1]->nodeA->GeID();
        nodeB_id = fem[id_el-1]->nodeB->GeID();        
        
        XbmXa = X.segment((nodeB_id-1)*3+1 -1,3) - X.segment((nodeA_id-1)*3+1 -1,3);
        
        fint = fem[id_el-1]->fint;
        
        onetol =  1.0/fem[id_el-1]->l_act;                 //     1/l
        
        dl_dU.head(3)        = -fem[id_el-1]->R.block(1-1,1-1,3,1);
        dl_dU.segment(7-1,3) =  fem[id_el-1]->R.block(1-1,1-1,3,1);
        
        de1 = (-onetol*onetol)*( XbmXa * dl_dU.transpose());    //
        de1 += onetol*de1_part1;
        
        
        //===   de3_du === 0    TEMPORARY HAS TO BE FIXED!!!!
        de3 = MatrixXdDiff::Zero(3,12);
        
        // de2_du =    de3_du X e1 + e3 X de1_du
        // HNOT EFFICIENT IN THIS WAY!
        for (int i=1; i<= 12; i++)
        {
            de1_i = de1.block(1-1,i-1,3,1);
            e3 = fem[id_el-1]->R.block(1-1,3-1,3,1);
            de2.block(1-1,i-1,3,1) = e3.cross(de1_i);
        }
        
        // ====== Krot
        
        Krot.block(1-1,1-1,3,12)  =  de1*fint(1-1)  + de2*fint(2-1)  + de3*fint(3-1) ;
        Krot.block(4-1,1-1,3,12)  =  de1*fint(4-1)  + de2*fint(5-1)  + de3*fint(6-1) ;
        Krot.block(7-1,1-1,3,12)  =  de1*fint(7-1)  + de2*fint(8-1)  + de3*fint(9-1) ;
        Krot.block(10-1,1-1,3,12) =  de1*fint(10-1) + de2*fint(11-1) + de3*fint(12-1) ;
        
        // ================= > insert in the right position
        
        
        for (int jjj=1; jjj<= 2; jjj++)
        {
            if (jjj==1) {dof_jjj  =  (nodeA_id-1)*6 +1;} else {dof_jjj  =  (nodeB_id-1)*6 +1;} 
            for (int kkk=1; kkk<= 2; kkk++)
            {
                if (kkk==1) {dof_kkk  =  (nodeA_id-1)*6 +1;} else {dof_kkk  =  (nodeB_id-1)*6 +1;} 
                
                Ksys.block(dof_jjj-1,dof_kkk-1,6,6) += Krot.block((jjj-1)*6+1 -1,(kkk-1)*6+1 -1,6,6) ;
                
            }
        }
        
    }  // end loop on the FE
    
}




/*===================================================
 //      Evaluate the residual
 //===================================================
 This subroutine:
 (a) evaluates the residual
 (b) applies the B.C.
 WARNING: (Fext and Fint need to be updated before)    */

void CStructure::EvalResidual()
{
    Residual = Fext - Fint;
   
    std::cout<< "Fext = \n" << Fext.transpose() << std::endl;
    //std::cout<< "Fext_red = \n" << Fext.transpose()*KRBE << std::endl;
    std::cout<< "Fint = \n" << Fint.transpose() << std::endl;
    cout << "Residual = \n" <<Residual.transpose()<< endl;
    
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

void CStructure::SolveLinearStaticSystem(int iIter)
{
    std::cout << "-->  Solving Linear System, "  << std::endl;
    cout << "Ksys = \n" <<Ksys << endl;    
    dU = Ksys.fullPivHouseholderQr().solve(Residual);
    std::cout << "dU (after) = \n" << dU << std::endl;
    addouble relative_error = (Ksys*dU -Residual).norm() / Residual.norm(); // norm() is L2 norm
    //std::cout<< "Ksys = \n" << Ksys << std::endl;
    std::cout<< "Residual = \n" << Residual << std::endl;
    std::cout << "The relative error is:\n" << relative_error << std:: endl;
    if (relative_error > 1.0e-7)
    {
        std::cout << "Solution of Linear System not precise enough!" << std:: endl;
    	throw std::exception();
    }
/*
// Debug
    if (iIter ==0)
    {
    dU =  VectorXdDiff::Zero(18);
    dU(9-1) = pow(10,-6);
    
    std::ofstream file("./dU.dat");
    if (file.is_open())
    {        
        file  <<  dU <<  endl;
    }
    
    std::ofstream file1("./Ktang.dat");
    cout << "Ksys = \n" << Ksys << endl;
    if (file1.is_open())
    {        

        file1  <<  Ksys << endl;
    }    
      
    std::ofstream file3("./Residual_red.dat");
    if (file3.is_open())
    {        

        file3  <<  Residual << endl;
    }    
    
      
     
    }
    
    if (iIter ==1)
    {
    std::ofstream file5("./Residual_iter1.dat");
    if (file5.is_open())
    {        

        file5  <<  Residual << endl;
    }       
          
    }    
 */
    //	Decomposition  	                   Method     Requirements 	Speed 	Accuracy
    //	PartialPivLU 	             partialPivLu()  Invertible 	   ++ 	+
    //	FullPivLU 	                    fullPivLu() 	None 	       - 	+++
    //	HouseholderQR 	             householderQr() 	None 	        ++ 	+
    //	ColPivHouseholderQR 	colPivHouseholderQr() 	None 	         + 	++
    //	FullPivHouseholderQR 	fullPivHouseholderQr() 	None 	         - 	+++
    //	LLT 	                          llt() 	Positive definite  +++ 	+
    //	LDLT 	                         ldlt() Positive or negative semidefinite 	+++ 	++
    
    
}

void CStructure::SolveLinearStaticSystem_RBE2(int iIter)
{ 
    // Debug
    cout << "Ksys = \n" <<Ksys << endl; 
     
    //
    
    
    std::cout << "-->  Reducing Linear System (RBE2...), "  << std::endl;
    Ksys_red = KRBE.transpose()*Ksys*KRBE - KRBE_ext;
    cout << "Ksys_red = \n" <<Ksys_red << endl;
    Residual_red = KRBE.transpose()* Residual; 
    std::cout << "-->  Solving Linear System, "  << std::endl;
    
    std::cout << "Residual red = \n" << Residual_red << endl;
    
    dU_red = Ksys_red.fullPivHouseholderQr().solve(Residual_red);
    //std::cout << "dU_red (after) = \n" << dU_red << std::endl;
    addouble relative_error = (Ksys_red*dU_red -Residual_red).norm() / Residual_red.norm(); // norm() is L2 norm
    //std::cout<< "Ksys = \n" << Ksys << std::endl;
    //std::cout<< "Residual = \n" << Residual.norm() << std::endl;
    std::cout << "The relative error is:\n" << relative_error << std:: endl;
    if (relative_error > 1.0e-7)
    {
        std::cout << "Solution of Linear System not precise enough!" << std:: endl;
    	throw std::exception();
    }
    
    std::cout << "-->  Expanding Linear System (RBE2...), "  << std::endl;
    //Caution, at this point RBE2 slave displacements are still linear
    dU = KRBE*dU_red;
    std::cout << "KRBE.transpose() = \n" << KRBE.transpose() << std::endl;
    std::cout << "KRBE_ext = \n" << KRBE_ext << std::endl;
    std::cout << "dU (after) = \n" << dU << std::endl;
    std::cout << "dU_red (after) = \n" << dU_red << std::endl;
     /*
    // Debug
    if (iIter ==0)
    {
    dU_red =  VectorXdDiff::Zero(18);
    dU_red(8 -1) = pow(10,-6);
    dU = KRBE*dU_red;
    
    cout << "dU = \n" << dU << endl;
    
    std::ofstream file("./dU_red.dat");
    if (file.is_open())
    {        

        file  <<  dU_red  << endl;
    }
    
    std::ofstream file1("./Ktang_red.dat");
    if (file1.is_open())
    {        

        file1  <<  KRBE.transpose()*Ksys*KRBE  << endl;
    }    
    
    std::ofstream file2("./KRBE_ext.dat");
    if (file2.is_open())
    {        

        file2  <<  KRBE_ext  << endl;
    }  
    
    std::ofstream file3("./Residual_red.dat");
    if (file3.is_open())
    {        

        file3  <<  Residual_red  << endl;
    }    
    
    std::ofstream file4("./Ksys_red.dat");
    if (file4.is_open())
    {        

        file4  <<  Ksys_red  << endl;
    }   
    
    std::ofstream file7("./KRBE.dat");
    if (file7.is_open())
    {        

        file7  <<  KRBE  << endl;
    }     
    
    std::ofstream file8("./Ktang.dat");
    if (file7.is_open())
    {        

        file8  <<  Ksys << endl;
    }     
     
    }
    
    if (iIter ==1)
    {
    std::ofstream file5("./Residual_red_iter1.dat");
    if (file5.is_open())
    {        

        file5  <<  Residual_red  << endl;
    }       
    
    std::ofstream file6("./Ktang_red_iter1.dat");
    if (file6.is_open())
    {        

        file6  <<  KRBE.transpose()*Ksys*KRBE  << endl;
    }      
    }
    */
    //VectorXdDiff Fext_red = KRBE.transpose()*Fext;
    //std::cout << "-->  F_ext "  <<Fext << std::endl;

//	Decomposition  	                   Method     Requirements 	Speed 	Accuracy
//	PartialPivLU 	             partialPivLu()  Invertible 	   ++ 	+
//	FullPivLU 	                    fullPivLu() 	None 	       - 	+++
//	HouseholderQR 	             householderQr() 	None 	        ++ 	+
//	ColPivHouseholderQR 	colPivHouseholderQr() 	None 	         + 	++
//	FullPivHouseholderQR 	fullPivHouseholderQr() 	None 	         - 	+++
//	LLT 	                          llt() 	Positive definite  +++ 	+
//	LDLT 	                         ldlt() Positive or negative semidefinite 	+++ 	++

    
}

void CStructure::SolveLinearStaticSystem_RBE2_penalty(int iIter)
{
    std::cout << "-->  Solving Linear System with penalty method for rigid constraints, "  << std::endl;
    //cout << "Ksys = \n" <<Ksys << endl;    
    Ksys = Ksys + K_penal;
    Residual = Residual- V_penal;
    
    dU = Ksys.fullPivHouseholderQr().solve(Residual);
    std::cout << "dU (after) = \n" << dU << std::endl;
    addouble relative_error = (Ksys*dU -Residual).norm() / Residual.norm(); // norm() is L2 norm
    //std::cout<< "Ksys = \n" << Ksys << std::endl;
    std::cout<< "Residual = \n" << Residual << std::endl;
    std::cout << "The relative error is:\n" << relative_error << std:: endl;
    if (relative_error > 1.0e-7)
    {
        std::cout << "Solution of Linear System not precise enough!" << std:: endl;
    	throw std::exception();
    }

}

/*===================================================
 *            Update Coordinates
 * ==================================================
 This member function upates the coordinates XYZ
 (expressed in global reference system) of the finite element nodes.
 
 */

void CStructure::UpdateCoord()
{
    std::cout << "-->  Update Global Coordinates "  << std::endl;
    
    /* We have the X array, we need to add the displacements referred to the pre-last displacement local reference system.
     Thus, this operation need to eb done before the rottion matrix is updated. */
    
    //int n_segdofs;  int n_segXdofs;
    //int indm;
    
    //
    //int indx = 1;    // position in the Xarray of the first node of current segment
    int posX = 1;    // current  position in the X array
    int posU = 1;    // current position in the U array
    
    /*  DEBUG  */
    VectorXdDiff DX;
    
    DX = VectorXdDiff::Zero(nNode*3);
    /**/
    
    // Cumulative displacement update 
    U +=dU;
    
    
    // Browsing all the nodes of the current segment
    for (int i_node=1-1; i_node<=nNode-1 ; i_node++)  // Here we need to check through the connectivity (nfem is not related to the number of nodes)
    {
        
        DX.segment(posX-1,3) = dU.segment(posU-1,3);  //  *
        X.segment(posX-1,3) += DX.segment(posX-1,3);
        
        
        posX += 3;
        posU += 6;
    }
    //std::cout << "X = \n" << X <<std::endl;
    
    ////std::ofstream myfile3 ("./output/echo_dX.out", std::ios_base::out | std::ios_base::app);
    //myfile3 <<  DX << std::endl;
    
}
/*===================================================
 *            Update Coordinates RBE2
 * ==================================================
 */
void CStructure::UpdateCoord_RBE2(int iIter)
{
    
    std::cout << "-->  Update RBE2 Slave coordinates "  << std::endl;    
    
    int iRBE2;
    int idMaster, idSlave;
    VectorXdDiff Master = VectorXdDiff::Zero(3);
    VectorXdDiff Slave = VectorXdDiff::Zero(3);
    VectorXdDiff versor = VectorXdDiff::Zero(3);
    VectorXdDiff Slave_up = VectorXdDiff::Zero(3);
    
    for(iRBE2 = 0; iRBE2 < nRBE2; iRBE2++){
        idMaster = RBE2[iRBE2]->node_master-> GeID();
        idSlave = RBE2[iRBE2]->node_slave-> GeID();
        Master = X.segment((idMaster-1)*3+1 -1,3);
        Slave  = X.segment((idSlave-1)*3+1 -1,3);
        RBE2[iRBE2]->axis_vector_old = RBE2[iRBE2]->axis_vector;
        if (iIter !=0) {
        versor = (Slave - Master)/ (Slave - Master).norm();        
        RBE2[iRBE2]->axis_vector = versor*RBE2[iRBE2]->l_rigid;
        Slave_up = Master + RBE2[iRBE2]->axis_vector;
        X.segment((idSlave-1)*3+1 -1,3) = Slave_up;
        }
        else
        {
            RBE2[iRBE2]->axis_vector = (Slave - Master);
        }
        
    }
    
}

/*===================================================
 *            Update Rigid Constraints (RBE2)
 * ==================================================
 */

void CStructure::UpdateRigidConstr(int iIter)
{
    
    //if (iIter!=1){ 
    UpdateCoord_RBE2(iIter);
    //}    
        
    for(int iRBE2 = 0; iRBE2 < nRBE2; iRBE2++){
        RBE2[iRBE2]->UpdateKinemMatirx();
    }
    
}


//===================================================
//     Initialize  Coordinates
//===================================================
/* This member function initialize the coordinates of the fem nodes. It should be changed to a more general function reading the grid ecc ecc
 *
 *
 */
void CStructure::InitialCoord()
{
    std::cout << "---------  Resetting Initial Coordinate Values "  << std::endl;
    // Here we should refer directly to the node objects
    X  = VectorXdDiff::Zero(nNode*3);
    X0 = VectorXdDiff::Zero(nNode*3);
    
    int posX = 1;    // current  position in the X array
    int count = 0;   // number of fe upstream the node
    
    //Browse the nodes    (again this is not related to the number of fem elements)
    for (int id_node=1-1; id_node<= nNode -1; id_node++)
    {
        
        
        for (int iDim=0; iDim < 3; iDim++) {
            X(posX+iDim-1) = node[id_node]->GetCoordinate(iDim);
            X0(posX+iDim-1) = node[id_node]->GetCoordinate0(iDim);
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
    std::cout << "-->  Updating Length "  << std::endl;
    //
    
    
    int nodeA_id = 0;
    int nodeB_id = 0;
    
    Vector3dDiff Xa = Vector3dDiff::Zero();
    Vector3dDiff Xb = Vector3dDiff::Zero();
    Vector3dDiff temp = Vector3dDiff::Zero();
    
    for (int id_fe=1; id_fe<=nfem; id_fe++)
    {
        nodeA_id = fem[id_fe-1]->nodeA->GeID();
        nodeB_id = fem[id_fe-1]->nodeB->GeID();
        // Remember nodes ID is sequential  node_ini =   (nodeA_id-1)*3+1
        
        Xa.head(3) = X.segment((nodeA_id-1)*3+1  -1,3);
        Xb.head(3) = X.segment((nodeB_id-1)*3+1  -1,3);
        
        temp = Xb - Xa;
        fem[id_fe-1]->l_prev = fem[id_fe-1]->l_act;
        fem[id_fe-1]->l_act = temp.norm();
    }
    
}

/*===================================================
 *  n Update Rotation Matrix
 *===================================================
 Given the incremental displacement and the previous (cumulative) Rotation Matrix,
 
 (a)  the (cumulative) rotation matrix
 (b)  incremental rotation matrix is updated
 */

void CStructure::UpdateRotationMatrix()
{
    
    //To be generalized for the considered connectivity
    
    //=============   Updating Rotation Matrix   ======================
    
    std::cout << "-->  Updating Rotation Matrix "  << std::endl;
    
    // dX_AB is a 6 array, dU_AB is 12 array.
    // First/Last 6 entries are first/last node's dofs current coordinates/displ.  of current finite element
    VectorXdDiff dU_AB = VectorXdDiff::Zero(12);
    VectorXdDiff  X_AB = VectorXdDiff::Zero(6);
    
    
    
    int nodeA_id = 0;
    int nodeB_id = 0;        
    
    // This has to be done for every finite element
    for (int i_fe=1; i_fe<=nfem; i_fe++)
    {
        
        nodeA_id = fem[i_fe-1]->nodeA->GeID();
        nodeB_id = fem[i_fe-1]->nodeB->GeID();
        // Remember nodes ID is sequential  node_ini =   (nodeA_id-1)*3+1 or  (nodeA_id-1)*6+1
        
        X_AB.head(3) = X.segment((nodeA_id-1)*3+1 -1,3)  ;     // They are already in the local CS (but not the updated final one).
        
        X_AB.tail(3) = X.segment((nodeB_id-1)*3+1 -1,3)  ;     // They are already in the local CS (but not the updated final one).
        
        dU_AB.head(6) = dU.segment((nodeA_id-1)*6+1 -1,6)  ;   // They are already in the local CS (but not the updated final one).
	
        dU_AB.tail(6) = dU.segment((nodeB_id-1)*6+1 -1,6)  ;   // They are already in the local CS (but not the updated final one).
        
        fem[i_fe-1]->EvalRotMat(dU_AB,X_AB);     // Calling the coordinate update routine
        
        //cout << "X_AB = \n" << X_AB <<endl;
        //cout << "dU_AB = \n" << dU_AB <<endl;        
    }
    
    
    
#ifdef DEBG
    for (int i_fe=1; i_fe<=nfem; i_fe++)
    {
        //std::ofstream myfile4 ("./output/echo_R_Re.out", std::ios_base::out | std::ios_base::app);
        //		myfile4  << fem[i_fe-1]->Rrig.block(0,0,3,3) << std::endl;
        //		myfile4  << fem[i_fe-1]->R.block(0,0,3,3)  << std::endl;
    }
    
#endif
    
    
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
    
    std::cout << "-->  Updating Internal Forces "   << std::endl;
    
    // dU si the incremental displacement   || Can be faster, since we retrieve dU from xi*phi
    // Need to evaluate the displacements in the new reference system.
    // Re is the matrix whoch rotates from one to the other one.
    
    int nodeA_id = 0;
    int nodeB_id = 0;  
    
    // Nodal full Rotation PseudoVectors and Matrices
    Vector3dDiff pseudo_A = Vector3dDiff::Zero();
    Vector3dDiff pseudo_B = Vector3dDiff::Zero();
    
    Matrix3dDiff Rnode_A = Matrix3dDiff::Zero();
    Matrix3dDiff Rnode_B = Matrix3dDiff::Zero();
    
    // Nodal ELASTIC Rotation Matrix
    Matrix3dDiff Rel_A = Matrix3dDiff::Zero();
    Matrix3dDiff Rel_B = Matrix3dDiff::Zero();
    
    Matrix3dDiff  Rreduc = Matrix3dDiff::Zero();;
    Matrix3dDiff  Rtransp = Matrix3dDiff::Zero();;
    Matrix3dDiff  R_rigtransp;
    
    int tot_dofs= nNode*6;
    
    // Element's level  incremental   forces/elastic displ
    VectorXdDiff duel     = VectorXdDiff::Zero(12);
    
    // Nodal vecotr of internal forces
    Fint = VectorXdDiff::Zero(nNode*6);    // VERY IMPORTANT
    
    /*-------------------------------
     //     LOOPING FINITE ELEMENTS
     * -------------------------------*/
    
    for (int id_fe=1;     id_fe <= nfem ; id_fe++)
    {
        nodeA_id = fem[id_fe-1]->nodeA->GeID();
        nodeB_id = fem[id_fe-1]->nodeB->GeID();
        // Remember nodes ID is sequential  node_ini =   (nodeA_id-1)*3+1 or  (nodeA_id-1)*6+1
        
        Vector3dDiff node1_disp = dU.segment((nodeA_id-1)*6+1 -1,3);
        Vector3dDiff node2_disp = dU.segment((nodeB_id-1)*6+1 -1,3);
        
        /*----------------------------
         //      TRANSLATIONAL PART
         * ---------------------------*/
        // Relative displacement of the second node is only along the new axis direction
        
        duel(7-1) = fem[id_fe-1]->l_act - fem[id_fe-1]->l_prev;
        
        
        /*----------------------------
         *       ROTATIONAL PART
         * ----------------------------*/
        // (a) incremental pseudo-vector is in global CS
        
        // Extracting incremental rotation
        pseudo_A = dU.segment((nodeA_id-1)*6+4 -1,3);
        pseudo_B = dU.segment((nodeB_id-1)*6+4 -1,3);
        
        
        // (b) transforming nodal pesudo-vector in  Rotation/Transformation Matrix
        /*CAREFULL, this rotation does not directly lead from global to nodal triad. But is is
         * an INCREMENTAL rotation given in global coordinates. Which means that, it should be augmented with
         * the rotation from global to old_local.
         * Rnodal = Rnode_A*Rprev  */
        
        PseudoToRot(pseudo_A , Rnode_A);
        PseudoToRot(pseudo_B , Rnode_B);
        
        
        // (C) Using identity      Rnode*Rprev = R*Relastic
        /* Relastic = R'*Rnode_A*Rprev  */
        //
        Rreduc = fem[id_fe-1]->R.block(0,0,3,3);
        Rtransp = Rreduc.transpose();
        Rel_A = Rtransp  *  Rnode_A  * fem[id_fe-1]->Rprev.block(0,0,3,3);
        Rel_B = Rtransp  *  Rnode_B  * fem[id_fe-1]->Rprev.block(0,0,3,3);
        
        // (c) Transforming in pseudo-vector, since infinitesimal (elastic), the components are independent
        RotToPseudo(pseudo_A , Rel_A);
        RotToPseudo(pseudo_B , Rel_B);
        
        
        duel.segment(4 -1,3)  = pseudo_A;
        duel.segment(10 -1,3) = pseudo_B;
        
        
        //
#ifdef DEBG
        //std::ofstream echo_dUel ("./output/echo_dUel.out", std::ios_base::out | std::ios_base::app);
        echo_dUel << duel << std::endl;
#endif
        
        // Icncrementing the cumulative elastic displacements
        // phi is the deformational state vector
        //
        // eps = {    DL,    DTheta ,  Theta_y_el_B ,  Theta_z_el_B , Theta_y_el_A,  Theta_z_el_A)
        
        fem[id_fe-1]->eps(1-1) += duel( 7-1);
        fem[id_fe-1]->eps(2-1) += duel( 10-1) - duel( 4-1);
        fem[id_fe-1]->eps(3-1) += duel( 11-1);
        fem[id_fe-1]->eps(4-1) += duel( 12-1);
        fem[id_fe-1]->eps(5-1) += duel( 5-1);
        fem[id_fe-1]->eps(6-1) += duel( 6-1);
        
        // Constitutive relation between deformational and tensional state
        //
        // phi = tensional state = { N ,  Mt , MBy , MBz , MAy , M_Az }
        
        fem[id_fe-1]->phi =  fem[id_fe-1]->Kprim*fem[id_fe-1]->eps;
        
        MatrixXdDiff Na = MatrixXdDiff::Zero(6,6);
        MatrixXdDiff Nb = MatrixXdDiff::Zero(6,6);
        
        fem[id_fe-1]->EvalNaNb(Na , Nb);
        
        // Updating  cumulative internal forces
        
        fem[id_fe-1]->fint.segment(1-1,6) =  Na.transpose()*fem[id_fe-1]->phi;
        fem[id_fe-1]->fint.segment(7-1,6) =  Nb.transpose()*fem[id_fe-1]->phi;
        
        
        
#ifdef DEBG
        //std::ofstream echo_eps ("./output/echo_eps.out", std::ios_base::out | std::ios_base::app);
        echo_eps << fem[id_fe-1]->eps << std::endl;
        //std::ofstream echo_phi ("./output/echo_phi.out", std::ios_base::out | std::ios_base::app);
        echo_phi << fem[id_fe-1]->phi << std::endl;
        //std::ofstream echo_fint ("./output/echo_fint.out", std::ios_base::out | std::ios_base::app);
        echo_fint << fem[id_fe-1]->fint << std::endl;
#endif
        
        
        // Contribution to the NODAL Internl Forces ARRAY
        
        Fint.segment((nodeA_id-1)*6+1 -1,6)   +=  fem[id_fe-1]->R*  fem[id_fe-1]->fint.segment(1-1,6);
        Fint.segment((nodeB_id-1)*6+1 -1,6) +=  fem[id_fe-1]->R*  fem[id_fe-1]->fint.segment(7-1,6);
        
        
        
    }   // end loop inside the nodes
    
    
    
    //
#ifdef DEBG
    //std::ofstream echo_Fint ("./output/echo_Fint.out", std::ios_base::out | std::ios_base::app);
    echo_Fint << Fint << std::endl;
#endif
    
}
//
//
//
//===================================================
//     TOOLS: Echoes COORDINATES
//===================================================
/*
 * Outputs the Coordinates in Global Ref. System
 */
void CStructure::EchoCoord()
{
    
#ifdef DEBG
    //std::ofstream Xcoord ("./output/echo_Xcoord.out", std::ios_base::out | std::ios_base::app);
    Xcoord <<  X << std::endl;
    
#endif
}


//
/********************************************
 *
 *  Outputs the displacements
 *
 *******************************************/

void CStructure::EchoDisp()
{
    
#ifdef DEBG
    //std::ofstream udisp ("./output/echo_udisp.out", std::ios_base::out | std::ios_base::app);
    udisp <<  dU << std::endl;
#endif
}
//
//===================================================
//      Outputs the Residual
//===================================================

void CStructure::EchoRes()
{
#ifdef DEBG
    //std::ofstream echo_Residual ("./output/echo_Residual.out", std::ios_base::out | std::ios_base::app);            //  Residual in GLOBAL REF
    echo_Residual <<  Residual << std::endl;
#endif
    
}

//===================================================
//      Outputs the External Forces
//===================================================

void CStructure::EchoFext()
{
#ifdef DEBG
    //std::ofstream echo_Fext ("./output/echo_Fext.out", std::ios_base::out | std::ios_base::app);            //  Residual in GLOBAL REF
    echo_Fext <<  Fext << std::endl;
#endif
    
}

//===================================================
//      Outputs the K Matrix
//===================================================

void CStructure::EchoMatrixK()
{
#ifdef DEBG
    //std::ofstream echo_MatrixK ("./output/echo_MatrixK.out", std::ios_base::out | std::ios_base::app);            //  Residual in GLOBAL REF
    
    //	MatrixXdDiff KK = MatrixXdDiff::Zero(12,12);
    //	KK = Ksys.block(1-1,1-1,6,6);
    echo_MatrixK  <<  Ksys.block(1-1,1-1,12,12)  << std::endl;
    
#endif
}

