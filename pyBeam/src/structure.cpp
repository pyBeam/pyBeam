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

CStructure::CStructure(CInput *input, CElement **element)
{

	FollFlag = input->Get_FollowerFlag();

	DOF = input->Get_nDOF();   // To be replaced             

	// Links to the finite-element object
	fem = element;

	// Resizes and zeros the M matrices
	nfem = input->Get_nFEM();             // to be replaced
        
        // Get the constrain matrix [NODE_ID DOF_ID]
        Constr_matrix = input->GetConstrMatrix();
        

	// Resizes and zeros the K matrices
	Ksys.resize((nfem+1)*6,(nfem+1)*6);
	Ksys = MatrixXdDiff::Zero((nfem+1)*6,(nfem+1)*6);

	M.resize((nfem+1)*6,(nfem+1)*6);
	M = MatrixXdDiff::Zero((nfem+1)*6,(nfem+1)*6);

	dU  = VectorXdDiff::Zero((nfem+1)*6);         // Whole system displacements

	X  = VectorXdDiff::Zero((nfem+1)*3);

	// Forces nodal Vector
	Ftip     =  Vector3dDiff::Zero();
        Fnom     =  VectorXdDiff::Zero((nfem+1)*6);
	Fext     =  VectorXdDiff::Zero((nfem+1)*6);
	Fint     =  VectorXdDiff::Zero((nfem+1)*6);
	Residual =  VectorXdDiff::Zero((nfem+1)*6);
			
}

CStructure::~CStructure(void)
{

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
void CStructure::AssemblyTang()
{
	//

	std::cout  << " Assembly Tangent Matrix"  << std::endl;

	int iii = 0; int dof = 0;   int dof_jjj = 0;   int dof_kkk = 0;
        int constr_dof_id;
	MatrixXdDiff Krotated = MatrixXdDiff::Zero(6,6);   // Matrice di appoggio

	// Setting to Zero the SYSTEM Stiffness
	Ksys = MatrixXdDiff::Zero((nfem+1)*6,(nfem+1)*6);

	// Element's contribution to Ktang
	MatrixXdDiff Ktang(12,12);  // 12 is the dimension of the Ktang element level
	/*------------------------------------
	 *    Cycle on the finite elements
	 *------------------------------------*/

	for (int id_el=1; id_el<= nfem; id_el++)
	{

		//  To evaluate the Tangent the Updated Elastic Matrix needs top be Updated

		fem[id_el-1]->ElementTang_Rao(Ktang);       //--> writes the fem[id_el].Ktang

		dof = (id_el-1)*6 ;   // first-1 dof

		for (int jjj=1; jjj<= 12; jjj+=6)   // Browses 1/7
		{
			dof_jjj  = dof + jjj;              // prev+1 / prev + 7
			for (int kkk=1; kkk<= 12; kkk+=6)
			{
				dof_kkk = dof + kkk;            // prev+1 / prev + 7

				// Rotates the element's SUBMATRIX tangent
				Krotated = (   fem[id_el-1]->R * Ktang.block(jjj-1,kkk-1,6,6)  ) * fem[id_el-1]->R.transpose() ;

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
            constr_dof_id = ( Constr_matrix(iii-1,1-1) -1) *6 + Constr_matrix(iii-1,2-1);
            Ksys.row(constr_dof_id-1) = VectorXdDiff::Zero((nfem+1)*6);
            Ksys.col(constr_dof_id-1) = VectorXdDiff::Zero((nfem+1)*6);
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

	int idofXini = 1;
	VectorXdDiff XbmXa = Vector3dDiff::Zero(3);

	int dof_jjj = 0; int dof_kkk = 0;


	for (int id_el=1; id_el<= nfem; id_el++)
	{

		XbmXa = X.segment(idofXini+3-1,3) - X.segment(idofXini-1,3);

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

		int dof = (id_el-1)*6 ;   // first-1 dof

		for (int jjj=1; jjj<= 12; jjj+=6)
		{
			dof_jjj  = dof + jjj;
			for (int kkk=1; kkk<= 12; kkk+=6)
			{
				dof_kkk = dof + kkk;


				Ksys.block(dof_jjj-1,dof_kkk-1,6,6) += Krot.block(jjj-1,kkk-1,6,6) ;

			}
		}

		idofXini+=3;

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

        int iii = 0; int constr_dof_id = 0;
        
	// BC on the residuals
        for (iii =1; iii<= Constr_matrix.rows(); iii++) {
            constr_dof_id = ( Constr_matrix(iii-1,1-1) -1) *6 + Constr_matrix(iii-1,2-1);
            Residual(constr_dof_id-1) = 0.0;
        }
        
}



/*===================================================
/      Solve linear static system
/===================================================
Solves the linear static problem.
It needs Ksys and Residual updated*/

void CStructure::SolveLinearStaticSystem()
{

	std::cout << "-->  Solving Linear System, "  << std::endl;
	
	dU = Ksys.fullPivHouseholderQr().solve(Residual);

	addouble relative_error = (Ksys*dU -Residual).norm() / Residual.norm(); // norm() is L2 norm
	std::cout << "The relative error is:\n" << relative_error << std:: endl;
    if (relative_error > 1.0e-7)
    {
    	std::cout << "Solution of Linear SYstem not precise enough!" << std:: endl;
    	throw std::exception();
    }

//	Decomposition  	                   Method     Requirements 	Speed 	Accuracy
//	PartialPivLU 	             partialPivLu()  Invertible 	   ++ 	+
//	FullPivLU 	                    fullPivLu() 	None 	       - 	+++
//	HouseholderQR 	             householderQr() 	None 	        ++ 	+
//	ColPivHouseholderQR 	colPivHouseholderQr() 	None 	         + 	++
//	FullPivHouseholderQR 	fullPivHouseholderQr() 	None 	         - 	+++
//	LLT 	                          llt() 	Positive definite  +++ 	+
//	LDLT 	                         ldlt() Positive or negative semidefinite 	+++ 	++


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

	int n_segdofs;  int n_segXdofs;
	int indm;

	//
	int indx = 1;    // position in the Xarray of the first node of current segment
	int posX = 1;    // current  position in the X array
	int posU = 1;    // current position in the U array

	/*  DEBUG  */
	VectorXdDiff DX;
	DX  = VectorXdDiff::Zero((nfem+1)*3);
	/**/

	// Browsing all the fem nodes of the current segment
	for (int i_node=1-1; i_node<=(nfem+1)-1 ; i_node++)  // Here we need to check through the connectivity (nfem is not related to the number of nodes)
	{

		DX.segment(posX-1,3) = dU.segment(posU-1,3);  //  *
		X.segment(posX-1,3) += DX.segment(posX-1,3);


		posX += 3;
		posU += 6;
	}


	////std::ofstream myfile3 ("./output/echo_dX.out", std::ios_base::out | std::ios_base::app);
	//myfile3 <<  DX << std::endl;

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

	X  = VectorXdDiff::Zero((nfem+1)*3);
	X0 = VectorXdDiff::Zero((nfem+1)*3);

	addouble le = fem[0]->le;

	int posX = 1;    // current  position in the X array
	int count = 0;   // number of fe upstream the node

	//Browse the nodes
	for (int id_node=1-1; id_node<= nfem + 1 -1; id_node++)
	{

		X(posX-1) = le*count;   // careful this accumulates the error
		X0(posX-1) = le*count;

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

	int node_ini = 1;


	Vector3dDiff Xa = Vector3dDiff::Zero();
	Vector3dDiff Xb = Vector3dDiff::Zero();
	Vector3dDiff temp = Vector3dDiff::Zero();

	for (int id_fe=1; id_fe<=nfem; id_fe++)
	{

		Xa.head(3) = X.segment(node_ini-1,3);
		Xb.head(3) = X.segment(node_ini+3-1,3);

		temp = Xb - Xa;
		fem[id_fe-1]->l_prev = fem[id_fe-1]->l_act;
		fem[id_fe-1]->l_act = temp.norm();

		node_ini += 3;

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

	int indu = 1;  int indX =1;

	// This has to be done for every finite element
	for (int i_fe=1; i_fe<=nfem; i_fe++)
	{
		X_AB.head(3) = X.segment(indX-1,3)  ;     // They are already in the local CS (but not the updated final one).
		indX += 3;
		X_AB.tail(3) = X.segment(indX-1,3)  ;     // They are already in the local CS (but not the updated final one).

		dU_AB.head(6) = dU.segment(indu-1,6)  ;   // They are already in the local CS (but not the updated final one).
		indu += 6;
		dU_AB.tail(6) = dU.segment(indu-1,6)  ;   // They are already in the local CS (but not the updated final one).

		fem[i_fe-1]->EvalRotMat(dU_AB,X_AB);     // Calling the coordinate update routine
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

	int dof_ini = 1;

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
	int tot_dofs= (nfem+1)*6;

	// Element's level  incremental   forces/elastic displ
	VectorXdDiff duel     = VectorXdDiff::Zero(12);

	// Nodal vecotr of internal forces
	Fint = VectorXdDiff::Zero((nfem+1)*6);    // VERY IMPORTANT

	/*-------------------------------
	//     LOOPING FINITE ELEMENTS
	 * -------------------------------*/

	for (int id_fe=1;     id_fe <= nfem ; id_fe++)
	{

		Vector3dDiff node1_disp = dU.segment(dof_ini-1,3);
		Vector3dDiff node2_disp = dU.segment(dof_ini+6-1,3);

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
		pseudo_A = dU.segment(dof_ini+3 -1,3);
		pseudo_B = dU.segment(dof_ini+9 -1,3);


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

		Fint.segment(dof_ini-1,6)   +=  fem[id_fe-1]->R*  fem[id_fe-1]->fint.segment(1-1,6);
		Fint.segment(dof_ini+6-1,6) +=  fem[id_fe-1]->R*  fem[id_fe-1]->fint.segment(7-1,6);

		dof_ini += 6;   // increase index in the complete U vector

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

