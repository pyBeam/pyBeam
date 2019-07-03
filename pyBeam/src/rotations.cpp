/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Tim Albring, Ruben Sanchez, Rocco Bombardieri, Rauno Cavallaro 
 * 
 * File developers: Rauno Cavallaro (Carlos III University Madrid)
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

#include "../include/rotations.h"


/**********************************
 *
 *        RotToPseudo
 *
 **********************************  */
/* This routine, given a rotation matrix,
transforms it in its pseudo-vector form */

void RotToPseudo(Vector3dDiff& pseudo , Matrix3dDiff R)
{
	addouble theta = 0.0;               // Angle of Rotation
	addouble rho[3] = {0.0,0.0,0.0};

	addouble fraction = (R.trace()  - 1.0)/2.0;

	if (abs(fraction) >= 1 - 5e-16)   // (FRACTION >= ONE)
		theta = 0.0;
	else
		theta = acos(fraction);

	if (theta>0.0)
	{
		rho[0] = R(3-1,2-1) - R(2-1,3-1);
		rho[1] = R(1-1,3-1) - R(3-1,1-1);
		rho[2] = R(2-1,1-1) - R(1-1,2-1);

		addouble norm_rho = 0.0 ;
		
		norm_rho = sqrt( pow(rho[0], 2.0) + 
		                 pow(rho[1], 2.0) + 
		                 pow(rho[2], 2.0));

		if (norm_rho <= 1.0e-20)

			pseudo  = Vector3dDiff::Zero(3);
		else

			for(unsigned short iVar = 0; iVar < 3; iVar++)
				pseudo(iVar) = theta * rho[iVar]/norm_rho;

		    pseudo(2-1) = (pseudo(2-1));  // tan(pseudo(2-1));
		    pseudo(3-1) = (pseudo(3-1));  // tan(pseudo(3-1));

	}
	else
	{
		pseudo  = VectorXdDiff::Zero(3);
	}

}

/**********************************
 *
 *      PseudotoRot
 *
 **********************************  */
/* This routine, given a pseudo-vector,
transforms it in its rotation matrix form */

void PseudoToRot(Vector3dDiff pseudo , Matrix3dDiff& R, int print)
{
	//const addouble pi = 2*acos(0.0);

	Vector3dDiff rot(0.0,0.0,0.0);

    //pseudo(2-1) = ( pseudo(2-1));  // atan( pseudo(2-1));
    //pseudo(3-1) = ( pseudo(3-1));  // atan( pseudo(3-1));

	addouble theta = pseudo.norm();
        
	if (theta != 0.0 )
	{
	rot = pseudo/pseudo.norm();
	}

	Matrix3dDiff SkewRot = Matrix3dDiff::Zero(3,3);
	SkewRot(2-1,1-1) = rot(3-1);    SkewRot(1-1,2-1) = -rot(3-1);
	SkewRot(3-1,1-1) =-rot(2-1);    SkewRot(1-1,3-1) =  rot(2-1);
	SkewRot(3-1,2-1) = rot(1-1);    SkewRot(2-1,3-1) = -rot(1-1);

	/* Rodriguez Matrix
	 *
	 for very small angles need to avoid divisions
	 Using:   sinx = x - x^3/3! + x^5/5!   and
	          cosx = 1 - x^2/2 + x^4/24  */

	if (theta < 1.0e-8)
	{

		R.block(0,0,3,3) =  MatrixXdDiff::Identity(3,3) +
				theta*(1-theta*theta/6.0 + pow(theta,4)/120.0)  * SkewRot +
				theta*theta*(0.5 - theta*theta/24.0 ) *SkewRot*SkewRot;

	}
	else
	{

		R.block(0,0,3,3) = MatrixXdDiff::Identity(3,3) +
				sin(theta) * SkewRot +
				(1 - cos(theta))*SkewRot*SkewRot;
	}


        if (print ==1)
        {
        std::cout << "\nTheta: " << theta << std::endl;      
        std::cout << "\nSkewRot: " << SkewRot << std::endl;      
        std::cout << "\nRot_Rodr: " << R << std::endl; 
           
        }
        
        
}

/**********************************
 *
 *      PseudotoRotDev
 *
 **********************************  */
/* This routine, given a pseudo-vector,
transforms it in the derivatives of the pseudo matrix with respect of each component of the pseudovector */
/*
void PseudoToRotDer(Vector3dDiff pseudo , Matrix3dDiff& dR_1, Matrix3dDiff& dR_2, Matrix3dDiff& dR_3)

{
    addouble aa = pseudo(0);
    addouble bb = pseudo(1);
    addouble cc = pseudo(2);
    addouble mod2 = aa*aa+bb*bb+cc*cc;
    
    addouble mod = sqrt(mod2);
       
    dR_1(0,0) = -aa*sin(mod)*1.0/mod-aa*sin(mod)*(bb*bb+cc*cc)*1.0/pow(mod2,3.0/2.0)-aa*(cos(mod)-1.0)*(bb*bb+cc*cc)*1.0/pow(mod2,2.0)*2.0;
    dR_1(0,1) = -(bb*(cos(mod)-1.0))/(mod2)-(aa*cc*cos(mod))/(mod2)+aa*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+(aa*aa)*bb*sin(mod)*1.0/pow(mod2,3.0/2.0)+(aa*aa)*bb*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_1(0,2) = -(cc*(cos(mod)-1.0))/(mod2)+(aa*bb*cos(mod))/(mod2)-aa*bb*sin(mod)*1.0/pow(mod2,3.0/2.0)+(aa*aa)*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+(aa*aa)*cc*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_1(1,0) = -(bb*(cos(mod)-1.0))/(mod2)+(aa*cc*cos(mod))/(mod2)-aa*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+(aa*aa)*bb*sin(mod)*1.0/pow(mod2,3.0/2.0)+(aa*aa)*bb*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_1(1,1) = -aa*sin(mod)*1.0/mod+(aa*(cos(mod)-1.0)*2.0)/(mod2)-aa*sin(mod)*(aa*aa+cc*cc)*1.0/pow(mod2,3.0/2.0)-aa*(cos(mod)-1.0)*(aa*aa+cc*cc)*1.0/pow(mod2,2.0)*2.0;
    dR_1(1,2) = -sin(mod)*1.0/mod-((aa*aa)*cos(mod))/(mod2)+(aa*aa)*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_1(2,0) = -(cc*(cos(mod)-1.0))/(mod2)-(aa*bb*cos(mod))/(mod2)+aa*bb*sin(mod)*1.0/pow(mod2,3.0/2.0)+(aa*aa)*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+(aa*aa)*cc*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_1(2,1) = sin(mod)*1.0/mod+((aa*aa)*cos(mod))/(mod2)-(aa*aa)*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_1(2,2) = -aa*sin(mod)*1.0/mod+(aa*(cos(mod)-1.0)*2.0)/(mod2)-aa*sin(mod)*(aa*aa+bb*bb)*1.0/pow(mod2,3.0/2.0)-aa*(cos(mod)-1.0)*(aa*aa+bb*bb)*1.0/pow(mod2,2.0)*2.0;'
            
            
            
    dR_2(0,0) = -bb*sin(mod)*1.0/mod+(bb*(cos(mod)-1.0)*2.0)/(mod2)-bb*sin(mod)*(bb*bb+cc*cc)*1.0/pow(mod2,3.0/2.0)-bb*(cos(mod)-1.0)*(bb*bb+cc*cc)*1.0/pow(mod2,2.0)*2.0;
    dR_2(0,1) = -(aa*(cos(mod)-1.0))/(mod2)-(bb*cc*cos(mod))/(mod2)+bb*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*(bb*bb)*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*(bb*bb)*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_2(0,2) = sin(mod)*1.0/mod+((bb*bb)*cos(mod))/(mod2)-(bb*bb)*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_2(1,0) = -(aa*(cos(mod)-1.0))/(mod2)+(bb*cc*cos(mod))/(mod2)-bb*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*(bb*bb)*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*(bb*bb)*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_2(1,1) = -bb*sin(mod)*1.0/mod-bb*sin(mod)*(aa*aa+cc*cc)*1.0/pow(mod2,3.0/2.0)-bb*(cos(mod)-1.0)*(aa*aa+cc*cc)*1.0/pow(mod2,2.0)*2.0;
    dR_2(1,2) = -(cc*(cos(mod)-1.0))/(mod2)-(aa*bb*cos(mod))/(mod2)+aa*bb*sin(mod)*1.0/pow(mod2,3.0/2.0)+(bb*bb)*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+(bb*bb)*cc*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_2(2,0) = -sin(mod)*1.0/mod-((bb*bb)*cos(mod))/(mod2)+(bb*bb)*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_2(2,1) = -(cc*(cos(mod)-1.0))/(mod2)+(aa*bb*cos(mod))/(mod2)-aa*bb*sin(mod)*1.0/pow(mod2,3.0/2.0)+(bb*bb)*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+(bb*bb)*cc*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_2(2,2) = -bb*sin(mod)*1.0/mod+(bb*(cos(mod)-1.0)*2.0)/(mod2)-bb*sin(mod)*(aa*aa+bb*bb)*1.0/pow(mod2,3.0/2.0)-bb*(cos(mod)-1.0)*(aa*aa+bb*bb)*1.0/pow(mod2,2.0)*2.0;'
            
            
            
    dR_3(0,0) = -cc*sin(mod)*1.0/mod+(cc*(cos(mod)-1.0)*2.0)/(mod2)-cc*sin(mod)*(bb*bb+cc*cc)*1.0/pow(mod2,3.0/2.0)-cc*(cos(mod)-1.0)*(bb*bb+cc*cc)*1.0/pow(mod2,2.0)*2.0;
    dR_3(0,1) = -sin(mod)*1.0/mod-((cc*cc)*cos(mod))/(mod2)+(cc*cc)*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_3(0,2) = -(aa*(cos(mod)-1.0))/(mod2)+(bb*cc*cos(mod))/(mod2)-bb*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*(cc*cc)*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*(cc*cc)*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_3(1,0) = sin(mod)*1.0/mod+((cc*cc)*cos(mod))/(mod2)-(cc*cc)*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*bb*cc*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_3(1,1) = -cc*sin(mod)*1.0/mod+(cc*(cos(mod)-1.0)*2.0)/(mod2)-cc*sin(mod)*(aa*aa+cc*cc)*1.0/pow(mod2,3.0/2.0)-cc*(cos(mod)-1.0)*(aa*aa+cc*cc)*1.0/pow(mod2,2.0)*2.0;
    dR_3(1,2) = -(bb*(cos(mod)-1.0))/(mod2)-(aa*cc*cos(mod))/(mod2)+aa*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+bb*(cc*cc)*sin(mod)*1.0/pow(mod2,3.0/2.0)+bb*(cc*cc)*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_3(2,0) = -(aa*(cos(mod)-1.0))/(mod2)-(bb*cc*cos(mod))/(mod2)+bb*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*(cc*cc)*sin(mod)*1.0/pow(mod2,3.0/2.0)+aa*(cc*cc)*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_3(2,1) = -(bb*(cos(mod)-1.0))/(mod2)+(aa*cc*cos(mod))/(mod2)-aa*cc*sin(mod)*1.0/pow(mod2,3.0/2.0)+bb*(cc*cc)*sin(mod)*1.0/pow(mod2,3.0/2.0)+bb*(cc*cc)*(cos(mod)-1.0)*1.0/pow(mod2,2.0)*2.0;
    dR_3(2,2) = -cc*sin(mod)*1.0/mod-cc*sin(mod)*(aa*aa+bb*bb)*1.0/pow(mod2,3.0/2.0)-cc*(cos(mod)-1.0)*(aa*aa+bb*bb)*1.0/pow(mod2,2.0)*2.0;'
            
            
}
*/