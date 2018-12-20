/*
 * pyBeam, a Beam Solver
 *
 * Copyright (C) 2018 Ruben Sanchez, Rauno Cavallaro
 * 
 * Developers: Ruben Sanchez (SciComp, TU Kaiserslautern)
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

#include "../include/Rotations.h"


/**********************************
 *
 *        RotToPseudo
 *
 **********************************  */
/* This routine, given a rotation matrix,
transforms it in its pseudo-vector form */

void RotToPseudo(Eigen::Vector3d& pseudo , Eigen::Matrix3d R)
{
	double theta = 0.0;               // Angle of Rotation
	double rho[3] = {0.0,0.0,0.0};

	double fraction = (R.trace()  - 1.0)/2.0;

	if (abs(fraction) >= 1 - 5e-16)   // (FRACTION >= ONE)
		theta = 0.0;
	else
		theta = acos(fraction);

	if (theta>0.0)
	{
		rho[0] = R(3-1,2-1) - R(2-1,3-1);
		rho[1] = R(1-1,3-1) - R(3-1,1-1);
		rho[2] = R(2-1,1-1) - R(1-1,2-1);

		double norm_rho = 0.0 ;
		
		norm_rho = sqrt( pow(rho[0], 2.0) + 
		                 pow(rho[1], 2.0) + 
		                 pow(rho[2], 2.0));

		if (norm_rho <= 1.0e-20)

			pseudo  = Eigen::Vector3d::Zero(3);
		else

			for(unsigned short iVar = 0; iVar < 3; iVar++)
				pseudo(iVar) = theta * rho[iVar]/norm_rho;

		    pseudo(2-1) = (pseudo(2-1));  // tan(pseudo(2-1));
		    pseudo(3-1) = (pseudo(3-1));  // tan(pseudo(3-1));

	}
	else
	{
		pseudo  = Eigen::VectorXd::Zero(3);
	}

}

/**********************************
 *
 *      PseudotoRot
 *
 **********************************  */
/* This routine, given a pseudo-vector,
transforms it in its  rotation matrix form */

void PseudoToRot(Eigen::Vector3d pseudo , Eigen::Matrix3d& R)
{
	//const double pi = 2*acos(0.0);


	Eigen::Vector3d rot(0.0,0.0,0.0);

    pseudo(2-1) = ( pseudo(2-1));  // atan( pseudo(2-1));
    pseudo(3-1) = ( pseudo(3-1));  // atan( pseudo(3-1));

	double theta = pseudo.norm();
	if (theta != 0.0 )
	{
	rot = pseudo/pseudo.norm();
	}

	Eigen::MatrixX3d SkewRot = Eigen::Matrix3d::Zero(3,3);
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

		R.block(0,0,3,3) =  Eigen::MatrixXd::Identity(3,3) +
				theta*(1-theta*theta/6.0 + pow(theta,4)/120.0)  * SkewRot +
				theta*theta*(0.5 - theta*theta/24.0 ) *SkewRot*SkewRot;

	}
	else
	{

		R.block(0,0,3,3) = Eigen::MatrixXd::Identity(3,3) +
				sin(theta) * SkewRot +
				(1 - cos(theta))*SkewRot*SkewRot;
	}


}
