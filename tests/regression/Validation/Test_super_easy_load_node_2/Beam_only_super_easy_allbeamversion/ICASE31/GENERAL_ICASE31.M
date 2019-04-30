

%####################################
%                                   #
%         FEM  PARAMETERS           #
%                                   #
%####################################


% Maximum number of layers. The number of layers used in the properties does not have to exceed this number!
NLMAX = 20;
%
% 1: the tangent matrix (in both linear and nonlinear cases)is forced to be symmetric. This instruction limitates numerical problems.
IFORCESYMMETRY = 0;
%*********************
%   IRIGID:
%     (0): elimination of the rigid body motion using  in the subroutine TRANSDISP
%          the matrix ROTN	 (which is the product between the matrices ROT and Ri (they correspond to
%          the normal ni). This is the original procedure developed by Levy and Gal.
%     (1): elimination of the rigid body motion using in the subroutine TRANSDISP the matrix ROTDEFORM,
%          which is relative to n'. This is my procedure
%     IRIGID is used inside the subroutine PUREROTATIONS
%0: elimination of the rigid body rotations by using Levy's method. IRIGID = 1: Demasi's method is used
%
IRIGID = 2;
%********************
%   KOUTPLANE:
%              (0) Mry and Mrx are set to be zero (see eq.53 of Gal,Levy Computer&Structures2005).
%              (1) Same as (0) with the cancellation of an additional term: er (see eq.53 of Gal,Levy Computer&Structures2005).
%              (2) KGSHELLOP is not modified (this leads to poor performances in Post-buckling regime)
%              (3) All the contribution of the Out of plane matrix (KGSHELLOP = [0]). This behaves better than 2 but is not satisfactory !!!!
%
%  Best performers are:   (1 - Reccommanded), (0), (3), (2) .
%
IKOUTPLANE = 1;  %
%
%**************************************
%   MASS_METHOD:
%               (0)  SEMI-INCONSISTENT
%               (1)  LUMPED
%               (2)  CONSISTENT
%               
ILUMPED = 1;

%**************************************
%   NUMBER of MODES to BE EVALUATED and Written in the output file:
%          
%               
NMODES = 5;
%**************************************
%   NUMBER of MODES to BE PLOTTED
%          
%               
NMODESPLOT = 4;
%**************************************
%   MAGNIFICATION of the MODE (in the plots):
%          
%               
MAGN = 10;  %30 % 50
%**************************************










