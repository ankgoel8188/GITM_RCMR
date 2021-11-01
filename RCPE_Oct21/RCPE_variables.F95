!----------------------------------------------------------------------------
! $Id: RCAC_routines.f95,v 2.1 2018/05/14$
!                       ,v 2.2 2021/10/26$
!
! Author: ANKIT GOEL, University of Michigan, 10/2021
!
! Comments: Contains definitions of the variables used to run RCPE algorithm 
!
!-----------------------------------------------------------------------------


Module RCPE_variables
	IMPLICIT NONE
	INTEGER :: lu, ly, lz, lv
	INTEGER :: ltheta, lphi
	INTEGER :: Nc
	DOUBLE PRECISION :: W_Rtheta, lambda
	
	INTEGER :: nf 					! Order of the FIR filter to be optimized
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Rtheta
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: g_int 						! Integrator. Size lz by 1
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi, phi_kron 				! Regressor phi. Size Nc(lu+lz)+lz by 1. Size lu by lu (Nc(lu+lz)+lz).
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PHI_window					! PHI buffer. pn*lu by ltheta
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u_window 					! Control buffer. lu by pn 
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: z_window 					! Performance buffer. lz by pc
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: N_filt, theta_star 		! Filter and controller gains. lz by nf*lu. ltheta by 1.
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_RLS

END MODULE RCPE_variables



