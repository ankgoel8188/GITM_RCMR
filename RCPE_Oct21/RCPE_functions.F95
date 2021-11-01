!----------------------------------------------------------------------------
! $Id: RCAC_routines.f95,v 2.1 2018/05/14$
!                       ,v 2.2 2021/10/26$
!
! Author: ANKIT GOEL, University of Michigan, 10/2021
!
! Comments: Following Routines to run RCAC one step algorithm are defined
!   initialize_RCPE  ::  Initializes all the variables used by RCAC one step 
!                   	 algorithm. Variables are defined in RCAC_variables.F95
!   identity    	 ::  Creates the identity matrix
!   weights     	 ::  Creates the weighting matrices
!   kronecker_AB	 ::  Calculates C = A kron B
!   RCPE 	 		 ::  Uses RCAC algorithm to find u(k)
!
!-----------------------------------------------------------------------------


subroutine initialize_RCPE(lu_val, ly_val, nf_val, lambda_val,W_Rtheta_val)
	use RCPE_variables
	implicit none

	INTEGER, intent(in) :: lu_val, ly_val, nf_val
	DOUBLE PRECISION, intent(in) :: W_Rtheta_val, lambda_val
	INTEGER :: ii = 0
	lu 	= lu_val
    ly 	= ly_val
    lz 	= ly
	nf 	= nf_val
	lambda 		= lambda_val
	W_Rtheta 	= W_Rtheta_val
	! print *, lambda, W_Rtheta
    Nc 	= 0

	lphi 	= Nc * (lu + lz) + lz 			! phi includes the integrator
	ltheta 	= lu * ( Nc * (lu + lz) + lz ) 	! phi includes the integrator
  
	
	! 16 Jan 2017: Following variables are allocated for MIMO Gf Optimization
	ALLOCATE(g_int (lz,1) )
	ALLOCATE(phi (Nc*(lu+lz)+lz,1) )
	ALLOCATE(phi_kron (lu, lu*Nc*(lu+lz)+lu*lz) )
	
	ALLOCATE(PHI_window (2*(Nc+1)*lu, ltheta))
	ALLOCATE(u_window (lu, 2*(Nc+1)))
	ALLOCATE(z_window (lz, 2*(Nc+1)))
		
	ALLOCATE(N_filt (lz, nf * lu))
	ALLOCATE(theta_star (ltheta, 1))
	
	ALLOCATE(P_RLS (ltheta, ltheta))
	ALLOCATE(I_lu (lu, lu))
	ALLOCATE(Phi_b_rr (nf*lu, ltheta))
	ALLOCATE(U_b_rr (nf*lu, 1))
	ALLOCATE(PHI_filt (lz, ltheta))
	ALLOCATE(PtimesPHI_filt (ltheta, lz))
	ALLOCATE(Gamma_P (lz, lz))
	
	g_int 		= 0
	PHI_window 	= 0
	u_window 	= 0
	z_window 	= 0 

	Phi_b_rr = 0
	U_b_rr   = 0
	
	N_filt 	= 0
	DO ii = 1,lz
		N_filt(ii,nf*lu) = -1
	END DO
	!N_filt = reshape( (/-.1, 0.0 ,0.0,-.1/), shape(N_filt))
	theta_star 	= 0
		
	! Initialize the covariance matrix
	CALL identity(P_RLS, ltheta,1.0/W_Rtheta)
	CALL identity(I_lu, lu, 1D0)

end subroutine initialize_RCPE

subroutine RCPE(u_out, theta_out, kk, u_in, z_in, y_in)
	USE RCPE_variables
  
	implicit none
	integer rr
	integer, intent(in) :: kk
	DOUBLE PRECISION, intent(in), DIMENSION(lu,1) :: u_in
	DOUBLE PRECISION, intent(in), DIMENSION(lz,1) :: z_in
	DOUBLE PRECISION, intent(in), DIMENSION(ly,1) :: y_in
	DOUBLE PRECISION, intent(out), DIMENSION(lu,1) :: u_out
	DOUBLE PRECISION, intent(out), DIMENSION(ltheta,1) :: theta_out
	
	
  
	
  
	! --------------------------------------------------------------
	! Regressor
	! --------------------------------------------------------------
	! Put data in control and performance buffers	
	u_window 		= cshift(u_window,-1,2)
	z_window 		= cshift(z_window,-1,2)	
	PHI_window 		= cshift(PHI_window,-lu,1)

	u_window(1:lu,1:1) = u_in !reshape(u_in, shape(u_window(1:lu,1)))
	z_window(1:lz,1:1) = z_in !reshape(z_in, shape(z_window(1:lz,1)))
  
	! Construct regressor, phi. Size Nc*(lu+lz)+lz
	g_int 						= g_int + z_in
	!phi(1:Nc*lu,1) 			= reshape( 1*u_window(:,1:Nc), shape(phi(1:Nc*lu,1)) )
	!phi(Nc*lu+1:Nc*lu+Nc*lz,1)	= reshape( 0*1+1*z_window(:,1:Nc), shape(phi(Nc*lu+1:Nc*lu+Nc*lz,1)) )	
	phi(Nc*lu+Nc*lz+1:,1)		= reshape( 1*g_int, shape(phi(Nc*lu+Nc*lz+1:,1)) )	
	
	! Construct regressor, PHI. Size lu by lu*(Nc*(lu+lz)+lz)
	CALL KRONECKER_AB(phi_kron, I_lu, lu, lu, transpose(phi),1, Nc*(lu+lz)+lz)
  
	PHI_window(1:lu,:) 	= phi_kron !reshape(z_in, shape(z_window(1:lz,1)))

	! --------------------------------------------------------------
	! Filter
	! --------------------------------------------------------------
	
	rr = 1
	Phi_b_rr = PHI_window(lu*rr+1:lu*rr+lu*nf,:)
	call vectorize(u_window(:,rr:rr+nf-1 ), lu, nf, U_b_rr)

	PHI_filt = matmul(N_filt, Phi_b_rr)
	U_filt = matmul(N_filt, U_b_rr)
    
    
	if ( kk> max(Nf+1, Nc+1) )  then
	       
		PtimesPhi_filt = matmul(P_RLS, transpose(PHI_filt))
		Gamma_P = 1/(lambda + matmul(PHI_filt, PtimesPhi_filt))
		P_RLS = P_RLS - matmul(matmul(PtimesPhi_filt,Gamma_P), transpose(PtimesPhi_filt))
		P_RLS = P_RLS/lambda

		theta_star 	= theta_star - matmul(PtimesPhi_filt, matmul(PHI_filt, theta_star)+z_in-U_filt)

		theta_out 	= theta_star
		u_out = matmul(phi_kron, theta_out)
	else 
		theta_out 	= 0
		u_out 		= u_in
	end if
  
	! call AppendMatrix2File(transpose(u_in), 1, lu,     'u_window11111111.dat')
	! call AppendMatrix2File(transpose(z_in), 1, lz,     'z_window11111111.dat')
	call AppendMatrix2File(transpose(u_out), 1, lu, 			'u_RCMR1111111111.dat')
	call AppendMatrix2File(transpose(theta_out), 1, ltheta,     'theta11111111111.dat')
	call writeMatrix2File(phi_kron,lu, lu*Nc*(lu+lz)+lu*lz, 	'phi_kron11111111.dat')
	call writeMatrix2File(P_RLS,ltheta, ltheta,       			'P_RLS11111111111.dat')
	! write(*,*) "RCMR output is ", kk, u_out, u_in
	!write(*,*) "Optimized estimator coefficients are ", theta_out
	!write(*,*) "Optimized filter is ", filter_out
	! write(*,*) kk, u_out, PHI_filt, U_b_rr,theta_star, lambda
end subroutine RCPE




subroutine identity(R, L, weight)
	use RCPE_variables
	implicit none    
	integer :: ii
	integer, intent(in) :: L
	DOUBLE PRECISION , intent(in)   :: weight
	DOUBLE PRECISION, DIMENSION(L,L) :: R
	R = 0
	DO ii = 1,L
		R(ii,ii) = weight
	END DO
end subroutine identity



! 19 Jan 2017: Following subroutine computes the kronecker product of two matrices 
SUBROUTINE KRONECKER_AB(C, A, Am, An, B, Bm, Bn)
	IMPLICIT NONE

	integer, intent(in) :: Am, An, Bm, Bn
	DOUBLE PRECISION, intent(in), DIMENSION(Am, An) :: A 
	DOUBLE PRECISION, intent(in), DIMENSION(Bm, Bn) :: B 
	DOUBLE PRECISION, intent(out), DIMENSION(Am*Bm,An*Bn) :: C
	integer :: ii, jj 

	DO ii = 1,Am
		DO jj = 1,An
			C(Bm*(ii-1)+1:Bm*ii, Bn*(jj-1)+1:Bn*jj) = A(ii, jj) * B
		END DO
	END DO
END SUBROUTINE KRONECKER_AB




SUBROUTINE vectorize(U, n, m, vecU)
	implicit none

	integer :: ii
	integer, intent(in) :: n, m
	double precision, intent (in), dimension(n,m) :: U
	double precision, intent (out), dimension(n*m,1) :: vecU
	Do ii = 0,m-1
		vecU(ii*n+1:ii*n+n,1) = U(:,ii+1)
	end do
end subroutine vectorize


subroutine writeMatrix(A,n)
	implicit none

	integer :: ii
	integer, intent(in) :: n
	double precision, intent (in), dimension(n,n) :: A
	do ii = 1,n
		write(*,*), A(ii,:)
	end do


end subroutine writeMatrix

subroutine writeMatrix2File(A,n,m, filename)
	implicit none

	integer :: ii
	integer, intent(in) :: n, m
	CHARACTER(20), intent(in)  :: filename
	double precision, intent (in), dimension(n,m) :: A

	!write(*,*), 'writing ' , filename
	open(unit = 8188, file = filename, action = 'write')
	do ii = 1,n
		write(8188,*), A(ii,:)
	end do

	close(8188)


end subroutine writeMatrix2File


subroutine AppendMatrix2File(A,n,m, filename)
	implicit none

	integer :: ii
	integer, intent(in) :: n, m
	CHARACTER(20), intent(in)  :: filename
	double precision, intent (in), dimension(n,m) :: A

	logical :: exist

	inquire(file=filename, exist=exist)
	if (exist) then
		open(unit = 8188, file = filename, status="old", position="append", action = 'write')
	else
		open(unit = 8188, file = filename, status="new", action = 'write')
	end if

	!write(*,*), 'writing ' , filename
	do ii = 1,n
		write(8188,*), A(ii,:)
	end do

	close(8188)


end subroutine AppendMatrix2File