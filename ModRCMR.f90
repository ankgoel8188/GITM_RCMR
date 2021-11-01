!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!------------------------------------------------------------------------------
! $Id: ModRCMR.f90,v 1.3 2013/10/24 18:40:31 agburr Exp $
!
! Author: Asad
!
! Comments: Common variables needed across GITM to perform RCMR data
!           assimilation, as well as subroutines to initialize the variables
!
! AGB 3/31/13: Removed unnecessary variables, changed names to be more
!              descriptive, added variables to allow flagging from UAM.in,
!              added initialization routines to allow flagging from UAM.in
! AGB 10/23/13: Adapted to allow driving of photoelectron heating efficiency
!------------------------------------------------------------------------------

module ModRCMR

  use ModGITM, only:iProc, nProcs, Sat_Loc
  use ModSatellites, only: SatAltDat
  use ModInputs, only: iCharLen_
  use ModTime	
  implicit none

  logical :: RCMRFlag = .false.
  integer :: rcmrStepToTurnOff = 1000000000000

  !! Variables declared by Asad's old 2 step code
  integer :: row, col, max_rows, max_cols, print_i, print_j, N, M, mm, dbuffer
  integer :: C_on, ustep, f_o_count, Pc, Pcc, lu2, l_dim, AllocateStatus, idty
  integer :: DeAllocateStatus, control_on, TRUTH_or_ID, lat_p, lon_p, alt_p
  integer :: s_dhat, ii, kkkk

  double precision :: eta, reg_val

  character (len=50) :: filename
  character (len=iCharLen_) :: RCMRInType, RCMROutType
	
  integer, dimension(1,1) :: dhat

  double precision, dimension(1,1) :: usum, UB, lambda, inp, y_k

  double precision, dimension(:), allocatable :: gathered, scattered
  double precision, dimension(:), allocatable :: gathered_sza, Sat_Proc
  double precision, dimension(:,:), allocatable :: P1, R2, T, theta1
  double precision, dimension(:,:), allocatable :: y_out, w, u, up, y0s, y, y0
  double precision, dimension(:,:), allocatable :: z, zp, zav
  double precision, dimension(:,:), allocatable :: diagn, y_mat, u_mat, z_mat


 
  !ANKIT: Variables below this line are for EDC RCMR
  double precision, dimension(:,:), allocatable :: TEC_true, TEC_lon, TEC_lat
  double precision, dimension(:,:), allocatable :: TEC_currentTime
  integer, dimension(:,:), allocatable :: TEC_step


  integer :: TEC_read_IOStatus = 0
  integer :: kkk = 1                   !counter for reading TEC data
  integer :: TimeArrayDummy = 0        !dummy variable to store timearray entries
  

  !! Following variables were defined by Asad, but are still used
  double precision :: Dts, Measure_Dts, scatter
  double precision, dimension(:,:), allocatable  :: u_out
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! RCAC One step variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: lv, RegZ, lz, lu, ly
  INTEGER :: ltheta, lphi
  INTEGER :: Nc
  integer :: dummy_int  ! to store tuning setting
  DOUBLE precision, dimension(20) :: dummy_real

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Rz, Ru, Ruf, Rtheta
  DOUBLE PRECISION :: W_Rz, W_Ru, W_Ruf, W_Rtheta

  ! Filter settings
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Nu, Du, Nz, Dz, Nphi, Dphi
  INTEGER :: nf_Nu, nf_Du, nf_Nz, nf_Dz, nf_Nphi, nf_Dphi
  INTEGER :: nf_z, nf_u
  ! Xbar, Xfbar are Buffers for filters
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Ubar, Ufbar, Zbar, Zfbar
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PHIbar, PHIfbar
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u_filtered, z_filtered, phi_filtered
  
  ! u_h, z_h and y_h are needed to construct phi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u_h, y_h, z_h
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PHI
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_reg
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_reg_row
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Vphi, Uphi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: VphiVec, UphiVec
  
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: theta_h
  
  DOUBLE PRECISION :: lambda1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PP, Tau, Ginv
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: XX, Zp1, Rp, Rbarinv
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Iltheta
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PpXxpTauInv, XxTheta, Rzp, PRdelta, dTheta
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: thetaout, uout         ! Theta(k) and u(k)
  
  ! Lapack inversion
  INTEGER :: INFO = 0 
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV, WORK
  double precision :: alpha                !ANKIT: multiplies z
  integer :: uFiltLength = 90              !Ankit:25Jan15 rcmr filter length

  ! Variables to read and write TEC data - ANKIT 1 January 2015
  integer :: points = 0
  real, allocatable , dimension(:,:):: TEC_location
  character (len=16), dimension(12) :: TEC_file
  integer :: TEC_proc = 0                  !Ankit19Feb2015: proc which outputs TEC data
  character(LEN = 15) :: filename1         !Ankit24Feb2015: name of tec output file
  character(LEN = 20) :: filenameRCMR1     !Ankit24Feb2015: name of tec output file
  !integer :: lz_tec = 1                    !ANkit25Feb2015: number of tec stations to be used in RCMR
  INTEGER, DIMENSION(:), ALLOCATABLE :: index_lz_tec !Ankit25Feb2015: Index of tec stations to be used in RCMR
  double precision :: Gf_gain = 1

  ! 15 Jan 2017: Following variables are defined for MIMO Gf Optimization
  ! 	Variables lz, lu, ly, ltheta, lphi, Nc are already defined with the 
  !	previous implementations of RCAC.
  INTEGER :: nf 					! Order of the FIR filter to be optimized
  INTEGER :: pc_GO, pn 				! Window sized
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: g_int 						! Integrator. Size lz by 1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_GO, phi_kron_GO 		! Regressor phi. Size Nc(lu+lz)+lz by 1. Size lu by lu (Nc(lu+lz)+lz).
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PHI_window					! PHI buffer. pn*lu by ltheta
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u_window 					! Control buffer. lu by pn 
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: z_window 					! Performance buffer. lz by pc
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Z_GO 							! lz by pc
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Zc, Z_hat_c 				! lz*pc by 1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PHI_b, U_b 				! nf*lu by ltheta. nf*lu by 1.
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Uhat_bar_unf, U_bar		! Used for filter optimization. nf*lu by pc*ltheta. nf*lu by pc
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PHI_bar_unf, U_bar_unf  	! Used for theta optimization. nf*lu*pc by ltheta. nf*lu*pc by 1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Nf_star, theta_star 		! OPtimized filter and controller gains. lz by nf*lu. ltheta by 1.
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Nf_star_pc, theta_star_pc 	! OPtimized filter and controller gains. pc*lz by pc*nf*lu. pc*ltheta by pc.
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PHI_bar_f, U_bar_f			! lz*pc by ltheta. lz*pc by 1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Uhat_bar, V_k 				! nf*lu by pc. Both.
  
	
  
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: I_pc						! pc by pc

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u_RCMR, u_avg

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RCMR_est, RCMR_est0 !ANK 12/04/17 Multiple estimates updated by RCMR
  
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_RLS

  real, DIMENSION(:), ALLOCATABLE :: Nf_Array

  real :: MinRho_R, MaxRho_R, MeanRho_R

  double precision, dimension(:), allocatable :: MinRho_true, MaxRho_true, MeanRho_true, Rho_step
  double precision :: real_dummy = 0.0
  
  !logical :: checkDifferencesFlag = .false.
  
contains 
  
  subroutine alloc_rcmr
    
    ! --------------------------------------------------------------
    ! Set the dimensional RCMR parameters
    ! --------------------------------------------------------------
    
    max_cols    = 1000000   !20 days * 24 hours * 60 min * 30 two seconds = 864000
    max_cols    = 3500000   !20 days * 24 hours * 60 min * 30 two seconds = 864000
    max_rows    = max_cols
    l_dim       = 600
    ly          = lz
    
    W_Rz  = 1.0
    W_Ru  = 0.0
    W_Ruf = 0.0
    
    IF (W_Ru == 0.0) THEN
       lv = lz
    ELSE
       lv = lz + lu
    ENDIF

    !! Filter allocation and initialization
    !! Following implementaion is for fixed Gf
    nf_Nz   = 3
    nf_Dz   = 3
    nf_Nu   = 3 
    nf_Du   = 3
    nf_Nphi = 3
    nf_Dphi = 3
    nf_z = 3
    nf_u = 3
    
    !Ankit 14 March 2015: Expanding filter size to implement Markov Parameters
    nf_z = 20
    nf_u = 20
    
    dummy_real = 0
    
    !! size calculations
    IF (RegZ == 1) THEN
       ltheta  = Nc * lu * (lu + lz)
       lphi    = Nc * (lu + lz)
    ELSE
       ltheta  = Nc * lu * (lu + ly)
       lphi    = Nc * (lu + ly)
    ENDIF
    
    !! Regressor buffer initialization
    Allocate(u_out      (lu, 1) )
    ALLOCATE(u_h        (lu, NC))
    ALLOCATE(z_h        (lz, NC+1))
    ALLOCATE(y_h        (ly, NC+1))
    ALLOCATE(theta_h    (ltheta, 2))
    ALLOCATE(phi_reg    (lphi, 1))
    ALLOCATE(Uphi       (lu, NC))
    ALLOCATE(UphiVec    (lu* NC, 1))
    ALLOCATE(PHI        (lu, lu*lphi))
    IF (RegZ == 1) THEN
       ALLOCATE(Vphi       (lz, NC))
       ALLOCATE(VphiVec    (lz* NC, 1))
    ELSE
       ALLOCATE(Vphi       (ly, NC))
       ALLOCATE(VphiVec    (ly* NC, 1))
    ENDIF
    u_h     = 0
    z_h     = 0
    y_h     = 0
    theta_h = 0
    phi_reg = 0
    Uphi    = 0
    Vphi    = 0
    UphiVec = 0
    VphiVec = 0
    
    !! Buffer values
    ALLOCATE(Ubar        (lu, nf_u))
    ALLOCATE(Ufbar       (lz, nf_u))
    
    ALLOCATE(Zbar        (lz, nf_z+1))
    ALLOCATE(Zfbar       (lz, nf_z))
    
    ALLOCATE(PHIbar      (lu*(nf_u+1), ltheta))
    ALLOCATE(PHIfbar     (lz*nf_u, ltheta))
    
    ALLOCATE(u_filtered  (lz, 1))
    ALLOCATE(z_filtered  (lz, 1))
    ALLOCATE(phi_filtered(lz, ltheta))
    
    Ubar    = 0
    Ufbar   = 0
    Zbar    = 0
    Zfbar   = 0
    PHIbar  = 0
    PHIfbar = 0
    u_filtered = 0
    z_filtered = 0
    phi_filtered = 0
    
    !! Initialize Cost function Weights
    ALLOCATE(Rz        (lz, lz))
    ALLOCATE(Ruf       (lz, lz))
    ALLOCATE(Ru        (lu, lu))
    ALLOCATE(Rtheta    (ltheta, ltheta))
    
    call identity(Rz, lz, W_Rz)
    call identity(Ru, lu, W_Ru)
    call identity(Ruf, lz, W_Ruf)
    call identity(Rtheta, ltheta, W_Rtheta)
    
    ALLOCATE(XX     (lv, ltheta))
    ALLOCATE(zp1    (lv, 1))
    ALLOCATE(Rp     (lv, lv))
    ALLOCATE(Rbarinv(lv, lv))
    ALLOCATE(Tau    (lv, lv))
    ALLOCATE(Ginv   (lv, lv))
    ALLOCATE(PpXxpTauInv (ltheta, lv))
    ALLOCATE(XxTheta (lv, 1))
    ALLOCATE(Rzp     (lv, 1))
    ALLOCATE(IPIV    (lv))
    ALLOCATE(WORK    (lv))
    ALLOCATE(PRdelta (ltheta, ltheta))
    ALLOCATE(dTheta  (ltheta, 1))
    ALLOCATE(Iltheta    (ltheta, ltheta))
    call identity(Iltheta, ltheta, 1.0_8)
    XX = 0     
    zp1 = 0
    Rp = 0
    Tau = 1
    Ginv = 0
    PpXxpTauInv = 0
    XxTheta = 0
    Rzp = 0
    
    !! Covariance Initialization
    ALLOCATE(PP     (ltheta, ltheta))
    call identity(PP, ltheta, 1/W_Rtheta)
    
    !! theta(k) and u(k) initialization
    ALLOCATE(thetaout  (ltheta, 1))
    ALLOCATE(uout      (lu,1))
    thetaout    = 0
    uout        = 0  
    
    ! --------------------------------------------------------------
    ! Allocate space for the RCMR variables
    ! --------------------------------------------------------------
    
    allocate(gathered(nProcs), STAT = AllocateStatus)
    allocate(gathered_sza(nProcs), STAT = AllocateStatus)
    allocate(scattered(nProcs), STAT = AllocateStatus)
    allocate(Sat_Proc(nProcs), STAT = AllocateStatus)
    allocate(T(lz,lu), STAT = AllocateStatus)
    allocate(y0s(nProcs,max_cols), STAT = AllocateStatus)
    allocate(w(1,max_rows), STAT = AllocateStatus)
    allocate(u(lu,max_rows), STAT = AllocateStatus)
    allocate(y_out(max_rows,1), STAT = AllocateStatus)
    allocate(y(lu,max_rows), STAT = AllocateStatus)
    allocate(y0(lz,max_rows), STAT = AllocateStatus)
    allocate(z(lz,max_rows), STAT = AllocateStatus)
    allocate(zp(lz,max_rows), STAT = AllocateStatus)
    allocate(up(lz,max_rows), STAT = AllocateStatus)
    allocate(zav(lz,max_rows), STAT = AllocateStatus)
    allocate(diagn(1,max_rows), STAT = AllocateStatus)
    allocate(z_mat(lz,l_dim), STAT = AllocateStatus)
    allocate(u_mat(lu,l_dim), STAT = AllocateStatus)
    allocate(y_mat(ly,l_dim), STAT = AllocateStatus)
    allocate(R2(lz,lz), STAT = AllocateStatus)
    allocate(P1(lu*Nc+ly*(Nc),lu*Nc+ly*(Nc)), STAT = AllocateStatus)
    allocate(theta1(1,lu*Nc+ly*(Nc)), STAT = AllocateStatus)
    
    !ANKIT: Variables below this line are for EDC RCMR
    !allocate(TEC_true(max_cols,lz))
    !allocate(TEC_step(max_cols,lz))
    !allocate(TEC_Lon(max_cols,lz))
    !allocate(TEC_Lat(max_cols,lz))
    !allocate(TEC_currentTime(max_cols,lz))
    
    !Ankit25Feb2015
    allocate(TEC_true(max_cols,points))
    allocate(TEC_step(max_cols,points))
    allocate(TEC_Lon(max_cols,points))
    allocate(TEC_Lat(max_cols,points))
    allocate(TEC_currentTime(max_cols,points))
    
    
    allocate(MinRho_true(max_cols))
    allocate(MaxRho_true(max_cols))
    allocate(MeanRho_true(max_cols))
    allocate(Rho_step(max_cols))
    
    !read(222, * ,IOSTAT=TEC_read_IOStatus) Dts
    !read(222, * ,IOSTAT=TEC_read_IOStatus) Measure_Dts
    !read(222, * ,IOSTAT=TEC_read_IOStatus) C_on
    !read(222, * ,IOSTAT=TEC_read_IOStatus) alpha
    !! write(*,*) "RCMR tuning settings are ", Nc, lambda1, W_Rtheta, Nu, Dts, Measure_Dts, C_on, alpha
    !read(222, * ,IOSTAT=TEC_read_IOStatus) uFiltLength
    !close(222)
    !write(*,*) "Read RCMR stuff from RCMR_tuning.in"
    
    
    ! 16 Jan 2017: Following variables are allocated for MIMO Gf Optimization
    ! nf                      = 1
    ! Nc                      = 10   
    pc_GO                   = 26000
    pn                      = pc_GO + nf + Nc - 1
    ltheta                  = lu * ( Nc * (lu + lz) + lz )
                 
    
    ALLOCATE(g_int (lz,1) )
    ALLOCATE(phi_GO (Nc*(lu+lz)+lz,1) )
    ALLOCATE(phi_kron_GO (lu, lu*Nc*(lu+lz)+lu*lz) )
    
    ALLOCATE(PHI_window (pn*lu, ltheta))
    ALLOCATE(u_window (lu, pn))
    ALLOCATE(z_window (lz, pc_GO))
    
    ALLOCATE(Z_GO (lz, pc_GO) )
    ALLOCATE(Zc (lz*pc_GO,1) )
    ALLOCATE(Z_hat_c (lz*pc_GO,1) )
    
    ALLOCATE(PHI_b (nf * lu, ltheta))
    ALLOCATE(U_b (nf * lu, 1))
    
    !ALLOCATE(Uhat_bar_unf (nf * lu, pc_GO * ltheta))
    !ALLOCATE(U_bar (nf * lu, pc_GO ))
    
    !ALLOCATE(PHI_bar_unf (pc_GO * nf * lu, ltheta))
    !ALLOCATE(U_bar_unf (pc_GO * nf * lu, 1))
    
    ALLOCATE(Nf_star (lz, nf * lu))
	ALLOCATE(theta_star (ltheta, 1))
    
    !ALLOCATE(Nf_star_pc (pc_GO*lz, pc_GO* nf * lu))
    !ALLOCATE(theta_star_pc (pc_GO*ltheta, pc_GO))
    
    !ALLOCATE(PHI_bar_f (lz*pc_GO, ltheta) )
    !ALLOCATE(U_bar_f (lz*pc_GO,1) ) 
    !ALLOCATE(Uhat_bar (nf*lu, pc_GO) )
    !ALLOCATE(V_k  (nf*lu, pc_GO) )  
    

    
    !ALLOCATE(I_pc (pc_GO,pc_GO) )
    !call identity(I_pc, pc_GO, 1D0) !Move this to initialization

    ALLOCATE(u_RCMR (lu, pc_GO) )
    ALLOCATE(u_avg  (lu, pc_GO) )

    ALLOCATE(P_RLS (ltheta,ltheta) )
    
    u_RCMR = 0.0
    u_avg  = 0.0
    
    
    PHI_window              = 0.0
    u_window                = 0.0
    z_window                = 0.0
    
    PHI_b                   = 0.0
    U_b                     = 0.0
    
    !Uhat_bar_unf            = 0.0
    !U_bar                   = 0.0
    
    !PHI_bar_unf             = 0.0
    !U_bar_unf               = 0.0
    
    Nf_star                 = 0.0
    theta_star              = 0.0
    
    !Nf_star_pc              = 0.0
    !theta_star_pc           = 0.0
    
    !PHI_bar_f               = 0.0
    !U_bar_f                 = 0.0
    !Uhat_bar                = 0.0
    !V_k                     = 0.0
    
    

    !Ankit18Sept2017: Filter initialized only once for ACS
    Nf_star = 0
    DO ii = 1,lz
       !Nf_star(iii,nf*lu-1+iii) = -1
       Nf_star(ii,nf*lu) = -1
    END DO
    
  end subroutine alloc_rcmr


  ! sub routine to generate EYE(size)
  subroutine identity(R, L, weight)
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
  
  

  ! 15 Jan 2017: Following subroutine computes the pinv of matrix A.
  SUBROUTINE PSUEDOINVERSE(Apinv, A, ltheta)
    use ModBlasLapack
    !use ModBlasLapack, only: LAPACK_getrf, LAPACK_getrs
    IMPLICIT NONE
    
    integer, intent(in) :: ltheta
    DOUBLE PRECISION, intent(in), DIMENSION(ltheta,ltheta) :: A 
    DOUBLE PRECISION, intent(out), DIMENSION(ltheta,ltheta) :: Apinv
    
    integer :: LWORK, INFO, ii
    DOUBLE PRECISION, DIMENSION(ltheta,ltheta) :: Sinv, U, VT
    DOUBLE PRECISION, DIMENSION(ltheta) :: S
    DOUBLE PRECISION, DIMENSION(5*ltheta+1) :: WORK
    
    LWORK = 5*ltheta+1
    
    U 		= 0
    VT 		= 0
    Sinv 	= 0
    
    !call DGESVD( 'A', 'A', ltheta, ltheta, A, ltheta, S, U, ltheta, VT, ltheta, WORK, LWORK, INFO )
    
    Do ii = 1,ltheta
       if (abs(S(1)/S(ii)) < 1D15) then
          Sinv(ii,ii) = 1/S(ii)
       endif
       !print *, S(1)/S(ii) , (S(1)/S(ii) < 1D15)
    end do
    Apinv = matmul(transpose(VT), matmul(Sinv,transpose(U)))
    
    if (0>1) then
       !write(*,*) 'Apinv is'
       Do ii = 1,ltheta
          print *, Apinv(ii,1:ltheta)
       end do
    endif
  END SUBROUTINE PSUEDOINVERSE
  
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
  
  
  subroutine init_markov_matrix
    
    ! T is a Markov parameter matrix that must be tuned for each type of
    ! assimilation
    
    if(RCMRInType == "RHO") then
       ! Markov matrix for terrestrial neutral mass density
       ! Assumes only one data input source
       
       if(RCMROutType == "F107") then
          T = reshape((/ 0.15 /), shape(T))
       else if(RCMROutType == "PHOTOELECTRON") then
          ! AGB 7/17/13: This is not settled yet
          write (*,*) "AGB RCMR WARNING: this is a test matrix"
          T = reshape((/ 0.15 /), shape(T))
       elseif(RCMROutType == "COND") then
         if (iProc == 0) then 
           write (*,*) "RCMRing Thermal conductivity"
         end if
       else
          write (*,*) "No Markov matrix for this output type: ", RCMROutType
          RCMRFlag = .false.
       end if
       
    elseif(RCMRInType == "TEC") then
       if(RCMROutType == "EDC") then
          write (*,*) "AGB RCMR WARNING: this is a test matrix"
          !T = reshape((/ 0.15 /), shape(T))  !Ankit 06July17: Removed since T is lz by lu. I dont want to keep initializing this. 
       endif
    else
       ! Markov matrix has not been established
       write (*,*) "No Markov matrix for this output type: ", RCMROutType
       RCMRFlag = .false.
    end if
    
    if (iProc == 0) then
      write (*,*) "Exiting init_markov_matrix"
    end if
  end subroutine init_markov_matrix
  
subroutine init_rcmr
  if (iProc == 0) then
    print *, "Now in init_rcmr"
  end if
  !SatAltDat = -1e32    ! AGB: changed from 1
  Sat_Loc   = 1
  
  ! --------------------------------------------------------------
  ! Set these RCMR parameters
  ! --------------------------------------------------------------

  TRUTH_or_ID = 1
  !  Dts         = 2.0  !ANKIT: Set from RCMR_tuning.in
  !  Measure_Dts = 60.0

  col         = 1
  eta         = 0.0
  reg_val     = 100.0
  !  C_on Ensures that estimates are not made before the model has settled
  !  C_on        = 120!1460
  dbuffer     = C_on-20 !100!1440
  lambda(1,1) = 0.9999
  s_dhat      = 1
  dhat(1,1)   = 1
  ustep       = 1



  do mm = 0,nProcs	
     if (iProc == mm) then
        lat_p = 9
        lon_p = 1
        alt_p = 36
     end if
  end do

  u_out      = 0.0
  Sat_Proc   = 0.0
  y_mat      = 0.0
  z_mat      = 0.0
  u_mat      = 0.0
  theta1     = 0.0
  control_on = 0

  


  u(:,:)       = 0.0
  gathered(:)  = 0.0
  scattered(:) = 0.0
  usum         = 0.0
  P1           = 0.0

  do idty = 1, Nc*(ly+lu)
     P1(idty,idty) = reg_val
  end do

  R2 = 0.0
  do idty = 1, lz
     R2(idty,idty) = 1.0
  end do

  !ANKIT: Variables below this line are for EDC RCMR
  TEC_true = 0
  TEC_step = 0
  TEC_lon  = 0
  TEC_lat  = 0
  TEC_currentTime = 0

  MinRho_true = 0
  MaxRho_true = 0
  MeanRho_true = 0

  
  !write(*,*) "Now reading TEC data "

  !ANKIT: This piece opens the following files, and 
  !       reads the data in to the arrays
  if ((iproc == 0 ) .and. ((RCMROutType == "EDC") .OR. (RCMROutType == "ManyEDC") ) ) then
     TEC_file(1) = "Data_TRUE_1.dat"
     TEC_file(2) = "Data_TRUE_2.dat"
     TEC_file(3) = "Data_TRUE_3.dat"
     TEC_file(4) = "Data_TRUE_4.dat"
     TEC_file(5) = "Data_TRUE_5.dat"
     TEC_file(6) = "Data_TRUE_6.dat"
     TEC_file(7) = "Data_TRUE_7.dat"
     TEC_file(8) = "Data_TRUE_8.dat"
     TEC_file(9) = "Data_TRUE_9.dat"
     TEC_file(10) = "Data_TRUE_10.dat"
     TEC_file(11) = "Data_TRUE_11.dat"
     TEC_file(12) = "Data_TRUE_12.dat"
          

     do ii = 1,points
        !write(*,*) "Reading TEC data from file ", iproc, ii, TEC_location(ii,1)*180/3.14, TEC_location(ii,2)*180/3.14
        open(unit = 22, file = TEC_file(ii), action = 'Read')
        do
           read(22, * ,IOSTAT=TEC_read_IOStatus) TEC_Step(kkk,ii), TEC_currentTime(kkk,ii), TimeArrayDummy, &
                TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, &
                TEC_lon(kkk,ii), TEC_lat(kkk,ii), TEC_true(kkk,ii)
           if (TEC_read_IOStatus < 0) exit
           kkk = kkk + 1
        end do
        do kkk=1,100
           write(*,*) TEC_step(kkk,ii),  TEC_currentTime(kkk,ii), TEC_lon(kkk,ii), TEC_lat(kkk,ii), TEC_true(kkk,ii)
        end do
        
        close(22)
        !write(*,*) "TEC data read from ", TEC_file(ii)
        kkk = 1
     end do
  elseif ((iproc == 0 ) .and. (RCMROutType == "COND") )  then
     !write(*,*) "Reading Density data from file Data_TRUE_1.dat on Processor", iproc
     open(unit = 22, file = 'Data_TRUE_1.dat', action = 'Read')
     do
        read(22, * ,IOSTAT=TEC_read_IOStatus) Rho_step(kkk), real_dummy, TimeArrayDummy, &
             TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, TimeArrayDummy, &
             MinRho_true(kkk), MaxRho_true(kkk), MeanRho_true(kkk)
        if (TEC_read_IOStatus < 0) exit
        kkk = kkk + 1
     end do
     do kkk=1,10
        write(*,*) Rho_step(kkk), (MinRho_true(kkk)), (MaxRho_true(kkk)), (MeanRho_true(kkk)), real_dummy
     end do

     close(22)
     !write(*,*) "Density data read from Data_TRUE_1.dat"
     
  endif
  

end subroutine init_rcmr


!Ankit:09May2015 This routine reads data from Data_TRUE_locations.in
!      Lat and Lon are stored in degrees, this routine converts them
!      to radians while storing them in TEC_Location(:,:)
subroutine read_tec_locations1 !Added 1 to the name. This routine is not called anywhere
	!use ModRCMR
	implicit none
	!integer, intent(in) :: RCMRFlag1
	integer :: read_stat = 0
	integer :: iii = 1

	open(unit = 111111, file = "Data_TRUE_locations.in", action = 'read')
	read(111111,*) TEC_proc
	read(111111,*) points
	allocate(TEC_location(points,2))
	do iii = 1,points !while (iii < points+1)
		read(111111, *, IOSTAT = read_stat) TEC_location(iii,1), TEC_location(iii,2)
		TEC_location(iii,1) = TEC_location(iii,1) * 3.141592653589793/180
		TEC_location(iii,2) = TEC_location(iii,2) * 3.141592653589793/180
	end do
	close(111111)
	!write(*,*) "> Read locations for TEC calculations"


	!  if (RCMRrun=='TRUTH') then
	!     iii = 1
	!     do while (iii < points+1)
	!        if (iii < 10) then
	!           write(filenameRCMR1, "(A14,I1,A4)") "Data_TRUE_RCMR_", iii, '.dat'
	!        else
	!           write(filenameRCMR1, "(A14,I2,A4)") "Data_TRUE_RCMR_", iii, '.dat'
	!        end if
	!        open(unit = 111111+iii, file = trim(filenameRCMR1), action = 'write')
	!        iii = iii +1
	!     end do
	!  elseif (RCMRrun=='ID') then
	!     do while (iii < points+1)
	!        if (iii < 10) then
	!           write(filename1, "(A9,I1,A4)") "Data_TRUE_", iii, '.dat'
	!        else
	!           write(filename1, "(A9,I2,A4)") "Data_TRUE_", iii, '.dat'
	!        end if
	!        open(unit = 111111+iii, file = trim(filename1), action = 'write')
	!        iii = iii +1
	!     end do
	!  end if

end subroutine read_tec_locations1

subroutine write_TEC_data1 !added 1 to the name. This routine doesnt work anyways
	implicit none
	logical:: exist
	real :: VTEC_interp, sza_test
	integer :: iii, ii_tec
	character (len=50) :: RCMRrun = 'nothing' ! this doesnt do anything
	!Ankit:09May2015: This piece can't be put in a subroutine, as the file
	!      pointer is lost with the subroutine. It can be done better, but 
	!      I dont know how to do it better.
	!      Do not waste time thinking why I couldnt put it in a subroutine
	!Ankit23May16: Put the TEC writing in a subroutine
	if (RCMRrun == 'ID') Then
		do iii=1,points! while (iii < points+1)
			if (iii < 10) then
				write(filenameRCMR1, "(A14,I1,A4)") "Data_TRUE_RCMR_", iii, '.dat'
			else
				write(filenameRCMR1, "(A14,I2,A4)") "Data_TRUE_RCMR_", iii, '.dat'
			end if

			inquire(file=filenameRCMR1, exist=exist)
			if (exist) then
				open(unit = 111111+iii, file = trim(filenameRCMR1), status="old", position="append", action = 'write')
			else
				open(unit = 111111+iii, file = trim(filenameRCMR1), status="new", action = 'write')
			end if

		end do
		elseif (RCMRrun == 'TRUTH') then
			do iii=1,points! while (iii < points+1)
				if (iii < 10) then
					write(filename1, "(A9,I1,A4)") "Data_TRUE_", iii, '.dat'
				else
					write(filename1, "(A9,I2,A4)") "Data_TRUE_", iii, '.dat'
				end if
				inquire(file=filename1, exist=exist)
				if (exist) then
					open(unit = 111111+iii, file = trim(filename1), status="old", position="append", action = 'write')
				else
					open(unit = 111111+iii, file = trim(filename1), status="new", action = 'write')
				end if
   
			end do
		end if


		!! Ankit 24Jan2015 - Added TEC writing at preset locations
		if ((iproc == TEC_proc) .AND. .true.) then
			do ii_tec = 1,points !while (ii_tec < points+1)
				call calc_single_vtec_interp(TEC_location(ii_tec,1), TEC_location(ii_tec,2), VTEC_interp)
				call get_sza(TEC_location(ii_tec,1), TEC_location(ii_tec,2), sza_test)
				write(111111+ii_tec , &
				"(I7, 1X, F15.4, 1X, I4, 1X, I2, 1X, I2, 1X, I2, 1X, I2, 1X,I2, 1X, I3, 1X, F9.3, 1X, F9.3, 1X, F12.8, 1X, F12.8)") &
				iStep, CurrentTime, iTimeArray(1), iTimeArray(2), iTimeArray(3), iTimeArray(4), iTimeArray(5), iTimeArray(6), &
				iTimeArray(7), TEC_location(ii_tec,1)*180/3.141592653589793, &
				TEC_location(ii_tec,2)*180/3.141592653589793 , VTEC_interp, sza_test
			end do
		end if


		!  close TEC files - ANKIT
		do while (ii_tec < points+1)
			close(111111+ii_tec)
			ii_tec = ii_tec+1
		end do

        
	end subroutine write_TEC_data1	
end module ModRCMR
