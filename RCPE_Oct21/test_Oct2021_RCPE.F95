!----------------------------------------------------------------------------
! $Id: RCAC_routines.f95,v 2.1 2018/05/14$
!                       ,v 2.2 2021/10/26$
!
! Author: ANKIT GOEL, University of Michigan, 10/2021
!
! Comments: Contains program runRCAC to test RCAC algorithm
!-----------------------------------------------------------------------------


PROGRAM runRCPE
    USE RCPE_variables
    IMPLICIT NONE
    
    
    integer :: kk = 0 
    integer :: MaxIter = 100
    DOUBLE PRECISION :: pi = 3.14159265359
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE  :: A, Ahat
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE  :: B
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE  :: C
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE  :: X, Xhat 
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE  :: Y, Yhat, Z 
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE  :: mu
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE  :: w
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE  :: theta_k
    
    !call initialize()     !! Initialize RCAC variables

    call initialize_RCPE(1, 1, 1, 1D0,1000000D0)

    ! Nf_star 	= RESHAPE( (/1/), shape(Nf_star) )
    
    ALLOCATE(A(2,2))
    ALLOCATE(Ahat(2,2))
    ALLOCATE(B(2,lu))
    ALLOCATE(C(ly,2))
    ALLOCATE(X(2,MaxIter))
    ALLOCATE(Y(ly,MaxIter))
    ALLOCATE(Xhat(2,MaxIter))
    ALLOCATE(Yhat(ly,MaxIter))
    ALLOCATE(Z(lz,MaxIter))
    ALLOCATE(mu(lu,MaxIter))
    ALLOCATE(w(1,1))
    ALLOCATE(theta_k (ltheta, 1))
    
    
    !! Define system here
    A   = reshape((/0.23, 0.1, 0.2 , 0.6/), shape(A))
    !A 	= reshape((/0.3, 0.0, 0.0 , 0.0/), shape(A))
    B 	= reshape((/0.9 ,0.3/), shape(B))
    C 	= reshape((/1.1 ,0.5/), shape(C))
    Ahat      = A
    Ahat(1,1) = 0
    
    X 	= 0
    Y 	= 0
    mu 	= 0
    Xhat = X
    Yhat = Y

    X(:,1) = 10
    Y(:,1) = matmul(C,X(:,1))
    
    DO kk = 2,MaxIter
      w           = 10+sin(2*pi/100*kk)

      X(:,kk)     = matmul(A,X(:,kk-1)) + matmul(B,w(:,1))
      Y(:,kk)     = matmul(C,X(:,kk))

      Xhat(:,kk)     = matmul(Ahat,Xhat(:,kk-1)) + matmul(B,w(:,1))
      Yhat(:,kk)     = matmul(C,Xhat(:,kk))
      
      Z(1:lz,kk)     = Y(1:ly,kk)-Yhat(1:ly,kk)
      
      call RCPE(mu(1:lu, kk), theta_k, kk, Ahat(1,1), Z(1:lz, kk), 0)
      Ahat(1,1)  = abs(mu(1, kk))
      print *, Ahat(1,1)
    end do  
    
    print *, 'All done!'
    
  END PROGRAM runRCPE
  