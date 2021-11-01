!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!----------------------------------------------------------------------------
! $Id: RCMR_routines.f90,v 1.5 2013/10/24 18:53:28 agburr Exp $
!
! Author: Asad, UMichigan, 2/2013
!
! Comments: Routines to run RCMR data assimilation
!
! AGB 3/31/13: Added comments, changed variable names to be more descriptive,
!              changed subroutine format for consistency with GITM files,
!              removed test output files, removed unused variables, streamlined
!              MPI calls.
! AGB 10/23/13: Adapted to allow driving of photoelectron heating efficiency
!-----------------------------------------------------------------------------
!
module arrayMod
  real, dimension(:), allocatable :: satLon, satLat, satAlt, satMeanDensity, &
                                     sat2Lon, sat2Lat, sat2Alt, sat2Density
  integer, dimension(:), allocatable :: satIndex, satYear, satMonth, satDay, &
                                        satHour, satMinute, satSecond
  character(len = 8), dimension(:), allocatable :: time
  character(len = 10), dimension(:), allocatable :: date  
end module arrayMod

subroutine run_RCMR
  use arrayMod
  use ModInputs
  use ModTime
  use ModGITM
  use ModMpi
  use ModRCMR
  use ModSatellites, only: SatCurrentDat, SatAltDat, nRCMRSat, RCMRSat
  use ModEUV, only: sza
  use EUA_ModMsis00, ONLY: gtd7
	
  implicit none

  integer :: status(MPI_STATUS_SIZE), iError, ii_loop 
  integer :: blockCheck = -1
  integer :: nLines = 0
  character(len = 21):: newFileName = "syntheticDataFile.txt"
  real :: output_est, ulimit, llimit, dummy_TEC_calculated, dummy
  double precision :: localVar
  double precision, dimension(:), allocatable :: TEC_calculated     !ANKIT: holds the value of TEC at AA
  double precision :: TEC_calculated_P1 = 0                         !ANKIT: holds the value of TEC at ND
  double precision :: u_sat_level = 0
  double precision, DIMENSION(lz,nf*lu) :: filter_GO
  double precision, DIMENSION(ltheta,1) :: theta_GO

  double precision, dimension(3) :: Rho_Calculated = 0
  integer :: returnProc = -1
 
  integer :: io = 0

  real :: diff = -1e32
  real :: diff2 = -1e32
  real :: diff3 = -1e32
  logical :: exist  
  double precision, dimension(90) :: MaxRho_90Min = 0
  
  integer :: iLat = -1
  integer :: iLon = -1
  integer :: iAlt = -1
  integer :: iBlock = -1 
  integer :: iiBlock = -1
  integer :: iiProc = -1  
  real :: SSLon, SSLat, ASLon, ASLat
  real :: rLon, rLat, rAlt
  real :: gitmDensity_SS = -1e32
  real :: msisDensity_SS = -1e32
  real :: gitmDensity_AS, msisDensity_AS

  real :: geo_lat, geo_lon, geo_alt, geo_lst
  real, dimension(1:2) :: msis_temp
  real, dimension(1:9) :: msis_dens
  real, dimension(7)  :: ap = 10.0
  logical :: foundLocation = .False.
  
  diff = -1e32
  diff2 = -1e32
  gitmDensity_SS = -1e32
  msisDensity_SS = -1e32
  gitmDensity_AS = -1e32
  msisDensity_AS = -1e32

  allocate(TEC_calculated(lz))
  !BP: Check number of lines in data file for RCMR, then allocated  variables with correct
  !    size. 
  if (nLines .eq. 0) then
    open (unit = 31, file = newFileName, status = "old")
    do
      read(31, *, iostat=io)
      if (io .ne. 0) then
        exit
      end if
      nLines = nLines + 1
    end do
  end if
  close(31)

  if (allocated(satYear) == .False.) then
    call arraySub(nLines)
  end if
  !!!!!!!!!!!!!!!!!!!!!

  call start_timing("RCMR")
  !call MPI_BARRIER(iCommGITM, iError)

  if ( (RCMROutType == "EDC") .OR. (RCMROutType == "ManyEDC") ) then
     if (iproc == TEC_proc) then
        do ii_loop = 1,lz
           call calc_single_vtec_interp(TEC_location(index_lz_tec(ii_loop),1), &
                TEC_location(index_lz_tec(ii_loop) ,2), TEC_calculated(ii_loop))  !Computed TEC at the TEC station
        end do
     endif
  elseif (RCMROutType == "COND") then
     !call calc_global_ave_rho(35,MinRho_R,MaxRho_R,MeanRho_R)
     !if (iproc == 0) then
        !Rho_Calculated(1) = MinRho_R
        !Rho_Calculated(2) = MaxRho_R

        !Simple queue to track previous 90 minute maxes
        !MaxRho_90min(2:90) = MaxRho_90min(1:89)
        !MaxRho_90min(1) = MaxRho_R
        !MaxRho_R = maxval(MaxRho_90min)
        !write(*,*) "Before printDensity", MaxRho_R


     !   Rho_Calculated(3) = MeanRho_R
        
     !end if
     !call MPI_Bcast(MaxRho_R, 1, MPI_DOUBLE_PRECISION, 0, iCommGITM, iError)
    
  else
     call MPI_BARRIER(iCommGITM, iError)
     localVar = SatAltDat(RCMRSat(1))
     call MPI_REDUCE(localVar, SatAltDat(RCMRSat(1)), 1, MPI_DOUBLE_PRECISION, &
          MPI_MAX, 0, iCommGITM, iError)
     if (iProc==0) then
        do mm=0,(nProcs-1)
           if (Sat_Proc(mm + 1) .ne. 0) then
              SatAltDat(RCMRSat(1)) = Sat_Proc(mm + 1)   !Ankit: Computed data at the satellite location
           end if
        end do
     end if
  endif

    
  !ANKIT: Added MPI_Bcast to update on all processors
  !ANKIT:5Feb15 - Added TEC bcasting from 12 locations
  !call MPI_Bcast( TEC_calculated(1:6)   , 6, MPI_DOUBLE_PRECISION,2, iCommGITM, iError)
  !call MPI_Bcast( TEC_calculated(7:12)  , 6, MPI_DOUBLE_PRECISION,3, iCommGITM, iError)
  call MPI_Bcast( TEC_calculated  , lz, MPI_DOUBLE_PRECISION, TEC_proc, iCommGITM, iError) !Ankit25Feb2015: Send TEC to all procs

  ! AGB: Initialize the output estimate
  if(RCMROutType == "F107") then
     output_est = f107_est
     llimit     = 70.0
     ulimit     = 400.0
  else if(RCMROutType == "PHOTOELECTRON") then
     output_est = PhotoElectronHeatingEfficiency_est
     llimit     = 0.02
     ulimit     = 0.2
  else if (RCMROutType == "EDC") then           !ANKIT: Eddy diffusion coefficient hard limits
     !output_est = EDC_est
     llimit     = 500
     ulimit     = 3000
  else if (RCMROutType == "ManyEDC") then           !ANKIT: Eddy diffusion coefficient hard limits
     !output_est = EDC_est
     llimit     = 500
     ulimit     = 3000
  elseif (RCMROutType == "COND") then
     llimit     = 0
     ulimit     = 1
  else
     write (*,*) "ERROR: unknown RCMR output type", RCMROutType
  end if
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !BP: Uncomment when not using MSIS for the two parameter run
  call readSyntheticData(nLines, lu)
  call printDensity(nLines, diff, diff2, blockCheck, lu)
  

  localVar = diff
  call MPI_REDUCE(localVar, diff, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, iCommGITM, iError)
  
  !if (lz == 2) then
  !  localVar = diff2
  !  call MPI_REDUCE(localVar, diff2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, iCommGITM, iError)
  !end if
  !!!!!!!!!!!!


  !-------------------------------------------------------------------
  !BP: 07/10/2020
  !Use MSIS subsolar and antisolar points to as truth data to compute
  !differences. 

  call get_subsolar(CurrentTime, VernalTime, SSLon, SSLat)
  !write(*,*) "Subsolar location:", SSLon*180.0/pi, SSLat*180.0/pi

  if (SSLon == 0) then
    SSLon = 1e-2
  endif

  if (SSLat == 0) then
    SSLat = 1e-2
  endif  

 ! call LocationIndex(SSLon, SSLat, iiBlock, iLon, iLat, rLon, rLat)   
 
 ! if (iiBlock .eq. 1) then
 !   call BlockAltIndex(400*1000.0, iiBlock, iLon, iLat, iAlt, rAlt)
 !   
    !MSIS rho
 !   geo_lon = SSLon*180.0/pi
 !   geo_lat = SSLat*180.0/pi
 !   geo_lon = mod(SSLon*180.0/pi + 360.0, 360.0)
 !
 !   if (geo_lat < -90.0) then
 !     geo_lat = -180.0-geo_lat
 !     geo_lon = mod(geo_lon+180.0,360.0)
 !   endif
 !   if (geo_lat >  90.0) then
 !     geo_lat =  180.0-geo_lat
 !     geo_lon = mod(geo_lon+180.0,360.0)
 !   endif
 !
 !   geo_alt = Altitude_GB(iLon, iLat, iAlt, iiBlock)/1000.0
 !   geo_lst = mod(utime/3600.0+geo_lon/15.0,24.0)
                                        
 !   CALL GTD7(iJulianDay,utime,geo_alt,geo_lat,geo_lon,geo_lst, &
 !             F107A,F107,AP,48,msis_dens,msis_temp)
 !   msisDensity_SS = msis_dens(6)
    !write(*,*) "utime", utime

 !   gitmDensity_SS = rLon*rLat*Rho(iLon, iLat, iAlt, iiBlock)       + &
 !                    (1-rLon)*rLat*Rho(iLon+1, iLat, iAlt, iiBlock) + &
 !                    rLon*(1-rLat)*Rho(iLon, iLat+1, iAlt, iiBlock) + &
 !                    (1-rLon)*(1-rLat)*Rho(iLon+1, iLat+1, iAlt, iiBlock)

    !BP: Uncomment this line when USING MSIS (important!)
    !diff = msisDensity_SS - gitmDensity_SS

    !if ( iTimeArray(6) == 0) then
    !  write(*,*) iTimeArray
    !  write(*,*) "Thermal Conductivity", ThermalConduction_AO2, ThermalConduction_AO, &
    !             ThermalConduction_s
    !  write(*,*) "GITM, MSIS, SS Diff:", gitmDensity_SS, msisDensity_SS, diff
    !endif
 ! endif

  !BP: Uncomment this line when USING MSIS (important!) 
  !localVar = diff 
  !call MPI_REDUCE(localVar, diff, 1, MPI_DOUBLE_PRECISION, &
  !                MPI_MAX, 0, iCommGITM, iError)
  !call MPI_Bcast(diff, 1, MPI_DOUBLE_PRECISION, 0, iCommGITM, iError)   

  !localVar = msisDensity_SS
  !call MPI_REDUCE(localVar, msisDensity_SS, 1, MPI_DOUBLE_PRECISION, &
  !                MPI_MAX, 0, iCommGITM, iError)
 
  !localVar = gitmDensity_SS
  !call MPI_REDUCE(localVar, gitmDensity_SS, 1, MPI_DOUBLE_PRECISION, &
  !                MPI_MAX, 0, iCommGITM, iError)
 
  !Find antisolar longitude
  if (SSLon .ge. pi) then 
    ASLon = SSLon - pi
  else
    ASLon = SSLon + pi
  endif
  
  !Antisolar latitude
  ASLat = -SSLat
  
  if (ASLon == 0) then
    ASLon = 1e-2
  endif

  if (ASLat == 0) then
    ASLat = 1e-2
  endif

  !call LocationIndex(ASLon, ASLat, iiBlock, iLon, iLat, rLon, rLat)
  !
  !
  !if (iiBlock .eq. 1) then
  !  call BlockAltIndex(400*1000.0, iiBlock, iLon, iLat, iAlt, rAlt)
    
    !MSIS rho                                                                                 
  !  geo_lon = ASLon*180.0/pi
  !  geo_lat = ASLat*180.0/pi


  !  geo_alt = Altitude_GB(iLon, iLat, iAlt, iiBlock)/1000.0
  !  geo_lst = mod(utime/3600.0+geo_lon/15.0,24.0)

    ! Call MSIS (results will be im mks units)                                                
  !  CALL GTD7(iJulianDay,utime,geo_alt,geo_lat,geo_lon,geo_lst, &
  !            F107A,F107,AP,48,msis_dens,msis_temp)
  !  msisDensity_AS = msis_dens(6)
 
  !  gitmDensity_AS = rLon*rLat*Rho(iLon, iLat, iAlt, iiBlock)       + &
  !                   (1-rLon)*rLat*Rho(iLon+1, iLat, iAlt, iiBlock) + &
  !                   rLon*(1-rLat)*Rho(iLon, iLat+1, iAlt, iiBlock) + &
  !                   (1-rLon)*(1-rLat)*Rho(iLon+1, iLat+1, iAlt, iiBlock)
    !BP: Uncomment this line when USING MSIS (important!) 
    !diff2 = msisDensity_AS - gitmDensity_AS              
 
    !if (iTimeArray(5) == 0 .and. iTimeArray(6) == 0) then
    !  write(*,*) "AS Diff:", gitmDensity_AS, msisDensity_AS, diff2
    !endif  
  !endif
  
  !when using MSIS

  !if (lz == 2) then
  !  localVar = diff2

  !  call MPI_REDUCE(localVar, diff2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, iCommGITM, iError)
  !end if


  ! Begin the first RCMR loop
  if (mod((istep-1)*Dts,Measure_Dts) == 0.0 .AND. TRUTH_or_ID==1 .AND. iProc == 0) then
     if ( (RCMROutType == "EDC") .OR. (RCMROutType == "ManyEDC") ) then
     !if (RCMROutType == "EDC") then                                                 !ANKIT 6 Oct 2014
        do ii_loop = 1,lz
           zp(ii_loop,ustep) = alpha*(TEC_calculated(ii_loop) - TEC_true(iStep,index_lz_tec(ii_loop))) !Ankit25Feb2015: added choice of tec station
           !write(*,*) "Step is        : ",  iStep, index_lz_tec(ii_loop)
           !write(*,*) "True TEC is    : ",  TEC_true(iStep,index_lz_tec(ii_loop))
           !write(*,*) "Computed TEC is: ",  TEC_calculated(ii_loop)
           !write(*,*) "Z is           : ",  zp(:,ustep)
        end do
     elseif  (RCMROutType == "COND") then         
        if (lz==1) then
            write(*,*) "Adding", diff, "to zp"
            zp(1, ustep) = diff*1.0e10  

        elseif (lz==2) then
          zp(1,ustep) = 1.0e10*diff 
          zp(2,ustep) = 1.0e10*diff2 

        elseif (lz == 3) then
           zp(1,ustep) = 1.0e10*(Rho_Calculated(1) - MinRho_true(iStep)) 
           zp(2,ustep) = 1.0e10*(Rho_Calculated(2) - MaxRho_true(iStep)) 
           zp(3,ustep) = 1.0e10*(Rho_Calculated(3) - MeanRho_true(iStep))
        end if
     else
        zp(:,uStep) = SatAltDat(RCMRSat(1)) - SatCurrentDat(RCMRSat(1))
     endif
     
     ! Set the averaged error
     if ( (RCMROutType == "EDC") .OR. (RCMROutType == "ManyEDC") ) then
        !if (RCMROutType == "EDC") then          !ANKIT 4 Nov 2014
        y(:,ustep) = TEC_calculated         !ANKIT5Feb2015
        if (uStep <= 120) then
           zav(:,ustep) = sum(zp(:,1:uStep))/uStep
        else
           zav(:,ustep) = sum(zp(:,uStep-120+1:uStep))/120
        end if
        z(:,uStep) = zav(:,uStep)           !Ankit 25Jan15 - Adds moving average filter to Z
        z(:,uStep) = zp(:,uStep)            !Ankit 25Jan15 - Removes the averaging from Z
     else
        if (uStep <= 120) then
           zav(:,ustep) = sum(zp(:,1:uStep))/uStep
        else
           zav(:,ustep) = sum(zp(:,uStep-120+1:uStep))/120
        end if
        y(1,uStep) = 1.0e12*zav(1,uStep)
        z(:,uStep) = 1.0e12*zav(1,uStep)
     endif

          
     if (RCMROutType == "EDC") then          ! ANKIT 6 Nov 2014
        ! CALL RCMR_onestep(u_out, ustep, EDC_est, zp(:,ustep), y(:,ustep))
        ! CALL RCMR_onestep(u_out, ustep, EDC_est, z(:,ustep), y(:,ustep), output_est)
        ! CALL RCAC_GfOptim(u_out_GO, theta_out_GO, filter_out, kk, u_in, z_in, y_in)
        CALL RCAC_GfOptim(u_out, theta_GO, filter_GO, ustep, EDC_est, zp(:,ustep-1), zp(:,ustep-1))

     elseif (RCMROutType == "ManyEDC") then          ! ANKIT 12/04/2017
        
        !CALL RCAC_GfOptim(u_out, theta_GO, filter_GO, ustep, RCMR_est, zp(:,ustep-1), zp(:,ustep-1))
        CALL RCAC_GfOptim(u_out, theta_GO, filter_GO, ustep, u_out, zp(:,ustep-1), zp(:,ustep-1))

     elseif (RCMROutType == "COND") then          ! ANKIT 03/06/2019
        if (iProc == 0) then
          write(*,*) uStep, u_out
        endif

        CALL RCAC_GfOptim(u_out, theta_GO, filter_GO, ustep, u_out, zp(:,ustep-1), zp(:,ustep-1))

        do ii_loop = 1,lu 
           u_out(ii_loop,1) = max(u_out(ii_loop,1),-1e-3) !previously -1e-3 when using abs()
           u_out(ii_loop,1) = min(u_out(ii_loop,1),1e-3)
        enddo
        if (iProc == 0) then
          write(*,*) uStep, u_out
        endif
     else
        !CALL RCMR_Function(llimit, ulimit, Nc, ustep, s_dhat, lz, lu, lu2, ly, &
        !     l_dim, Pc, Pcc, control_on, dbuffer, C_on, dhat, eta, lambda, usum, &
        !     T,  R2, y_mat, z_mat, u_mat, P1, u_out, theta1, UB)
     endif

     
     !up(:,ustep) = u_out(:,1)
     u_RCMR(:,ustep) = u_out(:,1)
     
     u_avg(:,ustep) = u_out(:,1)      !! Ankit:01/23/2019. Removed filtering
     

     if (RCMROutType == "ManyEDC") then !ANK:12/04/2018

        RCMR_est(1) = abs(u_avg(1,ustep))
        RCMR_est(2) = abs(u_avg(2,ustep))
        
        NorthPoleEddyCoef = RCMR_est0(1) + RCMR_est(1)
        SouthPoleEddyCoef = RCMR_est0(2) + RCMR_est(2)
        !write(*,*) "RCMR u_RCMR: ", u_RCMR(:,ustep)
        !write(*,*) "RCMR u_avg : ", u_avg(:,ustep)
        !write(*,*) "N/S EDC is : ", NorthPoleEddyCoef, SouthPoleEddyCoef

     elseif (RCMROutType == "COND") then !ANK:03/06/2019
        if (lu==1) then
           !RCMR_est(1) = 1000*abs(u_avg(1,ustep))
           RCMR_est(1) = abs(u_avg(1,ustep))

           !BP: Uncomment to let RCMR work again!
           ThermalConduction_AO = RCMR_est0(1) + RCMR_est(1)
           
        elseif (lu==2) then
           !Original 
           !RCMR_est(1) = abs(u_avg(1,ustep))
           
           !Two parameter A/s                                                                              
           RCMR_est(1) = abs(u_avg(1,ustep))*1000.0
           RCMR_est(2) = abs(u_avg(2,ustep))
       
           !ThermalConduction_AO = RCMR_est0(1) + RCMR_est(1)
           ThermalConduction_s = RCMR_est0(1) + RCMR_est(1)
           ThermalConduction_AO2 = RCMR_est0(2) + RCMR_est(2)
           
           !For the A and S estimation, we want to maintain a ratio between the 
           !atomic and molecular conductivity constants. This is the ratio from 
           !the Pavlov 2017 expressions for the individual species. ~0.73
           ThermalConduction_AO = ThermalConduction_AO2/0.73
        elseif (lu==3) then
           RCMR_est(1) = abs(u_avg(1,ustep))
           RCMR_est(2) = abs(u_avg(2,ustep))
           RCMR_est(3) = 1000*abs(u_avg(3,ustep))
           
           ThermalConduction_AO = RCMR_est0(1) + RCMR_est(1)
           ThermalConduction_AO2 = RCMR_est0(2) + RCMR_est(2)
           ThermalConduction_s = RCMR_est0(3) + RCMR_est(3)
        end if

        !BP: unedit eventually   
        !write(*,*) "RCMR u_RCMR: ", u_RCMR(:,ustep)
        !write(*,*) "RCMR u_avg : ", u_avg(:,ustep)
        if (iTimeArray(5) == 0 .and. iTimeArray(6) == 0) then
          write(*,*) "TCond is   : ", ThermalConduction_AO, ThermalConduction_AO2, ThermalConduction_s
        end if
                
     else
        
        if (ustep > C_on) then
           output_est = u_avg(1, ustep)
        else if (ustep <= C_on) then
           u_avg(:,ustep) = output_est
        end if
        
        if(RCMROutType == 'F107') then
           f107_est  = output_est
           f107a_est = f107_est
        else if(RCMROutType == "PHOTOELECTRON") then
           PhotoElectronHeatingEfficiency_est = output_est
        elseif(RCMROutType == "EDC") then !ANKIT 6 Oct 2014
           EDC_est = output_est
        end if
     endif
     ustep = ustep + 1
     
  end if

  ! END of the first RCMR loop
  
  
  if (RCMROutType == "ManyEDC") then !ANK:12/04/2017
     write(*,*) "Here:", RCMROutType, lu
     do ii_loop = 1,lu
        if (iProc==0) then
           scattered(:) = RCMR_est(ii_loop) + 1*RCMR_est0(ii_loop)
           !scattered = min(max(scattered,llimit),ulimit)
           !write(*,*) "Scattering : ", ii_loop, RCMR_est(ii_loop) + RCMR_est0(ii_loop)
        end if
        call MPI_SCATTER(scattered, 1, mpi_double_precision, scatter, 1, &
             mpi_double_precision, 0, iCommGITM, iError)
        if (ii_loop == 1) then
           NorthPoleEddyCoef = scatter
        elseif (ii_loop == 2) then
           SouthPoleEddyCoef = scatter
        end if
        !write(*,*) "Checking parameter passing", iProc, ii_loop, scatter, RCMR_est, NorthPoleEddyCoef, SouthPoleEddyCoef
     enddo
  elseif (RCMROutType == "COND") then !ANK:03/06/2019
     do ii_loop = 1,lu
        if (iProc==0) then
           !scattered(:) = RCMR_est(1) + 1*RCMR_est0(1)
           scattered(:) = RCMR_est(ii_loop) + 1*RCMR_est0(ii_loop)
        end if
        call MPI_SCATTER(scattered, 1, mpi_double_precision, scatter, 1, &
             mpi_double_precision, 0, iCommGITM, iError)
     
        if (ii_loop == 1) then
!          BP: Uncomment to let RCMR work again!
           ThermalConduction_AO = scatter
        elseif (ii_loop == 2) then
           ThermalConduction_AO2 = scatter
           !If I want to do the ThermalConduction_O2/0.72 adjustment, do I need to add in
           !this commented line:
          ThermalConduction_AO = ThermalConduction_AO2/0.73
         
        elseif (ii_loop == 3) then
           ThermalConduction_AO = scatter
        end if
     end do
    
  else
          
     ! Send the F10.7 Estimate to all the different processors
     ! AGB Question: why use MPI_SCATTER and not MPI_Bcast?
     if (iProc==0) then
        scattered(:) = output_est
     end if
     
     call MPI_SCATTER(scattered, 1, mpi_double_precision, scatter, 1, &
          mpi_double_precision, 0, iCommGITM, iError)
     
     if(RCMROutType == 'F107') then
        f107_est  = scatter
        f107a_est = f107_est
     else if(RCMROutType == "PHOTOELECTRON") then
        PhotoElectronHeatingEfficiency_est = scatter
     else if(RCMROutType == "EDC") then !ANKIT 6 Oct 2014
        EDC_est = scatter
     end if
     
  endif
  
  !BP
  !do iiProc = 0,nProcs
  !  if (iProc == iiProc) then
  !    write(procString,10) iProc
  !    10    format (I2)
  !    fileNameTC = procString // "_TC.txt"
 ! 
 !     inquire(file = fileNameTC, exist = exist)
 !     if (exist) then
 !       open(34, file = fileNameTC, status = "old", position = "append", &
 !           action = "write")
 !     else
 !       open(34, file = fileNameTC, status = "new", action = "write")
 !     end if
 !    
 !     write(34,*) ThermalConduction_s, ThermalConduction_AO2, ThermalConduction_AO
 !   endif
 ! enddo

 ! close(34) 
 
  call end_timing("RCMR")
end subroutine run_RCMR


! 16 Jan 2017: Following subroutine is RCAC with MIMO Gf Optimization
subroutine RCAC_GfOptim(u_out_GO, theta_out_GO, filter_out, kk, u_in, z_in, y_in)
  use ModRCMR
  use ModBlasLapack, only: LAPACK_getrf, LAPACK_getrs

  implicit none
  integer iii, jjj, rr
  integer, intent(in) :: kk
  DOUBLE PRECISION, intent(in), DIMENSION(lu,1) :: u_in
  DOUBLE PRECISION, intent(in), DIMENSION(lz,1) :: z_in
  DOUBLE PRECISION, intent(in), DIMENSION(ly,1) :: y_in
  DOUBLE PRECISION, intent(out), DIMENSION(lu,1) :: u_out_GO
  DOUBLE PRECISION, intent(out), DIMENSION(ltheta,1) :: theta_out_GO
  DOUBLE PRECISION, intent(out), DIMENSION(lz,nf*lu) :: filter_out

  DOUBLE PRECISION, DIMENSION(1,1) :: lambda_k

  DOUBLE PRECISION, DIMENSION(nf*lu,ltheta) :: Phi_b_rr
  DOUBLE PRECISION, DIMENSION(nf*lu,1) :: U_b_rr
  DOUBLE PRECISION, DIMENSION(lz,kk) :: z_hat_rr
  DOUBLE PRECISION, DIMENSION(lz,ltheta) ::  PHI_filt
  DOUBLE PRECISION, DIMENSION(lz,1) ::  U_filt
  DOUBLE PRECISION, DIMENSION(ltheta,lz) :: PtimesPHI_filt
  DOUBLE PRECISION, DIMENSION(lz,lz) :: Gamma_P, I_lz, Gamma_P_lu, Gamma_P_inv
  
  DOUBLE PRECISION, DIMENSION(lu,lu) :: I_lu

  integer, dimension(lz) :: IPIV_lz
  

  z_hat_rr = 0
  Phi_b_rr = 0
  U_b_rr   = 0
  
  CALL identity(I_lu, lu, 1D0)
  CALL identity(I_lz, lz, 1D0)
  

  ! Put data in control and performance buffers	
  u_window 	= cshift(u_window,-1,2)
  z_window 	= cshift(z_window,-1,2)	
  PHI_window 	= cshift(PHI_window,-lu,1)

  u_window(1:lu,1:1) = u_in !reshape(u_in, shape(u_window(1:lu,1)))
  z_window(1:lz,1:1) = z_in !reshape(z_in, shape(z_window(1:lz,1)))
  
  ! Construct regressor, phi. Size Nc*(lu+lz)+lz
  g_int 				= g_int + z_in
  phi_GO(1:Nc*lu,1) 			= reshape( 1*u_window(:,1:Nc), shape(phi_GO(1:Nc*lu,1)) )
  phi_GO(Nc*lu+1:Nc*lu+Nc*lz,1)	        = reshape( 0*1+1*z_window(:,1:Nc), shape(phi_GO(Nc*lu+1:Nc*lu+Nc*lz,1)) )	
  phi_GO(Nc*lu+Nc*lz+1:,1)		= reshape( 1*g_int, shape(phi_GO(Nc*lu+Nc*lz+1:,1)) )	

  
  
  ! Construct regressor, PHI. Size lu by lu*(Nc*(lu+lz)+lz)
  !CALL kronecker(phi_kron_GO, phi_GO)
  CALL KRONECKER_AB(phi_kron_GO, I_lu, lu, lu, transpose(phi_GO),1, Nc*(lu+lz)+lz)
  
  PHI_window(1:lu,:) = phi_kron_GO !reshape(z_in, shape(z_window(1:lz,1)))
  PHI_b 	= PHI_window(lu+1:(nf+1)*lu,:);
  U_b 		= reshape( u_window(:,1:nf) , shape(U_b) 	)
  Zc  		= reshape(z_window, shape(Zc)) 	!First lz components are z(kk-1)
  Z_GO   	= z_window	                 		!First column is z(kk-1)

  if (kk==1) then
     CALL identity(P_RLS, ltheta,1.0/W_Rtheta)
  endif

  Nf_star = reshape( Nf_Array, shape(Nf_star))
  lambda_k = 0.9999
  
  if (kk> max(max(Nf+1, Nc+1),C_on) )  then
      
      rr = 1
      Phi_b_rr = PHI_window(lu*rr+1:lu*rr+lu*nf,:)
      call vectorize(u_window(:,rr:rr+nf-1 ), lu, nf, U_b_rr)

      PHI_filt = matmul(Nf_star, Phi_b_rr)
      U_filt = matmul(Nf_star, U_b_rr)

      PtimesPhi_filt = matmul(P_RLS, transpose(PHI_filt))

      if (.False.) THEN
         Gamma_P = lambda1*I_lz + matmul(PHI_filt, PtimesPhi_filt)
      else
         Gamma_P = I_lz + matmul(PHI_filt, PtimesPhi_filt)
      end if
      Gamma_P_lu =Gamma_P
      CALL LAPACK_getrf(lz, lz, Gamma_P_lu, lz ,IPIV_lz, INFO)
      CALL identity(Gamma_P_inv, lz, 1D0)
      CALL LAPACK_getrs('n', lz, lz, Gamma_P_lu, lz, IPIV_lz, Gamma_P_inv, lz, INFO)

      !Gamma_P = 1/(lambda1 + matmul(PHI_filt, PtimesPhi_filt))
      !write(*,*) 'Gamma and Gamma inverse are: ', Gamma_P, Gamma_P_inv
      P_RLS = P_RLS - matmul(matmul(PtimesPhi_filt,Gamma_P_inv), transpose(PtimesPhi_filt))
      P_RLS = P_RLS/lambda1

      theta_star = theta_star - matmul(PtimesPhi_filt, matmul(PHI_filt, theta_star)+z_in-U_filt)

      theta_out_GO 	= theta_star
      filter_out 	= Nf_star
      J_GO 		= 0
      Z_hat_c 	        = 0

      u_out_GO = matmul(phi_kron_GO, theta_out_GO)
  else 
      theta_out_GO 	= 0
      filter_out 	   = 0
      u_out_GO 		= u_in
  end if
  
  !BP: Turn these back on! (4/13/2020) 
  call AppendMatrix2File(transpose(u_in), 1, lu,     'u_window11111111.dat')
  call AppendMatrix2File(transpose(z_in), 1, lz,     'z_window11111111.dat')
  call AppendMatrix2File(transpose(u_out_GO), 1, lu, 'u_RCMR1111111111.dat')
  call AppendMatrix2File(transpose(theta_out_GO), 1, ltheta,     'theta11111111111.dat')
  call AppendMatrix2File(filter_out,lz, nf*lu,    'filter1111111111.dat')
  call AppendMatrix2File(lambda_k, 1, 1,             'lambda_k11111111.dat')
  call writeMatrix2File(phi_kron_GO,lu, lu*Nc*(lu+lz)+lu*lz, 'phi_kron_GO11111.dat')
  call writeMatrix2File(A_theta,ltheta, ltheta,     'A_theta111111111.dat')
  call writeMatrix2File(A_thetapinv,ltheta, ltheta, 'A_thetapinv11111.dat')
  call writeMatrix2File(P_RLS,ltheta, ltheta,       'P_RLS11111111111.dat')
  !write(*,*) "z,lambda", z_in, lambda_k
  !write(*,*) "RCMR output is ", kk, u_out_GO, u_in
  !write(*,*) "Optimized estimator coefficients are ", theta_out_GO
  !write(*,*) "Optimized filter is ", filter_out
  
end subroutine RCAC_GfOptim



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
     !write(*,*), A(ii,:)
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


subroutine readSyntheticData(nLines, lu)
  use arrayMod
  use ModGITM
  use ModTime
  
  implicit none
  integer :: nLines, n, i, lu
  logical :: IsFirstTime = .true.
  character(len = 21) :: newFileName = "syntheticDataFile.txt"

  if (IsFirstTime ) then
    n = nLines
     open (unit = 31, file = newFileName, status = "old")
  
     ! DATA FORMAT - first line should be
     !1833  2002-12-31 00:05:00 1.352269e-12 2.84807e-12 2.1001695e-12 0.00036 0.00056 0.69 
  
     !write(*,*) "Reading synthetic data..."
     do i = 1, n
        if (lu == 1) then
          read(31, *) satIndex(i), date(i), time(i), satLon(i), satLat(i), &
                      satAlt(i), satMeanDensity(i)
        else if (lu == 2) then
          read(31, *) satIndex(i), date(i), time(i), satLon(i), satLat(i), &
                      satAlt(i), satMeanDensity(i), sat2Lon(i), sat2Lat(i), &
                      sat2Alt(i), sat2Density(i)
        end if

        read(date(i)(1:4), '(i4)') satYear(i)
        read(date(i)(6:7), '(i2)') satMonth(i)
        read(date(i)(9:10), '(i2)') satDay(i)
        
        read(time(i)(1:2), '(i2)') satHour(i)
        read(time(i)(4:5), '(i2)') satMinute(i)
        read(time(i)(7:8), '(i2)') satSecond(i)
     end do
  
     close(31)

     !if (iProc == 0) then
     !  write(*,*) date(1), date(n), satLon(1), satLon(n)
     !end if

     IsFirstTime = .false.

  endif

  return

end subroutine readSyntheticData

subroutine printDensity(nLines, diff, diff2, iiBlock, lu)
  use ModTime
  use ModGITM
  use arrayMod
  use ModMPI
  !use ModRCMR, only: MaxRho_R, MinRho_R, MeanRho_R 
  real :: LatFind, LonFind, AltFind
  real :: rLon, rLat, rAlt
  real :: gitmDensity
  integer :: iiLat, iiLon, iiAlt
  integer :: loc = -99
  integer :: i, lu
  integer :: n
  integer, intent(out) :: iiBlock
  !integer, intent(out) :: returnProc, returnProc2
  logical :: exist
  real, intent(out) :: diff, diff2
  real :: percentDifference

  n = nLines  
  diff = -1e32
  diff2 = -1e32
    
  do i = 1, n
    if ((satYear(i) .eq. iTimeArray(1))   .and. &
        (satMonth(i) .eq. iTimeArray(2))  .and. &
        (satDay(i) .eq. iTimeArray(3))    .and. &
        (satHour(i) .eq. iTimeArray(4))   .and. &
        (satMinute(i) .eq. iTimeArray(5)) .and. &
        (satSecond(i) .eq. iTimeArray(6))) then

      loc = i
      exit
    endif
  end do  

  if (loc .eq. -99 .and. i .eq. n) then
    write(*,*) "ERROR: Time not found in satellite data"
    return
  endif
  
  !write(*,*) "Location found", loc
  ! ------------------------------------------------------
  !For first satellite (CHAMP)
  LonFind = satLon(loc)*pi/180.0
  LatFind = satLat(loc)*pi/180.0
  call LocationIndex(LonFind, LatFind, iiBlock, iiLon, iiLat, rLon, rLat)

  !Finds the satellite location and finds the corresponding mass density
  !provided by GITM and prints it if it's within 5 seconds of CHAMP
  !write(*,*) iiBlock
    
  if (iiBlock .eq. 1) then
    !if (iTimeArray(5) == 10) then
    !  write(*,*) "iProc", iProc
    !endif
    if (satHour(loc) .eq. iTimeArray(4) .and. &
        satMinute(loc) .eq. iTimeArray(5) .and. &
        satSecond(loc) .eq. iTimeArray(6)) then
      write(*,*) "GITM's time:", iTimeArray(2), "/", iTimeArray(3), "/", &
      iTimeArray(1), iTimeArray(4), ":", iTimeArray(5), ":", iTimeArray(6)
      
      !write(*,*) "Data file time:", satMonth(loc), "/", &
      !           satDay(loc), "/", satYear(loc), satHour(loc), ":", &
      !           satMinute(loc), ":", satSecond(loc)

      !write(*,*) "At this time, CHAMP location is:", satLon(loc), &
      !           satLat(loc), satAlt(loc) 
      

      !write(*,*) "iiBlock:", iiBlock
      !write(*,*) "iiLon:", iiLon
      !write(*,*) "iiLat:", iiLat
      !write(*,*) "rLon:", rLon
      !write(*,*) "rLat:", rLat

      AltFind = satAlt(loc)*1000.0
      call BlockAltIndex(AltFind, iiBlock, iiLon, iiLat, iiAlt, rAlt)
      !write(*,*) "iiAlt:", iiAlt
      
      
      
      !write(*,*) "GITM's found location:", Longitude(iiLon, iiBlock)*180/pi, &
      !           Latitude(iiLat, iiBlock)*180/pi, &
      !           Altitude_GB(iiLon, iiLat, iiAlt, iiBlock)      

      !GITM's density at satellite's location
      gitmDensity = rLon*rLat*rAlt*Rho(iiLon, iiLat, iiAlt, iiBlock)+ &
                    (1-rLon)*rLat*rAlt*Rho(iiLon+1, iiLat, iiAlt, iiBlock) + &
                    rLon*(1-rLat)*rAlt*Rho(iiLon, iiLat+1, iiAlt, iiBlock) + &
                    (1-rLon)*(1-rLat)*rAlt*Rho(iiLon+1, iiLat+1, iiAlt, iiBlock) + &
                    rLat*rLon*(1-rAlt)*Rho(iiLon, iiLat, iiAlt+1, iiBlock)+ &
                    (1-rLon)*rLat*(1-rAlt)*Rho(iiLon+1, iiLat, iiAlt+1, iiBlock) + &
                    rLon*(1-rLat)*(1-rAlt)*Rho(iiLon, iiLat+1, iiAlt+1, iiBlock) + &
                    (1-rLon)*(1-rLat)*(1-rAlt)*Rho(iiLon+1, iiLat+1, iiAlt+1, iiBlock)

      !write(*,*) "GITM's density:", gitmDensity
      !write(*,*) "Truth data density:", satMeanDensity(loc)
      
      !Populate the difference array, z
      diff = satMeanDensity(loc) - gitmDensity
      !write(*,*) "Diff in function:", diff
      !inquire(file = "RCMR_locations.txt", exist = exist)
      !if (exist) then
      !  open(33, file = "RCMR_locations.txt", status = "old", position = "append", &
      !       action = "write")
      !else
      !  open(33, file = "RCMR_locations.txt", status = "new", action = "write")
      !end if

      !111 FORMAT(I4, 1X, I2, 1X, I2, 1X, &
      !           I2, 1X, I2, 1X, I2, 1X, &
      !           F8.2, 1X, F8.2, 1X, F8.2)
      !write(33,111) iTimeArray(1), iTimeArray(2), &
      !              iTimeArray(3), iTimeArray(4), &
      !              iTimeArray(5), iTimeArray(6), &
      !              satLon(loc), satLat(loc), satAlt(loc)

      !close(33)
      
    endif
  endif 

  ! ------------------------------------------------------                                   
  !For second satellite GRACE 
  if (lu == 2) then                                  
    LonFind = sat2Lon(loc)*pi/180.0
    LatFind = sat2Lat(loc)*pi/180.0
    call LocationIndex(LonFind, LatFind, iiBlock, iiLon, iiLat, rLon, rLat)
    if (iiBlock .eq. 1) then
      if (satHour(loc) .eq. iTimeArray(4) .and. &
          satMinute(loc) .eq. iTimeArray(5) .and. &
          satSecond(loc) .eq. iTimeArray(6)) then

        !write(*,*) "At this time, GRACE location is:", sat2Lon(loc), &
        !           sat2Lat(loc), sat2Alt(loc)

        !write(*,*) "iiBlock:", iiBlock
        !write(*,*) "iiLon:", iiLon
        !write(*,*) "iiLat:", iiLat
        !write(*,*) "rLon:", rLon
        !write(*,*) "rLat:", rLat

        AltFind = sat2Alt(loc)*1000.0
        call BlockAltIndex(AltFind, iiBlock, iiLon, iiLat, iiAlt, rAlt)
        !write(*,*) "iiAlt:", iiAlt

        !write(*,*) "GITM's found location:", Longitude(iiLon, iiBlock)*180/pi, &
        !           Latitude(iiLat, iiBlock)*180/pi, &
        !           Altitude_GB(iiLon, iiLat, iiAlt, iiBlock)
      
        !GITM's density at satellite's location                                                  
        gitmDensity = rLat*rLon*rAlt*Rho(iiLon, iiLat, iiAlt, iiBlock)+ &
                      (1-rLon)*rLat*rAlt*Rho(iiLon+1, iiLat, iiAlt, iiBlock) + &
                      rLon*(1-rLat)*rAlt*Rho(iiLon, iiLat+1, iiAlt, iiBlock) + &
                      (1-rLon)*(1-rLat)*rAlt*Rho(iiLon+1, iiLat+1, iiAlt, iiBlock) + &
                      rLat*rLon*(1-rAlt)*Rho(iiLon, iiLat, iiAlt+1, iiBlock)+ &
                      (1-rLon)*rLat*(1-rAlt)*Rho(iiLon+1, iiLat, iiAlt+1, iiBlock) + &
                      rLon*(1-rLat)*(1-rAlt)*Rho(iiLon, iiLat+1, iiAlt+1, iiBlock) + &
                      (1-rLon)*(1-rLat)*(1-rAlt)*Rho(iiLon+1, iiLat+1, iiAlt+1, iiBlock)

        !write(*,*) "GITM's density:", gitmDensity
        !write(*,*) "Truth data density:", sat2Density(loc)

        !Populate the difference array, z                                                        
        diff2 = sat2Density(loc) - gitmDensity
      !returnProc2 = iProc
      endif
    endif
  endif

end subroutine printDensity

subroutine arraySub(nLines)
  use arrayMod

  !write(*,*) "Inside arraySub()"
  allocate(satIndex(nLines), &
             satYear(nLines), &
             satMonth(nLines), &
             satDay(nLines), &
             satHour(nLines), &
             satMinute(nLines), &
             satSecond(nLines), &
             satMeanDensity(nLines), &
             sat2Density(nLines), &
             satLat(nLines), sat2Lat(nLines), &
             satLon(nLines), sat2Lon(nLines), &
             satAlt(nLines), sat2Alt(nLines), &
             time(nLines), date(nLines))

end subroutine arraySub
