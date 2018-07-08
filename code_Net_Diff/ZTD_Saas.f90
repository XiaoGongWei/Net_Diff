! -------------------------------------------------------------------------
! Purpose  : This Function is to Compute ZTD using Saastamounien model
!            
!            
! 
! Author   : Junping Chen   17.Dec.2008
! Reference: EPOS
! Changes  : 
!     Changed by Yize Zhang 5 Sempember,2013, used for personal Visual Studio Fortran   
!############################ end_of_header ##############################################

 subroutine ztd_saas(p,t,rh,rlat,height,zdd,zwd)
! ------------------------------------------------------------------------------

! use MOD_define_type   Changed by Yize Zhang
 implicit none

! real (r8_) :: ePressure !water vapour pressure in hPa
 real (8) :: p      !  total pressure, mbars
 real (8) :: t      !  temperature, degree C
 real (8) :: rh     !  relative humity
 real (8) :: rlat   !  geodetic latitude, radians
 real (8) :: height !  station height, m
 real (8) :: zdd    ! zenith dry delay in meters
 real (8) :: zwd    ! zenith wet delay in meters
!
! local variable
 real (8) :: e, tk, fact

 e=rh  ! gpt2 model
 if (rh.eq.0.6d0) e = rh * 6.11d0 * 10.d0 ** (7.5d0 * t / (t + 2.373d+02))  ! other model
!
!! Temperature in Kelvin
 tk = t + 273.15
!
!! ellipsoidal/elevation-dependent variation 
 fact = 1.d0-0.266d-2 * cos(2.0d0*rlat)-0.28d-3 * height*1d-3
!
!! Zenith dry delay, meters
 zdd = 0.22768d-2 * p /fact  
!
!! Zenith wet delay, meters
 zwd = 0.22768d-2 * (.1255d+4/tk+.5d-1) * e/fact

 return
 end
