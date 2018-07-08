C     Last change:  OO   20 Nov 2000   11:22 am
      SUBROUTINE KLOBUCHAR(fi,lambda,elev,azimuth,tow,alfa,beta,dIon1)
C     ==================================================================
C     Subroutine for computing an Ionospheric range correction for the *
C     GPS L1 frequency from the parameters broadcasted in the GPS      *
C     Navigation Message.                                              *
C     ==================================================================
C     References:                                                      *
C     Klobuchar, J.A., (1996) "Ionosphercic Effects on GPS", in        *
C       Parkinson, Spilker (ed), "Global Positioning System Theory and *
C       Applications, pp.513-514.                                      *
C     ICD-GPS-200, Rev. C, (1997), pp. 125-128                         *
C     NATO, (1991), "Technical Characteristics of the NAVSTAR GPS",    *
C       pp. A-6-31   -   A-6-33                                        *
C     ==================================================================
C     Author : Ola Ovstedal, Department of Mapping Sciences            *
C                            Agricultural University of Norway         *
C     Date   : 06.03.2000                  Last modified : 20.11.2000  *
C     ==================================================================
C     Input :                                                          *
C       fi            : Geodetic latitude of receiver          (deg)   *
C       lambda        : Geodetic longitude of receiver         (deg)   *
C       elev          : Elevation angle of satellite           (deg)   *
C       azimuth       : Geodetic azimuth of satellite          (deg)   *
C       tow           : Time of Week                           (sec)   *
C       alfa(4)       : The coefficients of a cubic equation           *
C                       representing the amplitude of the vertical     *
C                       dalay (4 coefficients - 8 bits each)           *
C       beta(4)       : The coefficients of a cubic equation           *
C                       representing the period of the model           *
C                       (4 coefficients - 8 bits each)                 *
C     Output:                                                          *
C       dIon1         : Ionospheric slant range correction for         *
C                       the L1 frequency                       (metre) *
C     ==================================================================
C     Called by       :                                                *
C                                                                      *
C     Subroutines called:                                              *
C           Internal:                                                  *
C           External:                                                  *
C     Include files   :                                                *
C     ==================================================================
C     Changes:                                                         *
C     20.11.2000 : Different constants assigned in subroutine so that  *
C                  the subroutine is self-contained (for GPS-TOOLBOX)  *
C     ==================================================================
      IMPLICIT NONE
C
      REAL*8 alfa(4),beta(4),fi,lambda,elev,e,azimuth,a,tow,t,
     +       dIon1,psi,lat_i,long_i,lat_m,sf,PER,AMP,x,deg2semi,
     +       semi2rad,deg2rad,pi,c


      pi       =  3.1415926535898D0             ! WGS-84 value of pi
      c        =  2.99792458D8                  ! speed of light
      deg2semi =  1.D0/180.D0                   ! degees to semisircles
      semi2rad =  pi                            ! semisircles to radians
      deg2rad  =  pi/180.D0                     ! degrees to radians

      a = azimuth*deg2rad                       ! asimuth in radians
      e = elev*deg2semi                         ! elevation angle in
                                                ! semicircles

      psi = 0.0137D0 / (e+0.11D0) - 0.022D0     ! Earth Centered angle

      lat_i = fi*deg2semi + psi*DCOS(a)         ! Subionospheric lat
      IF(lat_i .GT. 0.416D0) THEN
        lat_i = 0.416D0
      ELSEIF(lat_i .LT. -0.416) THEN
        lat_i = -0.416D0
      ENDIF

                                                ! Subionospheric long
      long_i = lambda*deg2semi + (psi*DSIN(a)/DCOS(lat_i*semi2rad))

                                                ! Geomagnetic latitude
      lat_m = lat_i + 0.064D0*DCOS((long_i-1.617D0)*semi2rad)

      t = 4.32D4*long_i + tow                    ! Local time
      t = DMOD(t,86400.D0)                      ! Seconds of day
      IF(t .GT. 86400.D0)  t = t - 86400.D0
      IF(t .LT. 0.D0)      t = t + 86400.D0

      sf = 1.D0 + 16.D0*(0.53D0-e)**3           ! Slant factor

                                                ! Period of model
      PER = beta(1) + beta(2)*lat_m + beta(3)*lat_m**2 +beta(4)*lat_m**3

      IF(PER .LT. 72000.D0) PER = 72000.D0

      x = 2.D0*pi*(t-50400.D0)  /  PER          ! Phase of the model
                                                ! (Max at 14.00 =
                                                ! 50400 sec local time)

                                                ! Amplitud of the model
      AMP = alfa(1) + alfa(2)*lat_m + alfa(3)*lat_m**2 +alfa(4)*lat_m**3
      IF(AMP .LT. 0.D0) AMP = 0.D0

                                                ! Ionospheric corr.
      IF(DABS(x) .GT. 1.57D0) THEN
        dIon1 = sF * (5.D-9)
      ELSE
        dIon1 = sF * (5.D-9 + AMP*(1.D0 - x*x/2.D0 + x*x*x*x/24.D0))
      ENDIF

      dIon1 = c * dIon1
      RETURN
      END
