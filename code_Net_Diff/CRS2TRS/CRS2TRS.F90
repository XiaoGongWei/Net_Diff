! =========================CRS2TRS=====================
! This subroutine is used to calculate the transfer matrix between CRS
! & and TRS given a MJD time and the EOP parameters (including Xp,Yp
! & dUT1, DX00, DY00)
! The reference frame is IAU 2000A

! Adapt from SOFA

! Inputs:
!     Leap_Sec  leap seconds
!     MJD        Modified Julian Day, in GPST
!     Xp
!     Yp
!     dUT1       UT1-UTC
!     DX00     
!     DY00

! Outputs:
!     RC2T       rotation matrix from CRS (Celestrial Reference system) 
!                    to TRS (Terrestial Reference System)
!     RT2C       rotation matrix from TRS (Terrestial Reference System)
!                    to CRS (Celestrial Reference system) 

! Written by: Yize Zhang
! =====================End of Header====================

subroutine CRS2TRS(Leap_sec, MJD,Xp,Yp,dUT1,DX00,DY00,RC2T,RT2C)
implicit none
    ! Out:
        real(8) :: RC2T(3,3), RT2C(3,3)

    ! Variables
    !Arcsecond to radians
    DOUBLE PRECISION AS2R
    PARAMETER (AS2R=4.848136811095359935899141D-6)
    
    ! 2Pi
    DOUBLE PRECISION D2PI
    PARAMETER ( D2PI = 6.283185307179586476925287D0 )
    
    real(8) :: Leap_sec, MJD
    INTEGER J,I  ! IY, IM, ID, IH, MIN, 
    DOUBLE PRECISION XP, YP, DUT1  ! SEC, 
    DOUBLE PRECISION DDP80, DDE80, DX00, DY00, DX06, DY06
    DOUBLE PRECISION DJMJD0, DATE, TIME, UTC, CurrDAT
    DOUBLE PRECISION TAI, TT, TUT, UT1, RP(3,3), DP80, DE80
    DOUBLE PRECISION DPSI, DEPS, EPSA, RN(3,3), RNPB(3,3)
    DOUBLE PRECISION EE, GST, RC2TI(3,3), RPOM(3,3)
    DOUBLE PRECISION RC2IT(3,3), X, Y, S
    DOUBLE PRECISION RC2I(3,3), ERA, DP00, DE00, RB(3,3)
    DOUBLE PRECISION RPB(3,3), V1(3), V2(3), DDP00, DDE00
    
    DOUBLE PRECISION iau_OBL80, iau_EQEQ94, iau_ANP, iau_GMST82
    DOUBLE PRECISION iau_ERA00, iau_SP00, iau_EE00, iau_GMST00
    DOUBLE PRECISION iau_S06
   
   
    ! Polar motion (arcsec->radians).
        XP = XP * AS2R
        YP = YP * AS2R
    
    ! CIP offsets wrt IAU 2000A (mas->radians).
        DX00 =DX00 * AS2R
        DY00 = DY00 * AS2R
        
    ! TT (MJD).
        DJMJD0=2400000.5d0
!        CALL iau_CAL2JD ( IY, IM, ID, DJMJD0, DATE, J )
!        TIME = ( 60D0*(60D0*DBLE(IH) + DBLE(MIN)) + SEC ) / 86400D0
!        UTC = DATE + TIME
!        call DAT(MJD, CurrDAT) ! GPST or UTC is both OK, it doesn't effect the value of CurrDAT
!        CALL iau_DAT ( IY, IM, ID, TIME, CurrDAT, J )
        !TAI = UTC + CurrDAT/86400D0
        TAI = MJD + 19.d0/86400.d0
        UTC = MJD - Leap_Sec/86400.d0
!        UTC=TAI - CurrDAT/86400.d0
        DATE=floor(UTC)      ! Used in subroutine iau_ERA00
        TIME=UTC-DATE   ! the decimal part of UTC
        TT = TAI + 32.184D0/86400D0   ! Used in subroutine  iau_XYS00A  and iau_POM00
        TUT = TIME + DUT1/86400D0   ! the decimal part of UTC, Used in subroutine  iau_ERA00  
!        UT1 = UTC+DUT1/86400d0 ! =DATE + TUT       ! Not used here
       
   ! CIP and CIO, IAU 2000A.
        CALL iau_XYS00A ( DJMJD0, TT, X, Y, S )
     ! Add CIP corrections.
        X = X + DX00
        Y = Y + DY00
     ! GCRS to CIRS matrix.
        CALL iau_C2IXYS ( X, Y, S, RC2I )
     ! Earth rotation angle
        ERA=iau_ERA00( DJMJD0+Date, TUT )
     ! Form celestial-terrestrial matrix (no polar motion yet)
        call iau_CR( RC2I, RC2TI )
        call iau_RZ( ERA, RC2TI )
     !  Polar motion matrix (TIRS --> ITRS, IERS 2003 )
        call iau_POM00( Xp, Yp, iau_SP00( DJMJD0, TT), RPOM )
     !  Form celestial-terrestrial matrix (including polar motion )
        call iau_RXR ( RPOM, RC2TI, RC2IT )
        
        RC2T=RC2IT
        RT2C=transpose(RC2T)
end subroutine
