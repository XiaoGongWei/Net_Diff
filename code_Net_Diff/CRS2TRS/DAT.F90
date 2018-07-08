! =====================DAT====================
! For a given UTC date, calculate delta(AT) = TAI-UTC.

! Adapt from SOFA, iau_DAT.for

! Note:
! 1)  A new version of this routine must be produced whenever 
!      & a new leap second is announced. 

! 2)  Since leap seconds is introduced after 1972,1,1, this subroutine is only
!      &  validity after 1972,1,1 (In fact, It's enough for GNSS data processing).
!      &  For a earlier time, please refer to SOFA, iau_DAT.for

! Input:
!      MJD      i    Modified Julian Day

! Output:
!     dDAT      d    TAI minus UTC, in seconds

! Called by:
!             C2R.F90

! Written by: Yize Zhang
! =================End of Header==============

subroutine DAT(MJD, dDAT)
implicit none
    ! intent in
    real(8):: MJD
    ! Intent out
    real(8) :: dDAT
    
      ! Local variables
      ! Release year for this version of iau_DAT
      integer, parameter :: IYV=2009
      ! Number of Delta(AT) changes (increase by 1 for each new leap second)
      integer,parameter :: NDat=28
      ! Dates (year, month) on which new Delta(AT) came into force
!      integer IDATE(2,NDAT)
      real(8) :: MJDS(NDAT)
     ! New Delta(AT) which came into force on the given dates
      real(8) DATS(NDAT)
      
      integer ::  i
      ! Dates( in MJD) and Delta(AT)s
!      DATA (IDATE(I, 1),I=1,2),DATS(1)  / 1960,  1,  1.4178180D0 /
!      DATA (IDATE(I, 2),I=1,2),DATS(2)  / 1961,  1,  1.4228180D0 /
!      DATA (IDATE(I, 3),I=1,2),DATS(3)  / 1961,  8,  1.3728180D0 /
!      DATA (IDATE(I, 4),I=1,2),DATS(4)  / 1962,  1,  1.8458580D0 /
!      DATA (IDATE(I, 5),I=1,2),DATS(5)  / 1963, 11,  1.9458580D0 /
!      DATA (IDATE(I, 6),I=1,2),DATS(6)  / 1964,  1,  3.2401300D0 /
!      DATA (IDATE(I, 7),I=1,2),DATS(7)  / 1964,  4,  3.3401300D0 /
!      DATA (IDATE(I, 8),I=1,2),DATS(8)  / 1964,  9,  3.4401300D0 /
!      DATA (IDATE(I, 9),I=1,2),DATS(9)  / 1965,  1,  3.5401300D0 /
!      DATA (IDATE(I,10),I=1,2),DATS(10) / 1965,  3,  3.6401300D0 /
!      DATA (IDATE(I,11),I=1,2),DATS(11) / 1965,  7,  3.7401300D0 /
!      DATA (IDATE(I,12),I=1,2),DATS(12) / 1965,  9,  3.8401300D0 /
!      DATA (IDATE(I,13),I=1,2),DATS(13) / 1966,  1,  4.3131700D0 /
!      DATA (IDATE(I,14),I=1,2),DATS(14) / 1968,  2,  4.2131700D0 /
      DATA MJDS(1),DATS(1) / 41317D0, 10D0 /  ! / 1972,  1, 10D0 /
      DATA MJDS(2),DATS(2) / 41499D0, 11D0 /  ! / 1972,  7, 11D0 /
      DATA MJDS(3),DATS(3) / 41683D0, 12D0 /  ! / 1973,  1, 12D0 /
      DATA MJDS(4),DATS(4) / 42048D0, 13D0 /  ! / 1974,  1, 13D0 /
      DATA MJDS(5),DATS(5) / 42413D0, 14D0 /  ! / 1975,  1, 14D0 /
      DATA MJDS(6),DATS(6) / 42778D0, 15D0 /  ! / 1976,  1, 15D0 /
      DATA MJDS(7),DATS(7) / 43144D0, 16D0 /  ! / 1977,  1, 16D0 /
      DATA MJDS(8),DATS(8) / 43509D0, 17D0 /  ! / 1978,  1, 17D0 /
      DATA MJDS(9),DATS(9) / 43874D0, 18D0 /  ! / 1979,  1, 18D0 /
      DATA MJDS(10),DATS(10) / 44239D0, 19D0 /  ! / 1980,  1, 19D0 /
      DATA MJDS(11),DATS(11) / 44786D0, 20D0 /  ! / 1981,  7, 20D0 /
      DATA MJDS(12),DATS(12) / 45151D0, 21D0 /  ! / 1982,  7, 21D0 /
      DATA MJDS(13),DATS(13) / 45516D0, 22D0 /  ! / 1983,  7, 22D0 /
      DATA MJDS(14),DATS(14) / 46247D0, 23D0 /  ! / 1985,  7, 23D0 /
      DATA MJDS(15),DATS(15) / 47161D0, 24D0 /  ! / 1988,  1, 24D0 /
      DATA MJDS(16),DATS(16) / 47892D0, 25D0 /  ! / 1990,  1, 25D0 /
      DATA MJDS(17),DATS(17) / 48257D0, 26D0 /  ! / 1991,  1, 26D0 /
      DATA MJDS(18),DATS(18) / 48804D0, 27D0 /  ! / 1992,  7, 27D0 /
      DATA MJDS(19),DATS(19) / 49169D0, 28D0 /  ! / 1993,  7, 28D0 /
      DATA MJDS(20),DATS(20) / 49534D0, 29D0 /  ! / 1994,  7, 29D0 /
      DATA MJDS(21),DATS(21) / 50083D0, 30D0 /  ! / 1996,  1, 30D0 /
      DATA MJDS(22),DATS(22) / 50630D0, 31D0 /  ! / 1997,  7, 31D0 /
      DATA MJDS(23),DATS(23) / 51179D0, 32D0 /  ! / 1999,  1, 32D0 /
      DATA MJDS(24),DATS(24) / 53736D0, 33D0 /  ! / 2006,  1, 33D0 /
      DATA MJDS(25),DATS(25) / 54832D0, 34D0 /  ! / 2009,  1, 34D0 /
      DATA MJDS(26),DATS(26) / 56109D0, 35D0 /  ! / 2012,  7, 35D0 /
      DATA MJDS(27),DATS(27) / 57204D0, 36D0 /  ! / 2015,  7, 36D0 /
      DATA MJDS(28),DATS(28) / 57754D0, 37D0 /  ! / 2017,  1, 37D0 /

      ! Find the most recent table entry.
      do i=NDAT,1,-1
          if (MJD >=MJDS(i) ) then
              dDAT=DATS(i)
              exit
          end if
      end do
      
    end subroutine