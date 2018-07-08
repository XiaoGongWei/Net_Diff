c****************************************************************************************
       SUBROUTINE TIDE(XSTA,iYR,iMONTH,iDAY,FHR,XSUN,XMON,DXTIDE)
C  
C PURPOSE    :  COMPUTATION OF TIDAL CORRECTIONS OF STATION DISPLACEMENTS  
C               CAUSED BY LUNAR AND SOLAR GRAVITATIONAL ATTRACTION  
C               (SEE IERS STANDARDS 1996)  
C         STEP 1 (HERE GENERAL DEGREE 2 AND 3 CORRECTIONS +  
C                CALL ST1IDIU + CALL ST1ISEM + CALL ST1L1) 
C         + STEP 2 (CALL STEP2DIU + CALL STEP2LON + CALL STEP2IDIU)  
C IT HAS BEEN DECIDED THAT THE STEP 3 NON-CORRECTION FOR PERMANENT TIDE 
C WOULD NOT BE APPLIED IN ORDER TO AVOID JUMP IN THE REFERENCE FRAME 
C (THIS STEP 3 MUST ADDED IN ORDER TO GET THE NON-TIDAL STATION POSITION  
C AND TO BE CONFORMED WITH THE IAG RESOLUTION.)  
C  
C    INPUT :  XSTA(I),I=1,2,3: GEOCENTRIC POSITION OF THE STATION (ITRF,
C                              CO-ROTATING FRAME) -- UNITS = METERS
C             XSUN(I),I=1,2,3: GEOC. POSITION OF THE SUN (ECEF FRAME) --
C                              UNITS = METERS
C             XMON(I),I=1,2,3: GEOC. POSITION OF THE MOON (ECEF FRAME) --
C                              UNITS = METERS
C             IYR : YEAR (UTC TIMESCALE)
C             IMONTH : MONTH (UTC TIMESCALE)
C             IDAY : DAY (UTC TIMESCALE)
C             FHR=hr+zmin/60.+sec/3600. : HR IN THE DAY 
C   OUTPUT :  DXTIDE(I),I=1,2,3: DISPLACEMENT VECTOR (GEOCENTRIC ITRF FRAME) --
C                                UNITS = METERS
C  
C SUBROUTINES CALLED  :  SPROD  
C                        ST1IDIU 
C                        ST1ISEM 
C                        ST1L1 
C                        STEP2DIU  
C                        STEP2LON 
C                        STEP2ILON  
C  
C AUTHOR IERS 1996 :  V. DEHANT, S. MATHEWS AND J. GIPSON
C    (TEST BETWEEN TWO SUBROUTINES) 
C AUTHOR IERS 2000 :  V. DEHANT AND S. MATHEWS
C    (TEST IN THE BERNESE PROGRAM BY C. BRUYNINX)
C  
C CREATED    :  96/03/23              LAST MODIFIED :  00/05/17 14:10  
C UPDATED    :2006/02/06  HEADER COMMENTS MODIFIED TO CLARIFY INPUT/OUTPUT
C                         UNITS AND SYSTEMS BY JIM RAY
C UPDATED    :2006/02/06  SUBROUTINE DUTC MODIFIED FOR LEAP SECOND ON
C                         2006.0 AND TO CORRECT do 5 i=1,87 from 84 to 87
C                         BY JIM RAY
C  
C  
C  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  
      DOUBLE PRECISION XSTA(3),XSUN(3),XMON(3),DXTIDE(3),XCORSTA(3)
      DOUBLE PRECISION H20,L20,H3,L3,H2,L2
      DOUBLE PRECISION mass_ratio_sun,mass_ratio_moon
C  
C NOMINAL SECOND DEGREE AND THIRD DEGREE LOVE NUMBERS AND SHIDA NUMBERS  
C ---------------------------------------------------------------------  
      DATA H20/0.6078D0/,L20/0.0847D0/,H3/0.292D0/,L3/0.015D0/
C  
C SCALAR PRODUCT OF STATION VECTOR WITH SUN/MOON VECTOR  
C -----------------------------------------------------  
      CALL SPROD(XSTA,XSUN,SCS,RSTA,RSUN)  
      CALL SPROD(XSTA,XMON,SCM,RSTA,RMON)  
      SCSUN=SCS/RSTA/RSUN  
      SCMON=SCM/RSTA/RMON
C   
C COMPUTATION OF NEW H2 AND L2  
C ----------------------------  
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA
      H2=H20-0.0006*(1.-3./2.*COSPHI**2)
      L2=L20+0.0002*(1.-3./2.*COSPHI**2)  
C
C P2-TERM  
C -------  
      P2SUN=3.*(H2/2.-L2)*SCSUN**2-H2/2.  
      P2MON=3.*(H2/2.-L2)*SCMON**2-H2/2.  
C  
C P3-TERM  
C -------  
      P3SUN=5./2.*(H3-3.*L3)*SCSUN**3+3./2.*(L3-H3)*SCSUN
      P3MON=5./2.*(H3-3.*L3)*SCMON**3+3./2.*(L3-H3)*SCMON

C  
C TERM IN DIRECTION OF SUN/MOON VECTOR  
C ------------------------------------  
      X2SUN=3.*L2*SCSUN  
      X2MON=3.*L2*SCMON  
      X3SUN=3.*L3/2.*(5.*SCSUN**2-1.)  
      X3MON=3.*L3/2.*(5.*SCMON**2-1.)
C  
C FACTORS FOR SUN/MOON  
C --------------------  
      MASS_RATIO_SUN=332945.943062d0
      MASS_RATIO_MOON=0.012300034d0
      RE =6378136.55d0
      FAC2SUN=MASS_RATIO_SUN*RE*(RE/RSUN)**3
      FAC2MON=MASS_RATIO_MOON*RE*(RE/RMON)**3
      FAC3SUN=FAC2SUN*(RE/RSUN)
      FAC3MON=FAC2MON*(RE/RMON)
C  
C TOTAL DISPLACEMENT  
C ------------------  
      DO 10 I=1,3  
      DXTIDE(I)=FAC2SUN*( X2SUN*XSUN(I)/RSUN + P2SUN*XSTA(I)/RSTA ) +
     1            FAC2MON*( X2MON*XMON(I)/RMON + P2MON*XSTA(I)/RSTA ) +  
     2            FAC3SUN*( X3SUN*XSUN(I)/RSUN + P3SUN*XSTA(I)/RSTA ) +   
     3            FAC3MON*( X3MON*XMON(I)/RMON + P3MON*XSTA(I)/RSTA )  
10    CONTINUE  
      call zero_vec8(xcorsta)
C  
C CORRECTIONS FOR THE OUT-OF-PHASE PART OF LOVE NUMBERS (PART H_2^(0)I  
C            AND L_2^(0)I )  
C FIRST, FOR THE DIURNAL BAND       

      CALL ST1IDIU(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)
      DO 11 I=1,3
      DXTIDE(I)=DXTIDE(I)+XCORSTA(I)  
11    CONTINUE
C  
C SECOND, FOR THE SEMI-DIURNAL BAND       
C 
      CALL ST1ISEM(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)
      DO 12 I=1,3
      DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
12    CONTINUE  
C  
C CORRECTIONS FOR THE LATITUDE DEPENDENCE OF LOVE NUMBERS (PART L^(1) )  
C   
      CALL ST1L1(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)
      DO 13 I=1,3  
      DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
13    CONTINUE    
C  
C CONSIDER CORRECTIONS FOR STEP 2  
C  
C CORRECTIONS FOR THE DIURNAL BAND:  
C 
C  FIRST, WE NEED TO KNOW THE DATE CONVERTED IN JULIAN CENTURIES 
c        
c   1) CALL THE SUBROUTINE COMPUTING THE JULIAN DATE 
c 
c      EXPRESSION OF THE HOURS, MINUTES AND SECONDES IN FRACTION OF DAY        
c 
      djuld=fjldy(iyr,imonth,iday,fhr)
C convert to centuries.
      T=(DJULD-2451545.)/36525.
c 
c   2) CALL THE SUBROUTINE COMPUTING THE CORRECTION OF UTC TIME  
c 
      CALL DUTC(IYR,IMONTH,IDAY,DTT) 
      fhr=fhr+dtt/3600.
C  
C  SECOND, WE CAN CALL THE SUBROUTINE STEP2DIU, FOR THE DIURNAL BAND CORRECTIONS,
C   (in-phase and out-of-phase frequency dependence):
C  
      CALL STEP2DIU(XSTA,FHR,T,XCORSTA)  
      DO 14 I=1,3  
      DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
14    CONTINUE  
C  
C  CORRECTIONS FOR THE LONG-PERIOD BAND,
C   (in-phase and out-of-phase frequency dependence):  
C
      call STEP2LON(XSTA,FHR,T,XCORSTA)
      DO 15 I=1,3  
      DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
15    CONTINUE  
C    
C CONSIDER CORRECTIONS FOR STEP 3  
C  
C UNCORRECT FOR THE PERMANENT TIDE  
C  
C      PI=3.141592654
C      SINPHI=XSTA(3)/RSTA  
C      COSPHI=dsqrt(XSTA(1)**2+XSTA(2)**2)/RSTA
C      COSLA=XSTA(1)/COSPHI/RSTA  
C      SINLA=XSTA(2)/COSPHI/RSTA  
C      DR=-DSQRT(5./4./PI)*H2*0.31460*(3./2.*SINPHI**2-0.5)
C      DN=-DSQRT(5./4./PI)*L2*0.31460*3.*COSPHI*SINPHI
C      DXTIDE(1)=DXTIDE(1)-DR*COSLA*COSPHI+DN*COSLA*SINPHI
C      DXTIDE(2)=DXTIDE(2)-DR*SINLA*COSPHI+DN*SINLA*SINPHI  
C      DXTIDE(3)=DXTIDE(3)-DR*SINPHI      -DN*COSPHI
C       
      RETURN  
      END  
C  *************************************************************  
c        
c    COMPUTING OF THE UTC TIME
c      
      subroutine dutc(iyear,imonth,iday,dut)
c  Made from the subroutine sent by Dennis McCarthy 
      implicit none
      INTEGER iyear,imonth,iday
      INTEGER i,j
      DOUBLE PRECISION dut,an
      double precision utc(2,150)
C     DIFFERENCES UTC-ET
C        
C       
       data ((utc(i,j), i= 1, 2), j= 1, 88) /
     1 1955.5d0,31.59d0,1956.5d0,32.06d0,1957.5d0, 31.81999d0, 
     4 1958.5d0,32.68999d0,1959.5d0,33.04999d0,1960.5d0,33.15999d0,
     7 1961.5d0,33.59d0,1962.0d0, 34.032d0, 1962.5d0, 34.235d0,
     x 1963.0d0,34.441d0,1963.5d0, 34.644d0, 1964.0d0, 34.95d0, 
     3 1964.5d0,35.28599d0,1965.0d0,35.72499d0,1965.5d0,36.15999d0,  
     6 1966.0d0,36.498d0,1966.5d0,36.96799d0,1967.0d0,37.444d0,
     9 1967.5d0,37.91299d0,1968.0d0,38.38999d0,1968.25d0,38.526d0, 
     2 1968.5d0,38.75999d0,1968.75d0,39.0d0,1969.0d0,39.23799d0,
     5 1969.25d0,39.472d0,1969.5d0,39.707d0,1969.75d0,39.94599d0, 
     8 1970.0d0,40.185d0,1970.25d0,40.41999d0,1970.5d0,40.65399d0,
     1 1970.75d0,40.89199d0,1971.0d0,41.131d0,1971.085d0,41.21099d0, 
     4 1971.162d0,41.284d0,1971.247d0,41.36399d0,1971.329d0,41.442d0,
     7 1971.414d0,41.52199d0,1971.496d0,41.59999d0,1971.581d0,41.68d0,
     * 1971.666d0,41.761d0,1971.748d0,41.838d0,1971.833d0,41.91899d0,
     3 1971.915d0,41.99599d0,1971.99999d0,42.18399d0,1972.0d0,
     . 42.18399d0,1972.49999d0,42.18399d0,1972.5d0,43.18399d0,
     6 1972.99999d0,43.18399d0,1973.0d0,44.18399d0,1973.99999d0,
     9 44.18399d0,1974.0d0,45.18399d0,1974.99999d0,45.18399d0,
     .   1975.0d0,
     2 46.18399d0,1975.99999d0,46.18399d0,1976.0d0,47.18399d0,
     5 1976.99999d0,47.18399d0,1977.0d0,48.18399d0,1977.99999d0,
     8 48.18399d0,1978.0d0,49.18399d0,1978.99999d0,49.18399d0,   
     1 1979.0d0,50.18399d0,1979.99999d0,50.18399d0,1980.0d0,
     . 51.18399d0,1981.49999d0,51.18399d0,1981.5d0,52.18399d0,
     4 1982.49999d0,52.18399d0,1982.5d0,53.18399d0,1983.49999d0,
     7 53.18399d0,1983.5d0,54.18399d0,1985.49999d0,54.18399d0,
     x 1985.5d0,55.18399d0,1987.99999d0,55.18399d0, 1988.0d0,
     3 56.18399d0,1989.99999d0,56.18399d0,1990.0d0,57.18399d0, 
     6 1990.99999d0,57.18399d0,1991.0d0,58.18399d0,1992.49999d0,
     . 58.18399d0, 1992.5d0,59.18399d0,1993.49999d0,59.18399d0,
     9 1993.5d0,60.18399d0,1994.49999d0,60.18399d0,1995.99999d0,
     2 61.18399d0,1996.0d0,62.18399d0,  
     5 1997.5d0,63.18399d0,1999.0d0,64.18399d0,  
     6 2006.0d0,65.18399d0,2007.0d0,65.18399d0/  
!
!!     an=iyear+imonth/12d0+iday/365.25d0 
      an=iyear+(imonth-1)/12d0+iday/365.25d0 
      if (an.lt.1955.5) then 
       !! print *,'THE LEAP SECONDS ARE NOT GIVEN BEFORE 1955. '
       !! PRINT *,'THE CORRECTIONS COMPUTED FOR THE PRIOR YEARS ARE '
       !! print *,'COMPUTED WITH UTC-ET=31.59 2 SECOND.'
       dut=31.59 
       goto 510  
      endif 
      if (an.gt.2008.0d0) then 
       !! print *,'THE LEAPS SECONDS ARE NOT GIVEN AFTER 2007 '
       !! PRINT *,'THE CORRECTIONS COMPUTED FOR THE LATER YEARS ARE '
       !! PRINT *,'COMPUTED WITH UTC-ET=65.184 SECONDS.'
       dut=65.184
       goto 510  
      endif 
      do 5 i=1,87
        if (an.gt.utc(1,i)) then
          dut=utc(2,i)
          goto 510
        endif
   5  continue
 510  return 
      end     
C********************************************************************
      DOUBLE PRECISION function fjldy(iyr,imon,iday,rhour)
C compute julian day given year,month,day, and real hours
C Uses algorthim of Jean Meeus in "Astronomical Algorithms"
C
C Written by
C    John Gipson
C    April 8, 1996
C This uses the astronomical convention in which the year which
C historians call -1 BC is represented as 0.
C Should work for all years > -4716.  Checked on the following dates:
C
C   16. April   -1410, 0.00 H is FJLDY = 1206160.5D0    C
C   31. January -1100, 0.00 H is FJLDY = 1319312.5D0    C
C   24. January -0800, 0.00 H is FJLDY = 1428880.5D0    C
C   17. January -0500, 0.00 H is FJLDY = 1538448.5D0    C
C   10. January -0200, 0.00 H is FJLDY = 1648016.5D0    C
C   03. January   100, 0.00 H is FJLDY = 1757584.5D0    C
C   29. February  400, 0.00 H is FJLDY = 1867216.5D0    C
C   20. December  699, 0.00 H is FJLDY = 1976720.5D0    C
C   15. February 1000, 0.00 H is FJLDY = 2086352.5D0    C
C   08. February 1300, 0.00 H is FJLDY = 2195920.5D0    C
C   11. February 1600, 0.00 H is FJLDY = 2305488.5D0    C
C   06. February 1900, 0.00 H is FJLDY = 2415056.5D0    C
C   01. January  1988, 0.00 H is FJLDY = 2447161.5D0    C
C   01. February 1988, 0.00 H is FJLDY = 2447192.5D0    C
C   29. February 1988, 0.00 H is FJLDY = 2447220.5D0    C
C   01. March    1988, 0.00 H is FJLDY = 2447221.5D0    C
C   01. February 2200, 0.00 H is FJLDY = 2524624.5D0    C
C   27. January  2500, 0.00 H is FJLDY = 2634192.5D0    C
C   23. January  2800, 0.00 H is FJLDY = 2743760.5D0    C
C   22. December 3002, 0.00 H is FJLDY = 2817872.5D0    C


      implicit none
      INTEGER iyr,imon,iday
      DOUBLE PRECISION rhour

      INTEGER jyr,jmon

      IF(imon .GT. 2) then
        jmon=imon
        jyr=iyr
      else
        jmon=imon+12
        jyr=iyr-1
      endif

C Factor of .1 in months and years insures that roundoff errors will not be a
C problem. For example 30.6*5=183.  Some machines will calculate this as
C 182.9999..., and truncate to 182.
      fjldy=INT(365.25d0*float(jyr+4716)+.1)
     >    +INT(30.6*float(jmon+1)+0.1) +iday-1524.5+rhour/24.

C
C  Check to see if gregorian or julian date. These handle leap years differently, 
C  as well as have an offset. Switchover date was Oct 4 1582. Next day was Oct 15 1582.
C                                     
      IF(iyr*10000.+imon*100.+iday .GT. 15821015.)
     >    fjldy=fjldy+float(2-INT(jyr/100)+INT(INT(jyr/100)/4))

      return
      end
C*********************************************************************************************
C 
C 
      SUBROUTINE SPROD(X,Y,SCAL,R1,R2) 
C 
C  COMPUTATION OF THE SCALAR-PRODUCT OF TWO VECTORS AND THEIR NORMS 
C 
C  INPUT :  X(I),I=1,2,3: COMPONENTS OF VECTOR X 
C           Y(I),I=1,2,3: COMPONENTS OF VECTOR Y 
C  OUTPUT :  SCAL: SCALAR PRODUCT OF X AND Y 
C            R1,R2  : LENGTHS OF THE TWO VECTORS X AND Y 
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(3),Y(3)
      R1=DSQRT(X(1)**2+X(2)**2+X(3)**2) 
      R2=DSQRT(Y(1)**2+Y(2)**2+Y(3)**2) 
      SCAL=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3) 
      RETURN 
      END 
C 
C-------------------------------------------------------------------------
C
      SUBROUTINE ST1L1(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)
C  
C THIS SUBROUTINE GIVES THE CORRECTIONS INDUCED BY THE LATITUDE DEPENDENCE  
C GIVEN BY L^(1) IN MAHTEWS ET AL (1991)  
C  
C       INPUT : XSTA,XSUN,XMON,FAC3SUN,FAC3MON  
C      OUTPUT : XCORSTA  
C  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C compute the norm of a vector
      DOUBLE PRECISION norm8
      DIMENSION XSTA(3),XSUN(3),XMON(3),XCORSTA(3)  
      DOUBLE PRECISION L1,L1D,L1SD
      DATA L1D/0.0012d0/,L1SD/0.0024d0/
      RSTA=norm8(xsta)
      SINPHI=XSTA(3)/RSTA  
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA  
      SINLA=XSTA(2)/COSPHI/RSTA  
      COSLA=XSTA(1)/COSPHI/RSTA  
      RMON=norm8(XMON)
      rsun=norm8(xsun)
C
C FOR THE DIURNAL BAND  
C  
      L1=L1D  
      DNSUN=-L1*SINPHI**2*FAC2SUN*XSUN(3)*(XSUN(1)*COSLA+XSUN(2)*SINLA)
     &            /RSUN**2
      DNMON=-L1*SINPHI**2*FAC2MON*XMON(3)*(XMON(1)*COSLA+XMON(2)*SINLA)
     &            /RMON**2
      DESUN=L1*SINPHI*(COSPHI**2-SINPHI**2)*FAC2SUN*XSUN(3)*
     1 (XSUN(1)*SINLA-XSUN(2)*COSLA)/RSUN**2
      DEMON=L1*SINPHI*(COSPHI**2-SINPHI**2)*FAC2MON*XMON(3)*
     1 (XMON(1)*SINLA-XMON(2)*COSLA)/RMON**2
      DE=3.*(DESUN+DEMON)  
      DN=3.*(DNSUN+DNMON)  
      XCORSTA(1)=-DE*SINLA-DN*SINPHI*COSLA  
      XCORSTA(2)=DE*COSLA-DN*SINPHI*SINLA  
      XCORSTA(3)=DN*COSPHI  
C   
C FOR THE SEMI-DIURNAL BAND  
C  
      L1=L1SD  
      COSTWOLA=COSLA**2-SINLA**2  
      SINTWOLA=2.*COSLA*SINLA  
      DNSUN=-L1/2.*SINPHI*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*
     1 COSTWOLA+2.*XSUN(1)*XSUN(2)*SINTWOLA)/RSUN**2
      DNMON=-L1/2.*SINPHI*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)*
     1 COSTWOLA+2.*XMON(1)*XMON(2)*SINTWOLA)/RMON**2
      DESUN=-L1/2.*SINPHI**2*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*
     1 SINTWOLA-2.*XSUN(1)*XSUN(2)*COSTWOLA)/RSUN**2
      DEMON=-L1/2.*SINPHI**2*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)*
     1 SINTWOLA-2.*XMON(1)*XMON(2)*COSTWOLA)/RMON**2
      DE=3.*(DESUN+DEMON)  
      DN=3.*(DNSUN+DNMON)  
      XCORSTA(1)=XCORSTA(1)-DE*SINLA-DN*SINPHI*COSLA  
      XCORSTA(2)=XCORSTA(2)+DE*COSLA-DN*SINPHI*SINLA  
      XCORSTA(3)=XCORSTA(3)+DN*COSPHI  
      RETURN  
      END  

C*************************************************************************
C     Last change:  VD   17 May 00   1:20 pm
C  THESE ARE THE SUBROUTINES FOR THE STEP2 OF THE TIDAL CORRECTIONS. 
C  THEY ARE CALLED TO ACCOUNT FOR THE FREQUENCY DEPENDENCE  
C  OF THE LOVE NUMBERS. 
C 
      SUBROUTINE STEP2DIU(XSTA,FHR,T,XCORSTA)  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XSTA(3),XCORSTA(3),DATDI(9,31)
      DOUBLE PRECISION deg2rad
      DATA deg2rad/0.0174532925d0/
      DATA ((DATDI(i,j),i=1,9),j=1,31)/  
     * -3., 0., 2., 0., 0.,-0.01, 0.0 , 0.0 , 0.0,   
     * -3., 2., 0., 0., 0.,-0.01, 0.0 , 0.0 , 0.0,   
     * -2., 0., 1.,-1., 0.,-0.02, 0.0 , 0.0 , 0.0,   
     * -2., 0., 1., 0., 0.,-0.08, 0.0 ,-0.01, 0.01,
     * -2., 2.,-1., 0., 0.,-0.02, 0.0 , 0.0 , 0.0,
     * -1., 0., 0.,-1., 0.,-0.10, 0.0 , 0.0 , 0.0,
     * -1., 0., 0., 0., 0.,-0.51, 0.0 ,-0.02, 0.03,
     * -1., 2., 0., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  0.,-2., 1., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  0., 0.,-1., 0., 0., 0.02, 0.0 , 0.0 , 0.0,
     *  0., 0., 1., 0., 0., 0.06, 0.0 , 0.0 , 0.0,
     *  0., 0., 1., 1., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  0., 2.,-1., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  1.,-3., 0., 0., 1.,-0.06, 0.0 , 0.0 , 0.0,
     *  1.,-2., 0.,-1., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  1.,-2., 0., 0., 0.,-1.23,-0.07, 0.06, 0.01,
     *  1.,-1., 0., 0.,-1., 0.02, 0.0 , 0.0 , 0.0,
     *  1.,-1., 0., 0., 1., 0.04, 0.0 , 0.0 , 0.0,
     *  1., 0., 0.,-1., 0.,-0.22, 0.01, 0.01, 0.0,
     *  1., 0., 0., 0., 0.,12.00,-0.80,-0.67,-0.03,
     *  1., 0., 0., 1., 0., 1.73,-0.12,-0.10, 0.0,
     *  1., 0., 0., 2., 0.,-0.04, 0.0 , 0.0 , 0.0, 
     *  1., 1., 0., 0.,-1.,-0.50,-0.01, 0.03, 0.0,
     *  1., 1., 0., 0., 1., 0.01, 0.0 , 0.0 , 0.0,
     *  0., 1., 0., 1.,-1.,-0.01, 0.0 , 0.0 , 0.0,
     *  1., 2.,-2., 0., 0.,-0.01, 0.0 , 0.0 , 0.0,
     *  1., 2., 0., 0., 0.,-0.11, 0.01, 0.01, 0.0,
     *  2.,-2., 1., 0., 0.,-0.01, 0.0 , 0.0 , 0.0,
     *  2., 0.,-1., 0., 0.,-0.02, 0.0 , 0.0 , 0.0,
     *  3., 0., 0., 0., 0., 0.0 , 0.0 , 0.0 , 0.0,
     *  3., 0., 0., 1., 0., 0.0 , 0.0 , 0.0 , 0.0/
      S=218.31664563D0+481267.88194D0*T-0.0014663889D0*T**2 
     1 +0.00000185139D0*T**3 
      TAU=fhr*15.D0+280.4606184D0+36000.7700536D0*T+0.00038793D0*T**2 
     1 -0.0000000258D0*T**3-S 
      PR=1.396971278*T+0.000308889*T**2+0.000000021*T**3 
     1 +0.000000007*T**4 
      S=S+PR 
      H=280.46645D0+36000.7697489D0*T+0.00030322222D0*T**2 
     1 +0.000000020*T**3-0.00000000654*T**4  
      P=83.35324312D0+4069.01363525D0*T-0.01032172222D0*T**2
     1 -0.0000124991D0*T**3+0.00000005263D0*T**4  
      ZNS=234.95544499D0 +1934.13626197D0*T-0.00207561111D0*T**2
     1 -0.00000213944D0*T**3+0.00000001650D0*T**4  
      PS=282.93734098D0+1.71945766667D0*T+0.00045688889D0*T**2 
     1 -0.00000001778D0*T**3-0.00000000334D0*T**4
C Reduce angles to between 0 and 360.
      s=  dmod(s,360.d0)
      tau=dmod(tau,360.d0)
      h=  dmod(h,360.d0)
      p=  dmod(p,360.d0)
      zns=dmod(zns,360.d0)
      ps=dmod(ps,360.d0)
c      WRITE(2,'(6f10.3)') tau,s,h,p,zns,ps

      RSTA=DSQRT(XSTA(1)**2+XSTA(2)**2+XSTA(3)**2)  
      SINPHI=XSTA(3)/RSTA  
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA  

      COSLA=XSTA(1)/COSPHI/RSTA
      SINLA=XSTA(2)/COSPHI/RSTA
      ZLA = DATAN2(XSTA(2),XSTA(1))
      DO 99 I=1,3  
      XCORSTA(I)=0.
99    continue
      DO 98 J=1,31
      THETAF=(TAU+DATDI(1,J)*S+DATDI(2,J)*H+DATDI(3,J)*P+
     1 DATDI(4,J)*ZNS+DATDI(5,J)*PS)*deg2rad
      DR=DATDI(6,J)*2.*SINPHI*COSPHI*SIN(THETAF+ZLA)+
     1 DATDI(7,J)*2.*SINPHI*COSPHI*COS(THETAF+ZLA)
      DN=DATDI(8,J)*(COSPHI**2-SINPHI**2)*SIN(THETAF+ZLA)+
     1 DATDI(9,J)*(COSPHI**2-SINPHI**2)*COS(THETAF+ZLA)
      DE=DATDI(8,J)*SINPHI*COS(THETAF+ZLA)+
     1 DATDI(9,J)*SINPHI*SIN(THETAF+ZLA)
      XCORSTA(1)=XCORSTA(1)+DR*COSLA*COSPHI-DE*SINLA  
     1 -DN*SINPHI*COSLA  
      XCORSTA(2)=XCORSTA(2)+DR*SINLA*COSPHI+DE*COSLA  
     1 -DN*SINPHI*SINLA  
      XCORSTA(3)=XCORSTA(3)+DR*SINPHI+DN*COSPHI  
98    CONTINUE   

      DO 97 I=1,3
      XCORSTA(I)=XCORSTA(I)/1000.
97    CONTINUE  
      RETURN  
      END  
C  
C  *************************************************************  
C
      SUBROUTINE STEP2LON(XSTA,FHR,T,XCORSTA)  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION deg2rad
      DOUBLE PRECISION XSTA(3),XCORSTA(3),DATDI(9,5)
      DATA deg2rad/0.0174532925d0/ 
      DATA ((DATDI(i,j),i=1,9),j=1,5)/  
     *   0, 0, 0, 1, 0,   0.47, 0.23, 0.16, 0.07,
     *   0, 2, 0, 0, 0,  -0.20,-0.12,-0.11,-0.05,
     *   1, 0,-1, 0, 0,  -0.11,-0.08,-0.09,-0.04,
     *   2, 0, 0, 0, 0,  -0.13,-0.11,-0.15,-0.07,
     *   2, 0, 0, 1, 0,  -0.05,-0.05,-0.06,-0.03/
C
      S=218.31664563D0+481267.88194D0*T-0.0014663889D0*T**2 
     1 +0.00000185139D0*T**3 
      PR=1.396971278*T+0.000308889*T**2+0.000000021*T**3 
     1 +0.000000007*T**4 
      S=S+PR 
      H=280.46645D0+36000.7697489D0*T+0.00030322222D0*T**2 
     1 +0.000000020*T**3-0.00000000654*T**4  
      P=83.35324312D0+4069.01363525D0*T-0.01032172222D0*T**2 
     1 -0.0000124991D0*T**3+0.00000005263D0*T**4  
      ZNS=234.95544499D0 +1934.13626197D0*T-0.00207561111D0*T**2
     1 -0.00000213944D0*T**3+0.00000001650D0*T**4  
      PS=282.93734098D0+1.71945766667D0*T+0.00045688889D0*T**2 
     1 -0.00000001778D0*T**3-0.00000000334D0*T**4  
      RSTA=DSQRT(XSTA(1)**2+XSTA(2)**2+XSTA(3)**2)  
      SINPHI=XSTA(3)/RSTA  
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA  
      COSLA=XSTA(1)/COSPHI/RSTA
      SINLA=XSTA(2)/COSPHI/RSTA
C reduce angles to between 0 and 360.
      s=  dmod(s,360.d0)
      tau=dmod(tau,360.d0)
      h=  dmod(h,360.d0)
      p=  dmod(p,360.d0)
      zns=dmod(zns,360.d0)
      ps=dmod(ps,360.d0)

      dr_tot=0.
      dn_tot=0.
      DO 99 I=1,3  
      XCORSTA(I)=0.
99    continue
      DO 98 J=1,5
      THETAF=(DATDI(1,J)*S+DATDI(2,J)*H+DATDI(3,J)*P+
     1 DATDI(4,J)*ZNS+DATDI(5,J)*PS)*DEG2RAD
      DR=DATDI(6,J)*(3.*SINPHI**2-1.)/2.*COS(THETAF)+
     1     DATDI(8,J)*(3.*SINPHI**2-1.)/2.*SIN(THETAF)
      DN=DATDI(7,J)*(COSPHI*SINPHI*2.)*COS(THETAF)+
     1     DATDI(9,J)*(COSPHI*SINPHI*2.)*SIN(THETAF)
      DE=0. 
      dr_tot=dr_tot+dr
      dn_tot=dn_tot+dn
      XCORSTA(1)=XCORSTA(1)+DR*COSLA*COSPHI-DE*SINLA  
     1 -DN*SINPHI*COSLA  
      XCORSTA(2)=XCORSTA(2)+DR*SINLA*COSPHI+DE*COSLA  
     1 -DN*SINPHI*SINLA  
      XCORSTA(3)=XCORSTA(3)+DR*SINPHI+DN*COSPHI  
98    CONTINUE   

      DO 97 I=1,3
      XCORSTA(I)=XCORSTA(I)/1000.
97    CONTINUE  
      RETURN  
      END  
C**************************************************************************************************
C-------------------------------------------------------------------------
C 
      SUBROUTINE ST1IDIU(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)
C  
C THIS SUBROUTINE GIVES THE OUT-OF-PHASE CORRECTIONS INDUCED BY 
C MANTLE INELASTICITY IN THE DIURNAL BAND 
C  
C       INPUT : XSTA,XSUN,XMON,FAC2SUN,FAC2MON  
C      OUTPUT : XCORSTA  
C  

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION norm8
      DIMENSION XSTA(3),XSUN(3),XMON(3),XCORSTA(3)  
      DATA DHI/-0.0025/,DLI/-0.0007/  
      RSTA=NORM8(XSTA)
      SINPHI=XSTA(3)/RSTA  
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA
      COS2PHI=COSPHI**2-SINPHI**2
      SINLA=XSTA(2)/COSPHI/RSTA  
      COSLA=XSTA(1)/COSPHI/RSTA  
      RMON=NORM8(XMON)
      RSUN=NORM8(XSUN)
      DRSUN=-3.*DHI*SINPHI*COSPHI*FAC2SUN*XSUN(3)*(XSUN(1)*
     1            SINLA-XSUN(2)*COSLA)/RSUN**2
      DRMON=-3.*DHI*SINPHI*COSPHI*FAC2MON*XMON(3)*(XMON(1)*
     1            SINLA-XMON(2)*COSLA)/RMON**2
      DNSUN=-3.*DLI*COS2PHI*FAC2SUN*XSUN(3)*(XSUN(1)*SINLA-
     1            XSUN(2)*COSLA)/RSUN**2
      DNMON=-3.*DLI*COS2PHI*FAC2MON*XMON(3)*(XMON(1)*SINLA-
     1            XMON(2)*COSLA)/RMON**2
      DESUN=-3.*DLI*SINPHI*FAC2SUN*XSUN(3)*
     1 (XSUN(1)*COSLA+XSUN(2)*SINLA)/RSUN**2
      DEMON=-3.*DLI*SINPHI*FAC2MON*XMON(3)*
     1 (XMON(1)*COSLA+XMON(2)*SINLA)/RMON**2
      DR=DRSUN+DRMON 
      DN=DNSUN+DNMON  
      DE=DESUN+DEMON 
      XCORSTA(1)=DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA  
      XCORSTA(2)=DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA  
      XCORSTA(3)=DR*SINPHI+DN*COSPHI  
      RETURN  
      END  
C
C-------------------------------------------------------------------------
      SUBROUTINE ST1ISEM(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)
C  
C THIS SUBROUTINE GIVES THE OUT-OF-PHASE CORRECTIONS INDUCED BY 
C MANTLE INELASTICITY IN THE DIURNAL BAND 
C  
C       INPUT : XSTA,XSUN,XMON,FAC2SUN,FAC2MON  
C      OUTPUT : XCORSTA  
C  

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION norm8
      DIMENSION XSTA(3),XSUN(3),XMON(3),XCORSTA(3)  
      DATA DHI/-0.0022/,DLI/-0.0007/  
      RSTA=NORM8(XSTA)
      SINPHI=XSTA(3)/RSTA  
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA
      SINLA=XSTA(2)/COSPHI/RSTA  
      COSLA=XSTA(1)/COSPHI/RSTA  
      COSTWOLA=COSLA**2-SINLA**2  
      SINTWOLA=2.*COSLA*SINLA  
      RMON=NORM8(XMON)
      RSUN=NORM8(XSUN)
      DRSUN=-3./4.*DHI*COSPHI**2*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*  
     1 SINTWOLA-2.*XSUN(1)*XSUN(2)*COSTWOLA)/RSUN**2
      DRMON=-3./4.*DHI*COSPHI**2*FAC2MON*((XMON(1)**2-XMON(2)**2)*  
     1 SINTWOLA-2.*XMON(1)*XMON(2)*COSTWOLA)/RMON**2
      DNSUN=3./2.*DLI*SINPHI*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*  
     1 SINTWOLA-2.*XSUN(1)*XSUN(2)*COSTWOLA)/RSUN**2
      DNMON=3./2.*DLI*SINPHI*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)*  
     1 SINTWOLA-2.*XMON(1)*XMON(2)*COSTWOLA)/RMON**2
      DESUN=-3./2.*DLI*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*  
     1 COSTWOLA+2.*XSUN(1)*XSUN(2)*SINTWOLA)/RSUN**2
      DEMON=-3./2.*DLI*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)*  
     1 COSTWOLA+2.*XMON(1)*XMON(2)*SINTWOLA)/RMON**2
      DR=DRSUN+DRMON 
      DN=DNSUN+DNMON  
      DE=DESUN+DEMON 
      XCORSTA(1)=DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA  
      XCORSTA(2)=DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA  
      XCORSTA(3)=DR*SINPHI+DN*COSPHI  
      RETURN  
      END  
C
C*********************************************************************
      DOUBLE PRECISION function norm8(a)
      DOUBLE PRECISION a(3)
      norm8=dSQRT(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      return
      end
C
C*********************************************************************
      subroutine zero_vec8(v)
      DOUBLE PRECISION v(3)
      v(1)=0.
      v(2)=0.
      v(3)=0.
      return
      end

