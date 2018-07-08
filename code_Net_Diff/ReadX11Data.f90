
! ==========Read  X11 Data(GPS)========
! Inputs:
!      Obsweek         Obs time of GPS week
!      Obssec            Obs time of GPS second
!
! Output:
!      ObsData          Including C1, P1, P2, P3, L1, L2, L3 data
!
! Written by: Yize Zhang
! ==============End of Header=============

subroutine ReadX11Data(Obsweek, Obssec, ObsData,sta)
use MOD_ObsData
use MOD_FileID
use MOD_constant
use MOD_VAR
implicit none
    type(type_ObsData) :: ObsData
    integer(1) :: Flag(6), lines
    integer :: Obsweek
    real(8) :: Obssec
    integer :: i, sta, j, PRN
    character*600 :: line
    integer :: year, mon, day, hour, min, week
    real(8) :: sow, MJD
    real(8) :: val(36), meanP

    if ( allocated(ObsData%PRN) ==.false.) then
    allocate(ObsData%PRN(32), ObsData%System(32), ObsData%C1(32),ObsData%P1(32),  &
                        ObsData%P2(32),ObsData%P3(32),ObsData%L1(32), ObsData%L2(32), ObsData%L3(32), &
                        ObsData%LLI1(32), ObsData%LLI2(32), ObsData%LLI3(32) )
    end if
    ObsData%C1=0.0d0
    ObsData%P1=0.0d0
    ObsData%P2=0.0d0
    ObsData%P3=0.0d0
    ObsData%L1=0.0d0
    ObsData%L2=0.0d0
    ObsData%L3=0.0d0
    ObsData%LLI1=0
    ObsData%LLI2=0
    ObsData%LLI3=0

    ObsData%PRNS=0
    ObsData%PRN=0
    ObsData%Flag=0
    ObsData%Clk_Bias=0.d0
    i=0
    ObsData%Flag=1   ! First set the flag false, in case that there is no data in this epoch
    do while(.true.)
            100 read(unit=ObsID(sta),fmt="(A)",end=200) line
            if ((len_trim(line)<324) .or. (len_trim(line)>325)) cycle
            if (index(line,'********')/=0) cycle

            read(line, *) week, sow,val(1:2), PRN
            if (week<1000) then
                week=week+1356
            end if
            if ((week-Obsweek)*604800.d0+sow-Obssec > 0.01d0) then ! If no data in this epoch
                backspace(ObsID(sta))
                exit
            elseif ((week-Obsweek)*604800.d0+sow-Obssec < -0.01d0) then
                cycle
            end if

            do j=1,i   ! delete repeat data
                if (PRN==ObsData%PRN(j)) goto 100
            end do
            ObsData%week=week
            ObsData%sow=Obssec
            ObsData%Flag=0
            call GPST2UTC(week, sow, MJD,ObsData%year,ObsData%mon,ObsData%day,ObsData%hour,ObsData%min,ObsData%sec)
            ObsData%PRNS=ObsData%PRNS+1
            i=i+1
            ObsData%System(i)='G'
            read(line,*)  val(1:4), ObsData%PRN(i), val(1), meanP, val(1:36)
            if (val(1)==0.d0) ObsData%P1(i)=MeanP+val(3)
            if (val(19)==0.d0) ObsData%P2(i)=MeanP+val(21)
            ObsData%L1(i)=val(5)
            ObsData%L2(i)=val(23)
    end do

    200 return ! When reach the end of the file
    
end subroutine