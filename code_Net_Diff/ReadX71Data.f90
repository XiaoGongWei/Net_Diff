
! ==========Read  X71 Data========
! Inputs:
!      Obsweek         Obs time of GPS week
!      Obssec            Obs time of GPS second
!
! Output:
!      ObsData          Including C1, P1, P2, P3, L1, L2, L3 data
!
! Written by: Yize Zhang
! ==============End of Header=============

subroutine ReadX71Data(Obsweek, Obssec, ObsData,sta)
use MOD_ObsData
use MOD_FileID
use MOD_constant
use MOD_VAR
implicit none
    type(type_ObsData) :: ObsData
    integer(1) :: Flag(15), lines
    integer :: Obsweek
    real(8) :: Obssec
    integer :: i, sta
    character*600 :: line
    integer :: year, mon, day, hour, min, week
    real(8) :: sow, MJD
    real(8) :: val(12), STD,rela,tide,eccen(3),com(3)
    real(8) :: temp(12)

    if ( allocated(ObsData%PRN) ==.false.) then
    allocate(ObsData%PRN(20), ObsData%System(20), ObsData%C1(20),ObsData%P1(20),  &
                        ObsData%P2(20),ObsData%P3(20),ObsData%L1(20), ObsData%L2(20), ObsData%L3(20), &
                        ObsData%LLI1(20), ObsData%LLI2(20), ObsData%LLI3(20) )
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
    ObsData%Flag=0
    ObsData%Clk_Bias=0.d0
    i=0
    ObsData%Flag=1   ! First set the flag false, in case that there is no data in this epoch
    do while(.true.)
            read(unit=ObsID(sta),fmt="(A)",end=200) line
            if (len_trim(line)/=559) cycle
            if (index(line,'********')/=0) cycle

            read(line, "(I4,F7.0)") week, sow
            week=week+1356
            if ((week-Obsweek)*604800.d0+sow-Obssec > 0.01d0) then ! If no data in this epoch
                backspace(ObsID(sta))
                exit
            elseif ((week-Obsweek)*604800.d0+sow-Obssec < -0.01d0) then
                cycle
            end if

            ObsData%week=week
            ObsData%sow=Obssec
            ObsData%Flag=0
            call GPST2UTC(week, sow, MJD,ObsData%year,ObsData%mon,ObsData%day,ObsData%hour,ObsData%min,ObsData%sec)
            ObsData%PRNS=ObsData%PRNS+1
            i=i+1
            ObsData%System(i)='C'
            read(line,'(18X,I2,12F14.3,48X,3F15.4,71X,12I2,10X,3I2,71X,3F9.3,9X,6F9.3)') ObsData%PRN(i), val(1:12), &
                      ObsData%L1(i),ObsData%L2(i),ObsData%L3(i),flag(1:12),flag(13:15),STD,rela,tide,eccen(1:3),com(1:3)
             ! std=0.d0
             ! eccen(1:3)=0.d0
            ObsData%STD(i)=std
            ObsData%rela(i)=rela
            if (IorQ==0) then
                ObsData%P1(i)=val(1)-std-rela-tide-eccen(1)
                ObsData%P2(i)=val(5)-std-rela-tide-eccen(2)
                ObsData%P3(i)=val(9)-std-rela-tide-eccen(3)
            else
                ObsData%P1(i)=val(3)-std-rela-tide-eccen(1)
                ObsData%P2(i)=val(7)-std-rela-tide-eccen(2)
                ObsData%P3(i)=val(11)-std-rela-tide-eccen(3)
            end if
            ObsData%L1(i)=ObsData%L1(i)+(-std-rela-tide-eccen(1))*10.23d6*152.6d0/c
            ObsData%L2(i)=ObsData%L2(i)+(-std-rela-tide-eccen(2))*10.23d6*118.0d0/c
            ObsData%L3(i)=ObsData%L3(i)+(-std-rela-tide-eccen(3))*10.23d6*124.0d0/c
            
            if (Orbit=="SP3")  then  ! 精密星历时改到质心上，不用，因为用了ant文件
!                ObsData%P1(i)=ObsData%P1(i)-com(1)  ! Process_Corr时与pcor保持一致
!                ObsData%P2(i)=ObsData%P2(i)-com(2)
!                ObsData%P3(i)=ObsData%P3(i)-com(3)
!                ObsData%L1(i)=ObsData%L1(i)-com(1)*10.23d6*152.6d0/c
!                ObsData%L2(i)=ObsData%L2(i)-com(2)*10.23d6*118.0d0/c
!                ObsData%L3(i)=ObsData%L3(i)-com(3)*10.23d6*124.0d0/c
            else  ! 广播星历时改到B3相心上
                ObsData%P1(i)=ObsData%P1(i)-com(1)+com(3)
                ObsData%P2(i)=ObsData%P2(i)-com(2)+com(3)
                ObsData%P3(i)=ObsData%P3(i)
                ObsData%L1(i)=ObsData%L1(i)+(com(3)-com(1))*10.23d6*152.6d0/c
                ObsData%L2(i)=ObsData%L2(i)+(com(3)-com(2))*10.23d6*118.0d0/c
                ObsData%L3(i)=ObsData%L3(i)
            end if
            
            if (IorQ==0) then
                if (flag(2)==1) ObsData%P1(i)=0.d0
                if (flag(6)==1) ObsData%P2(i)=0.d0
                if (flag(10)==1) ObsData%P3(i)=0.d0
            else
                if (flag(4)==1) ObsData%P1(i)=0.d0
                if (flag(8)==1) ObsData%P2(i)=0.d0
                if (flag(12)==1) ObsData%P3(i)=0.d0
            end if
            if (flag(13)==1) ObsData%L1(i)=0.d0
            if (flag(14)==1) ObsData%L2(i)=0.d0
            if (flag(15)==1) ObsData%L3(i)=0.d0
    end do

!    do while(.true.)
!            read(unit=ObsID(sta),fmt="(A)",end=200) line
!            if (index(line,'********') /= 0)  cycle
!        
!            read(line,*) week, sow
!            week=week+1356
!            if ((week-Obsweek)*604800.d0+sow-Obssec > 0.01d0) then ! If no data in this epoch
!                backspace(ObsID(sta))
!                exit
!            elseif ((week-Obsweek)*604800.d0+sow-Obssec < -0.01d0) then
!                cycle
!            end if
!
!            ObsData%week=week
!            ObsData%sow=Obssec
!            call GPST2UTC(week, sow, MJD,ObsData%year,ObsData%mon,ObsData%day,ObsData%hour,ObsData%min,ObsData%sec)
!            ObsData%PRNS=ObsData%PRNS+1
!            i=i+1
!            ObsData%System(i)='C'
!            read(line,*)  week, sow, temp(1:2), ObsData%PRN(i), val(1:12),temp(1:12), &
!                      ObsData%L1(i),ObsData%L2(i),ObsData%L3(i)
!            ObsData%P1(i)=val(2)
!            ObsData%P2(i)=val(6)
!            ObsData%P3(i)=val(10)
!    end do

    200 return ! When reach the end of the file
    
end subroutine