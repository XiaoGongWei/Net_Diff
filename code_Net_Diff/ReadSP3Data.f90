! ===============  ReadSP3Data  =============
!
! PURPOSE:
!     Read the GNS sp3 data for n epochs.
!
! INPUTS:
!      n              read n epochs of sp3 data
!
! OUTPUT:
!      SP3Data   the added sp3 data, in SP3Data module
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!===============  End of Header  ==============

subroutine ReadSP3Data(n)
use MOD_FileID
use MOD_SP3Data
use MOD_SP3Head
implicit none
    integer :: error
    integer ::  n,i,j,prn
    character(1) :: symbol, system
    character(80) :: line
    integer ::  y, mon,d, h, min, GPSweek
    real(8) :: s, GPSsec, Coor(3), Clk

    
    ! Read n epochs of SP3 Data
    Read_Epoch_Loop: do i=1,n
        read(unit=SP3ID,fmt="(A)",iostat=error) line
        if (error /=0) exit
        
        if (line(1:1)=="*") then
            ! ************Shift the data***********
            SP3Data%GPSweek=eoshift(SP3Data%GPSweek,shift=1,dim=1)
            SP3Data%GPSsec=eoshift(SP3Data%GPSsec,shift=1,dim=1)
            do j=1,GNum+RNum+CNum+NumE+JNum+INum
                SP3Data%Eph(j)%Coor=eoshift(SP3Data%Eph(j)%Coor,shift=1,dim=2,boundary=9999.0d0)
                SP3Data%Eph(j)%Clk=eoshift(SP3Data%Eph(j)%Clk,shift=1,boundary=9999.0d0)
            end do
            ! *************************************
            
             read(line,"(5X,5I3,F11.8)") y, mon, d, h, min, s
             call UTC2GPST(y, mon, d, h, min, s, SP3Data%GPSweek(10), SP3Data%GPSsec(10))
             j=1
             Read_PRN_Loop: do while (j<=SP3Head%PRNS)
                  read(unit=SP3ID,fmt="(A1,A1,I2,4F14.6)") symbol,system,prn, Coor, Clk
                   j=j+1
                  if (symbol /="P") cycle
                  if (system=="R") then
                      if (PRN>RNum) cycle
                      prn=prn+GNum
                  elseif (system=="C") then
                      if (PRN>CNum) cycle
                      prn=prn+GNum+RNum
                  elseif (system=="E") then
                      if (PRN>NumE) cycle
                      prn=prn+GNum+RNum+CNum
                  elseif (system=="J") then
                      if (PRN>JNum) cycle
                      prn=prn+GNum+RNum+CNum+NumE
                  elseif (system=="I") then
                      if (PRN>INum) cycle
                      prn=prn+GNum+RNum+CNum+NumE+JNum
                  elseif (system=="G") then
                      if (PRN>GNum) cycle
                      prn=prn
                  else
                      cycle
                  end if
                  !read(unit=SP3ID,fmt="(4F14.6)") Coor, SP3Data%Eph(prn)%Clk(10)
                  SP3Data%Eph(prn)%Coor(1:3,10)=Coor ! MATMUL(RT2C, Coor) ! transfer to TRS(terrestial reference system)
                  SP3Data%Eph(prn)%Clk(10)=Clk
             end do Read_PRN_Loop      
        end if   ! if (line(1:1)=="*")
        
    end do Read_Epoch_Loop

end subroutine