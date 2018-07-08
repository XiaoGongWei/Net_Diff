! ==========Read Multi Orbit==========
! Read the Multi Orbit of COMASS

! Input:
!         fid      FILE ID
!        PRN     PRN
!         n        epoches
!       
! WRITTEN BY: Yize Zhang
!==========End of Header==========

subroutine ReadMultiOrb(n,week1, sow1)
use MOD_SP3Data
implicit none
    integer :: error
    integer ::  n,i,j,PRN
    character(80) :: line
    integer :: week,week0,week1
    real(8) :: sow,sow0, sow1, Coor(3)

    week=week1
    sow=sow1
    ! Read n epochs of SP3 Data
    Read_Epoch_Loop: do i=1,n
        ! ************Shift the data***********
        SP3Data%GPSweek=eoshift(SP3Data%GPSweek,shift=1,dim=1)
        SP3Data%GPSsec=eoshift(SP3Data%GPSsec,shift=1,dim=1)
        do j=1,91
            SP3Data%Eph(j)%Coor=eoshift(SP3Data%Eph(j)%Coor,shift=1,dim=2,boundary=9999.0d0)
        end do
        ! *************************************

        SP3Data%GPSweek(10)=week
        SP3Data%GPSsec(10)=sow
        Read_PRN_Loop: do PRN=1,14
            100 read(unit=100+PRN,fmt="(A)",iostat=error) line
            if (error /=0) cycle

            read(line,"(I4,F7.0)") week0, sow0
            week0=week0+1356
            if ( (week0-week)*604800.d0+sow0-sow>0.d0 ) then
                backspace(100+PRN)
                cycle
            elseif ( (week0-week)*604800.d0+sow0-sow<0.d0 ) then
                goto 100
            end if
            
             read(unit=line,fmt="(11X,3F19.3)") Coor
             SP3Data%Eph(PRN+56)%Coor(1:3,10)=Coor/1000.d0

        end do Read_PRN_Loop 
     end do Read_Epoch_Loop
end subroutine