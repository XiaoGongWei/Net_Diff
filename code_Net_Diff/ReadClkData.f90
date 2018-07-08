!  ============ReadClkData==============
! Read the GNSS(GPS+GLONASS) clock data for one epoch.
!  Function:
!        Read 5-min interval GNSS satellite clock data
!
! Written by: Yize Zhang
!  ==========End of header==========

subroutine ReadClkData
use MOD_FileID
use MOD_ClkData
implicit none
    integer :: error
    integer ::  j
    character(len=2) :: clk_type
    character(len=4) :: name
    integer :: PRN
    integer :: year, mon, day, hour, min, min0
    real(8) :: sec, sec0
    real(8) :: clk
    

    read(unit=ClkID,fmt="(8X,I4,4I3,F10.6)",iostat=error)  year, mon, day, hour, min0, sec0
    if (error /=0) return
    
    ! ****Shift the clock data********
    do j=1,GNum+RNum+CNum+NumE+JNum+INum  ! Satellite clock velocity
        if ( (ClkData%AS(j)%Clk(1)/=9999.d0) .and. (ClkData%AS(j)%Clk(2)/=9999.d0)) then
            ClkData%AS(j)%ClkVel=(ClkData%AS(j)%Clk(2)-ClkData%AS(j)%Clk(1))/ &
            ((ClkData%GPSweek(2)-ClkData%GPSweek(1))*604800.d0+ClkData%GPSsec(2)-ClkData%GPSsec(1))
        end if
    end do
    ClkData%GPSweek(1)=ClkData%GPSweek(2)
    ClkData%GPSsec(1)=ClkData%GPSsec(2)
    do j=1,GNum+RNum+CNum+NumE+JNum+INum  ! AS
        ! First set the second clk 9999.0
        ClkData%AS(j)%Clk=eoshift(ClkData%AS(j)%Clk,shift=1,dim=1,boundary=9999.d0)
    end do
    ! *****************************
    
    call UTC2GPST(year, mon, day, hour, min0, sec0, ClkData%GPSweek(2), ClkData%GPSsec(2))
    backspace(ClkID)
    do while(.true.)
        read(unit=ClkID,fmt="(A2,1X,A4,1X,I4,4I3,F10.6,6X,E19.12)", iostat=error)  clk_type, name, year, mon, day, hour, min, sec, clk
        if (error /=0) exit
        if ((min /= min0) .or. (abs(sec-sec0)>0.1d0)) then
            backspace(ClkID)
            exit
        end if
        
        if (clk_type=="AS") then
            read(name(2:3),"(I2)") PRN
            if (name(1:1)=="G") then
                if (PRN>GNum) cycle
                ClkData%AS(PRN)%Clk(2)=clk
            else if (name(1:1)=="R") then
                if (PRN>RNum) cycle
                ClkData%AS(PRN+GNum)%Clk(2)=clk
            else if (name(1:1)=="C") then
                if (PRN>CNum) cycle
                ClkData%AS(PRN+GNum+RNum)%Clk(2)=clk
            else if (name(1:1)=="E") then
                if (PRN>NumE) cycle
                ClkData%AS(PRN+GNum+RNum+CNum)%Clk(2)=clk
            else if (name(1:1)=="J") then
                if (PRN>JNum) cycle
                ClkData%AS(PRN+GNum+RNum+CNum+NumE)%Clk(2)=clk
            else if (name(1:1)=="I") then
                if (PRN>INum) cycle
                ClkData%AS(PRN+GNum+RNum+CNum+NumE+JNum)%Clk(2)=clk
            end if  ! f (Name(1)=="G")
        end if   ! if (clk_type=="AS")
    end do  !  do while(.true.)
end subroutine