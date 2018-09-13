! ============== Cal_Sat_PosClk ================
!
! PURPOSE:
!        Entry of calculating the satellite position and satellitet clock.
! 
! INPUTS:
!          System                System flag, 'G',' R' or 'C'
!          Obsweek             obsweek
!          Obssec                obssec
!          PRN                     satellite number, integer
!          Coor                    Station Coordinate, unit in meter
!          t1                        signal transfer time using LC combined range, unit in second
!          Rela_flag            relativity flag (logical)
!                                     .true. : calculate relativity
!                                     .false. :  relativity=0
!
! OUTPUTS
!         Sat_Coor             satellite coordinate, unit in meter
!         Sat_Clk               satellite clock, unit in second
!         s                          distance between the satellite and the receiver, unit in meter
!         Rela                     relativity correction, unit in meter
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!    ===============End of header =================

subroutine Cal_Sat_PosClk(System, Obsweek, Obssec, PRN, Coor, t1, rela_flag,  Sat_Coor, Sat_Vel, Sat_Clk, s, Rela)
use MOD_Sp3Data
use MOD_ClkData
use MOD_FileID
use MOD_VAR
implicit none
    ! Intent in:
    character(1) :: System
    integer :: Obsweek, PRN
    real(8) :: Obssec, Coor(3), t1
    logical :: Rela_Flag
    ! Intent out:
    real(8) :: Sat_Coor(3), Sat_Vel(3), Sat_Clk, s , Rela
    ! Local variables
    integer :: error

    
    Sat_Coor=9999.d0
    Sat_Clk=9999.d0
    s=0.d0
    Rela=0.d0
    if (index(CLK,"CLK")/=0) then   ! If Precise Clock
        ! Detect whether to read 5-min interval clock data
        do while ( ( (Obsweek-ClkData%GPSweek(2))*604800.0d0+Obssec-ClkData%GPSsec(2) ) > 0.0d0)
            call ReadClkData()
            read(unit=ClkID,fmt="(2X)",iostat=error)  ! Check that if reach the end of the clock file
            backspace(ClkID)
            if (error /=0) exit
        end do
        if ( (( (Obsweek-ClkData%GPSweek(2))*604800.0d0+Obssec-ClkData%GPSsec(2) ) > 340.0d0)  .or. &
                (( (Obsweek-ClkData%GPSweek(1))*604800.0d0+Obssec-ClkData%GPSsec(1) ) < -340.0d0) ) then
                return
        end if
    end if
    if (index(Orbit,"SP3")/=0) then   ! If Precise Orbit
            ! Detect whether to read 15-min interval SP3 data
        do while ( ( (Obsweek-SP3Data%GPSweek(6))*604800.0d0+Obssec-SP3Data%GPSsec(6) ) > 0.0d0)
            call ReadSP3Data(1)
            read(unit=SP3ID,fmt="(2X)",iostat=error)  ! Check that if reach the end of the SP3 file
            if (error /=0) exit
            backspace(SP3ID)
!            call ReadMultiOrb(1,SP3Data%GPSweek(10), SP3Data%GPSsec(10)+300.d0)
        end do
        if ( (( (Obsweek-SP3Data%GPSweek(10))*604800.0d0+Obssec-SP3Data%GPSsec(10) ) > 40.0d0)  .or. &
                (( (Obsweek-SP3Data%GPSweek(1))*604800.0d0+Obssec-SP3Data%GPSsec(1) ) < -40.0d0) ) then
                return
        end if
    end if

    if ( (index(Orbit,"SP3")/=0) .and. (index(Clk,"SP3")/=0) ) then   ! If Precise Orbit
        call Cal_Sat_Pos(Obsweek, Obssec+ObsTime+SP3Time, PRN, Coor, t1, rela_flag , Sat_Coor, Sat_Vel, s, Rela, Sat_Clk)
        return
    end if

    if (index(Orbit,"SP3")/=0) then   ! If Precise Orbit
        call Cal_Sat_Pos_sp3(Obsweek, Obssec+ObsTime+SP3Time, PRN, Coor, t1, rela_flag, Sat_Coor, Sat_Vel, s, Rela)
    elseif  (index(Orbit,"BRD")/=0) then
        if (System=="G") then   ! GPS
            call Cal_Sat_Pos_n(System, Obsweek, Obssec+ObsTime, PRN, Coor, t1, Sat_Coor, Sat_Vel, s, Rela)
        elseif (System=="R") then   ! GLONASS
            call Cal_Sat_Pos_g(Obsweek, Obssec+ObsTime-leap_sec, PRN-GNum, Coor, t1,  Sat_Coor, Sat_Vel, s, Rela, Sat_Clk)
        elseif (System=="C") then   ! COMPASS
            call Cal_Sat_Pos_c(Obsweek, Obssec+ObsTime-14.d0, PRN-GNum-RNum, Coor, t1, Sat_Coor, Sat_Vel, s, Rela, Sat_Clk)
        elseif (System=="E") then   ! GALILEO
            call Cal_Sat_Pos_e(Obsweek, Obssec+ObsTime, PRN-GNum-RNum-CNum, Coor, t1, Sat_Coor, Sat_Vel, s, Rela)
        elseif (System=="J") then   ! QZSS
            call Cal_Sat_Pos_n(System, Obsweek, Obssec+ObsTime, PRN-GNum-RNum-CNum-NumE, Coor, t1, Sat_Coor, Sat_Vel, s, Rela)
        elseif (System=="I") then   ! IRNSS
            call Cal_Sat_Pos_n(System, Obsweek, Obssec+ObsTime, PRN-GNum-RNum-CNum-NumE-JNum, Coor, t1, Sat_Coor, Sat_Vel, s, Rela)
        else
            write(*,*) 'Error in other system by using navigation ephemeris.'
            pause
            stop
        end if
    end if
    if (index(Clk,"CLK")/=0) then   ! If Precise Clock
        call Cal_Sat_Clk(Obsweek, Obssec+ObsTime+SP3Time-t1,PRN, Sat_Clk)
    elseif  (index(Clk,"BRD")/=0) then
        if (System=="G") then   ! GPS
            call Cal_Sat_Clk_n(System, Obsweek, Obssec+ObsTime-t1, PRN, Sat_Clk)
        elseif (System=="R") then   ! GLONASS

        elseif (System=="C") then   ! COMPASS
             call Cal_Sat_Clk_c(Obsweek, Obssec+ObsTime-14.d0-t1, PRN-GNum-RNum, Sat_Clk)
        elseif (System=="E") then   ! GALILEO
            call Cal_Sat_Clk_e(Obsweek, Obssec+ObsTime-t1, PRN-GNum-RNum-CNum, Sat_Clk)
        elseif (System=="J") then   ! QZSS
            call Cal_Sat_Clk_n(System, Obsweek, Obssec+ObsTime-t1, PRN-GNum-RNum-CNum-NumE, Sat_Clk)
        elseif (System=="I") then   ! IRNSS
            call Cal_Sat_Clk_n(System, Obsweek, Obssec+ObsTime-t1, PRN-GNum-RNum-CNum-NumE-JNum, Sat_Clk)
        else
            write(*,*) 'Error in other system by using navigation ephemeris.'
            pause
            stop
        end if
    end if

end subroutine