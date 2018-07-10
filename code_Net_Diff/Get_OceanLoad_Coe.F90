! ==========  Get Ocean Load Coefficient  ==========
!
! PURPOSE:
!         Get the ocean load coefficient file (i.e. Tide.txt) ,
!  & at the given latitude and longitude of the station.
!
! INPUTS:
!    Lat              Latitude of the station, unit in degree
!    Lon             Longitude of the station, unit in degree
!    name          Name of the station
!
! OUTPUTS:
!    CMC           Correction for Mass Center
!    OLC            Ocean Load Coefficient
!    found          Whether the station OLC is found
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ===============  End of Header  ==============

subroutine Get_OceanLoad_Coe(Lat,Lon,Name,CMC,OLC, found)
use MOD_FileID
use MOD_FileDir
implicit none
    ! Intent in
    real(8) :: Lat, Lon
    character(4)  :: name
    ! Intent out
    real(8) :: CMC(11,6),OLC(11,6)
    ! Local Variables
    logical :: alive, found
    character(120) :: line
    integer :: i
    real(8) :: temp(11,6)=0.d0, tempLat, tempLon
    
    inquire(file=TideFile,exist=alive)
    if (.not. alive) then
        write(*,*) "Ocean Load Coefficient file: """//trim(TideFile)//""" doesn't exist!"
        pause
        stop
    end if
    OLCID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=OLCID, file=TideFile)
    
    do while(.true.)
        read(unit=OLCID,fmt='(A)',end=200) line
        if (index(line,'+cmc')/=0) then
            i=0
            do while(.true.)
                read(unit=OLCID,fmt='(A)',end=200) line
                if (line(1:1)=='*') cycle
                if (line(1:4)=='-cmc') exit
                i=i+1
                read(line(11:),*) temp(i,1:6)
            end do
            if (i/=11) then
                write(*,*) 'Error in: +/-cmc must have 11 lines'
                pause
                stop
            end if
        end if
        if (index(line,"+oload") /= 0)  exit
    end do
    
    !! convert CMC coefficients into OLC coefficients, angle in radians
    !! CMC components are still xyz instead of radial, west south for stations
    do i=1,11
        CMC(i,1)=dsqrt(temp(i,3)**2+temp(i,4)**2)
        CMC(i,2)=dsqrt(temp(i,5)**2+temp(i,6)**2)
        CMC(i,3)=dsqrt(temp(i,1)**2+temp(i,2)**2)
        CMC(i,4)=datan2(temp(i,4),temp(i,3))
        CMC(i,5)=datan2(temp(i,6),temp(i,5))
        CMC(i,6)=datan2(temp(i,2),temp(i,1))
    end do
    
    found=.false.
    do while (.true.)
        read(unit=OLCID,fmt='(A)',end=200) line
        if (index(line,'lon/lat')/=0) then
            read(line,'(50X,F9.4, F10.4)') TempLon, TempLat
            if (dabs(TempLat-Lat)<0.25d0 .and. dabs(TempLon-Lon)<0.25d0) then
                found=.true.
                do i=1,6
                    read(unit=OLCID,fmt=*,end=200) temp(1:11,i) 
                end do
                exit
            end if
        end if
    end do
    
    
    ! The local system is defined as ENU. The coefficients are for radial, west, south.
    !  We adapted the array from RWS to ENU. See A_RTK.
    OLC(1:11,1)= -temp(1:11,2)    ! amplitude, W->E, 2->1, minus
    OLC(1:11,4)= -temp(1:11,5)    ! phase
    OLC(1:11,2)= -temp(1:11,3)    ! amplitude, S->N, 3->2, minus
    OLC(1:11,5)= -temp(1:11,6)    ! phase
    OLC(1:11,3)= -temp(1:11,1)    ! amplitude, R->U, 1->3, minus
    OLC(1:11,6)= -temp(1:11,4)    ! phase
    
    200 close(OLCID)
    if (.not. found) then
        write(*,*) "*****WARNING: Ocean Load Coefficient not found for station "//name//". It will considered as zero."
!        pause
        return
    end if
    
end subroutine