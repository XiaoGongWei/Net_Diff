

module MOD_IONHead
    implicit none
    type type_IONhead
        real(8) :: Hgt1, Hgt2, dHgt
        real(8) :: Lat1, Lat2, dLat
        integer :: nLat
        real(8) :: Lon1, Lon2, dLon
        integer :: nLon
        integer :: EXPONENT, INTERVAL
        integer :: MAPS
    end type
    type(type_IONhead), save :: IONHead
end module

module MOD_IONData
    implicit none
    type type_IONData
        integer :: Week
        real(8) :: sow
        real(8), allocatable :: TEC(:,:)
    end type
    type(type_IONData), allocatable :: IONData(:)
end module


! ============= ReadIONEX =================
! PURPOSE:
!         This function is to read the IONEX format file
! 
! Written by: Yize Zhang
!============== End of Header =============

subroutine ReadIONEX(IONFile)
use MOD_FileID
use MOD_IONHead
use MOD_IONData
implicit none
    character(200) :: IONFile
    logical :: alive
    character(100) :: line
    character(len=30) keyword
    integer :: imap, i, j, k
    integer :: year, mon, Iday, hour, min, sec

    ! inquire that if the obs file exists
    inquire(file=IONFile,exist=alive)
    if (.not. alive) then
        write(*,*) "IONEX flie: """//trim(IONFile)//""" doesn't exist!"
        pause
        return
    end if
    
    IONID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=IONID, file=IONFile,action="read")
    
    ! read ION header
    do while(.true.)
        read(unit=IONID,fmt='(A)', end=200) line
        keyword=line(61:80)
        if (index(keyword,"HGT1 / HGT2 / DHGT") /= 0) then
            read(unit=line,fmt='(2X,3F6.1)') IONHead%Hgt1,  IONHead%Hgt2,  IONHead%dHgt
        elseif (index(keyword,"LAT1 / LAT2 / DLAT") /= 0) then
            read(unit=line,fmt='(2X,3F6.1)') IONHead%Lat1,  IONHead%Lat2,  IONHead%dLat
            IONHead%nLat=int(abs((IONHead%Lat1-IONHead%Lat2)/ IONHead%dLat))+1
        elseif (index(keyword,"LON1 / LON2 / DLON") /= 0) then
            read(unit=line,fmt='(2X,3F6.1)') IONHead%Lon1,  IONHead%Lon2,  IONHead%dLon
            IONHead%nLon=int(abs((IONHead%Lon1-IONHead%Lon2)/ IONHead%dLon))+1
        elseif (index(keyword,"EXPONENT") /= 0) then
            read(unit=line,fmt='(I6)') IONHead%EXPONENT
        elseif (index(keyword,"INTERVAL") /= 0) then
            read(unit=line,fmt='(I6)') IONHead%INTERVAL
        elseif (index(keyword,"# OF MAPS IN FILE") /= 0) then
            read(unit=line,fmt='(I6)') IONHead%Maps
        elseif (index(keyword,"END OF HEADER") /= 0) then
            exit
        end if
    end do
    
    allocate(IONData(IONHead%Maps))
    do i=1, IONHead%Maps
         allocate(IONData(i)%TEC(IONHead%nLat, IONHead%nLon))
    end do 

    ! read ION data
    do while(.true.)
        read(unit=IONID,fmt='(A)', end=200) line
        keyword=line(61:80)
        if (index(keyword,"START OF TEC MAP") /= 0) then
            read(unit=line,fmt='(I6)') imap
            i=0
        elseif (index(keyword,"EPOCH OF CURRENT MAP") /= 0) then
            read(unit=line,fmt='(5I6)') year, mon, Iday, hour, min !, sec
            call UTC2GPST(year, mon, Iday, hour, min, 0.d0, IONData(imap)%week,IONData(imap)%sow)
        elseif (index(keyword,"LAT/LON1/LON2/DLON/H") /= 0) then
            i=i+1
            do j=1,IONHead%nLon ! fix(IONHead%nLon/16.d0)+1
                k=mod(j,16)
                if (k==0) k=16
                if (k==1) then
                    read(unit=IONID,fmt='(A)', end=200) line
                end if
                read(line(k*5-4:k*5),'(F5.0)') IONData(imap)%TEC(i, j)
            end do
        elseif (index(keyword,"END OF TEC MAP") /= 0) then
            read(line,'(I6)') imap
            if (imap>=IONHead%maps) exit
        end if
    end do

    200 return
end subroutine
