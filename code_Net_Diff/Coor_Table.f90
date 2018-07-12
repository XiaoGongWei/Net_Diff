! ===========  Coor_Table  ============
!
! PURPOSE:
!        Read the coordinates of the stations and compute
!    the hydrostatic and wet zenith delay.
! 
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ============ End of Header ============

subroutine Coor_Table
use MOD_FileDir
use MOD_FileID
use MOD_CycleSlip
use MOD_Sta
use MOD_Ant
use MOD_Res
use MOD_ObsHead
use MOD_Var
implicit none
    logical :: alive
    character(300) :: line
    character(100) :: temp
    character(1) :: SKD
    integer(1) :: selectd, PRN
    integer :: N1,N2, Num, i,j, k
    real(8) :: Lat, Lon, Hgt,Coor(3)
    real(8) :: press ,tempre, rhumity, undu, MJD
    real(8) :: Lat2(2), Lon2(2), Hgt2(2), ah2(2), aw2(2),rh2(2), la2(2),undu2(2),press2(2), temp2(2),dT2(2),Tm2(2)
    
!    real(8) :: pi=3.14159265358979d0
    
    inquire(file=CoorTable,exist=alive)
    if (.not. alive) then
        write(*,*) "Coordinate table: """//trim(CoorTable)//""" doesn't exist!"
        pause
        stop
    end if
    CoorTableID=FileID_Mark
    FileID_Mark=FileID_Mark+1
    open(unit=CoorTableID, file=CoorTable)
    
    N1=0
    N2=0
    do while(.true.)
        read(unit=CoorTableID,fmt='(A)',end=200) line
        if (index(line,"!") /= 0)  line=line(1:index(line,"!")-1)
        if (len_trim(line)==0) cycle
        if ((index(line(1:1),"+") /= 0) .or. (index(line(1:1),"*") /= 0))   then
            cycle
        elseif (index(line(1:1),"-") /= 0)   then
            exit
        end if
        read(line,*) temp, SKD
        if ((SKD=="F") .or. (SKD=="f"))then  ! Station coordinate fixed
            N1=N1+1
        else if ((SKD=="K") .or. (SKD=="k"))then !  Station coordinate kinemetic
            N2=N2+1
        end if
    end do

    if (N2==0) then
        write(*,*) "Not found K in Coor_Table.txt, no station to be esimated."
        pause
        stop
    end if

    Num=N1+N2
    allocate(ObsID(Num))
    allocate(ObsHead(Num))
    STA%FixNum=N1
    STA%Num=Num
    if (N1>0) then
        allocate(Res(N1))
    else
         allocate(Res(1))
    end if
    allocate(STA%STA(Num))
    allocate(Ant(Num+GNum+RNum+CNum+NumE+JNum+INum))
!    allocate(STA%Name(Num), STA%SKD(Num))
!    allocate(STA%Coor(3,Num), STA%BLH(3,Num), STA%Rota(Num))
!    allocate(STA%Trop(Num))
    allocate(CycleSlip(Num))

    i=0
    j=0
    rewind(unit=CoorTableID)
    do while(.true.)
        read(unit=CoorTableID,fmt='(A)',end=200) line
        if (index(line,"!") /= 0)  line=line(1:index(line,"!")-1)
        if (len_trim(line)==0) cycle
        if ((index(line(1:1),"+") /= 0) .or. (index(line(1:1),"*") /= 0))   then
            cycle
        elseif (index(line(1:1),"-") /= 0)   then
            exit
        end if
        read(line,*) temp, SKD
        if ((SKD=="F") .or. (SKD=="f"))then  ! Station coordinate fixed
            i=i+1
            k=i
        else if ((SKD=="K") .or. (SKD=="k"))then !  Station coordinate kinemetic
            j=j+1
            k=N1+j
        end if  ! if ((SKD=="F") .or. (SKD=="f"))then  ! Statio
        read(line,*) STA%STA(k)%Name, STA%STA(k)%SKD,Coor, STA%STA(k)%NEU(3),STA%STA(k)%NEU(1),STA%STA(k)%NEU(2)  ! , STA%STA(k)%Ant_Type
        if (all(Coor/=0.d0)) then
            call XYZ2BLH(Coor(1), Coor(2), Coor(3),Lat,Lon,Hgt)  ! Lat and Lon is In degree
            STA%STA(k)%BLH(1:3)=(/Lat,Lon, Hgt/)
            STA%STA(k)%Rotation=reshape((/ -dsind(Lat)*dcosd(Lon),  -dsind(Lon),  dcosd(Lat)*dcosd(Lon),  &
                    -dsind(Lat)*dsind(Lon), dcosd(Lon) ,dcosd(Lat)*dsind(Lon),  dcosd(Lat) ,0D0, dsind(Lat)/), (/3,3/))
            STA%STA(k)%Coor(1:3)=Coor !+MATMUL(STA%STA(k)%NEU, STA%STA(k)%Rotation)
            STA%STA(k)%XYZ(1:3)=Coor
            STA%STA(k)%TrueCoor(1:3)=Coor
            STA%STA(k)%flag_InitialCoor=.true.
!            call XYZ2BLH(STA%STA(k)%Coor(1), STA%STA(k)%Coor(2), STA%STA(k)%Coor(3),Lat,Lon,Hgt)  ! Lat and Lon is In degree
!            STA%STA(k)%BLH(1:3)=(/Lat,Lon, Hgt/)

            ! Get the meteorological parameter 
            call UTC2MJD(int_year,1,int_doy, 12, 0, 0.d0, MJD)   ! In the middle of the day
            call gpt (MJD, Lat*pi/180.d0, Lon*pi/180.d0, Hgt, press ,tempre, rhumity, undu)
            if ( (cdattype(1:3)=='NOM') .or.  (trim(cdattype)=='GPT') ) then
                call gpt (MJD, Lat*pi/180.d0, Lon*pi/180.d0, Hgt, STA%STA(k)%TROP%press ,STA%STA(k)%TROP%temp, STA%STA(k)%TROP%e, STA%STA(k)%TROP%undu)
            elseif (cdattype(1:5)=='GPT2 ') then 
                Lat2(1)=Lat*pi/180.d0; Lon2(1)=Lon*pi/180.d0; Hgt2(1)=Hgt
                call gpt2(MJD, Lat2, Lon2, Hgt2, 1, 0,press2,temp2,undu2,rh2,ah2,aw2,undu2)
            elseif (cdattype(1:7)=='GPT2_1w') then 
                Lat2(1)=Lat*pi/180.d0; Lon2(1)=Lon*pi/180.d0; Hgt2(1)=Hgt
                call gpt2_1w(MJD, Lat2, Lon2, Hgt2, 1, 0,press2,temp2,dT2,Tm2,rh2,ah2,aw2,la2,undu2)
            elseif (cdattype(1:7)=='GPT2_5w') then 
                Lat2(1)=Lat*pi/180.d0; Lon2(1)=Lon*pi/180.d0; Hgt2(1)=Hgt
                call gpt2_5w(MJD, Lat2, Lon2, Hgt2, 1, 0,press2,temp2,dT2,Tm2,rh2,ah2,aw2,la2,undu2)
            else
                write(*,*) '***ERROR : unknown ZTD input meteo data : ',cdattype
                pause
                stop
            end if
            if (cdattype(1:4)=='GPT2') then
                STA%STA(k)%TROP%press=press2(1)
                STA%STA(k)%TROP%temp=temp2(1)
                STA%STA(k)%TROP%dT=dT2(1)
                STA%STA(k)%TROP%Tm=Tm2(1)
                STA%STA(k)%TROP%e=rh2(1)
                STA%STA(k)%TROP%ah=ah2(1)
                STA%STA(k)%TROP%aw=aw2(1)
                STA%STA(k)%TROP%la=la2(1)
            end if
            ! Calculate the tropsphere ZHD and ZWD
            if (cztd(1:5)=='EGNOS') then    ! EGNOS to get ZHD and ZTD
                call ZTD_EGNOS(Lat, Hgt, int_doy, STA%STA(k)%Trop%ZHD, STA%STA(k)%Trop%ZWD)
            else if (cztd(1:4)=='SAAS') then  ! SAAS
                call ZTD_SAAS(STA%STA(k)%TROP%press ,STA%STA(k)%TROP%temp, STA%STA(k)%TROP%e, Lat*pi/180.d0, Hgt, STA%STA(k)%TROP%ZHD,STA%STA(k)%TROP%ZWD)
                if (cdattype(7:7)=='w') then
                    call ZWD_Askne(STA%STA(k)%TROP%e,STA%STA(k)%TROP%Tm,STA%STA(k)%TROP%la,STA%STA(k)%TROP%ZWD)
                end if
            else if (cztd(1:5)=='UNB3 ') then  ! UNB3
                call ZTD_UNB3(Lat, Hgt, int_doy, STA%STA(k)%Trop%ZHD,STA%STA(k)%Trop%ZWD)
            elseif (cztd(1:5)=='UNB3m') then  ! UNB3m
                call ZTD_UNB3m(Lat, Hgt, int_doy, STA%STA(k)%Trop%ZHD,STA%STA(k)%Trop%ZWD)
            else
                write(*,*) '***ERROR : unknown ZTD model: ' ,cztd
                pause
                stop
            end if
        end if
    end do

    ! satellites select
    do while(.true.)
        read(unit=CoorTableID,fmt='(A)',end=200) line
        if (index(line,"!") /= 0)  line=line(1:index(line,"!")-1)
        if (len_trim(line)==0) cycle
        if ((index(line(1:10),"+satellite") /= 0) .or. (index(line(1:1),"*") /= 0))   then
            cycle
        elseif (index(line(1:10),"-satellite") /= 0)   then
            exit
        end if
        read(line,*) temp, selectd
        read(temp(2:3),'(I2)') PRN
        if (temp(1:1)=='G') then
            if (.not. SystemUsed(1)) cycle
        elseif (temp(1:1)=='R') then
            if (.not. SystemUsed(2)) cycle
            PRN=PRN+GNum
        elseif (temp(1:1)=='C') then
            if (.not. SystemUsed(3)) cycle
            PRN=PRN+GNum+RNum
        elseif (temp(1:1)=='E') then
            PRN=PRN+GNum+RNum+CNum
            if (.not. SystemUsed(4)) cycle
        elseif (temp(1:1)=='J') then
            if (.not. SystemUsed(5)) cycle
            PRN=PRN+GNum+RNum+CNum+NumE
        elseif (temp(1:1)=='I') then
            if (.not. SystemUsed(6)) cycle
            PRN=PRN+GNum+RNum+CNum+NumE+JNum
        else
            write(*,*) 'Format error in Coor_Table.txt, please check.'
            pause
            stop
        end if
        if (selectd==1) then  ! satellite used
            SatSelected(PRN)=1
        else  ! ! bad satellite
            SatSelected(PRN)=0
        end if
    end do

    200 close(CoorTableID)
    return
end subroutine