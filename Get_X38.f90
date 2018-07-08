!============ Get_X38 ==========
!
! PURPOSE:
!      This function is used to get X38
!
! INPUTS:
!       ObsWeek           
!       ObsSow
!
! Written by: Yize Zhang
!===========End of Header==========

subroutine Get_X38(ObsWeek, ObsSow) 
use MOD_Res
use MOD_VAR
use MOD_FileID
implicit none
    integer :: ObsWeek, Week
    real(8) :: ObsSow, Sow, tempsow1, tempsow2
    integer :: i, error
    character(100) :: line
    integer :: PRN, IQt, tempI, zone
    real(8) :: corr
    character(2)  :: code

    tempsow2=-1.d0
    Res(1)%N=0
    Res(1)%PRN=0
    Res(1)%Code=''
    Res(1)%L=0.d0
    Res(1)%CS=0
    do while(.true.)
        read(ResO_CID,fmt='(A)', iostat=error) line
        if (error /=0) exit
        read(line,*) week, sow, zone, IQt, tempI, PRN,code,corr
        week=week+1356
        tempsow1=sow
        if (delay/=0.d0) then
            sow=ceiling(sow/delay)*delay
        end if
        if ( (Week-ObsWeek)*604800.d0+(Sow-ObsSow)>0.d0) then
            backspace(ResO_CID)    ! Íí£¬ÍË³ö
            exit
        elseif ( (Week-ObsWeek)*604800.d0+(Sow-ObsSow)<0.d0) then
            cycle ! Ôç£¬¼ÌÐø
        end if

        if (tempsow1>tempsow2) then
            Res(1)%L=0.d0
            tempsow2=tempsow1
        end if
        if ((IQt==IorQ) .and. (index(line,'LC')/=0)) then    ! only for LC
            Res(1)%PRN(PRN*2)=PRN
            Res(1)%Code(PRN*2)=code
            Res(1)%L(PRN*2)=corr
            Res(1)%PRN(PRN*2-1)=PRN
            Res(1)%Code(PRN*2-1)='PC'
            Res(1)%L(PRN*2-1)=0.d0
            Res(1)%sow=sow
        end if
    end do
    Res(1)%N=70  ! only for beidou online
    
end subroutine