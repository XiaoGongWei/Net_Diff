!============ Get_Corr ==========
!
! PURPOSE:
!      This function is used to get corrction information 
!   of pseudorange and carrier phase 
!
! INPUTS:
!       ObsWeek           
!       ObsSow
!
! Written by: Yize Zhang
!===========End of Header==========

subroutine Get_Corr(ObsWeek, ObsSow) 
use MOD_Res
use MOD_FileID
implicit none
    integer :: ObsWeek, Week
    real(8) :: ObsSow, Sow
    integer :: i, error
    character(100) line
    integer :: PRN
    real(8) :: L, x
    character(2)  :: code
    integer :: CS, CS_Array(32)

    CS_Array=0
    Res(1)%N=0
    Res(1)%PRN=0
    Res(1)%Code=''
    Res(1)%L=0.d0
    i=0
    do while(.true.)
        read(ResO_CID,*, iostat=error) line, week, sow
        if (error /=0) exit
        if ( (Week-ObsWeek)*604800.d0+(Sow-ObsSow)>0.d0) then
            backspace(ResO_CID)
            exit
        elseif ( (Week-ObsWeek)*604800.d0+(Sow-ObsSow)<0.d0) then
            do while (.true.)
                read(ResO_CID,'(A)', iostat=error) line
                if (error /=0) exit
                if (index(line,'*')/=0) then
                    backspace(ResO_CID)    ! new epoch
                    exit
                end if
                read(line,*) Code, PRN, L,CS
                if (CS==1) CS_Array(PRN)=1
            end do
        else
            i=0
            do while (.true.)
                read(ResO_CID,'(A)', iostat=error) line
                if (error /=0) exit
                if (index(line,'*')/=0) then
                    backspace(ResO_CID)    ! new epoch
                    exit
                end if
                i=i+1
                read(line,*) Res(1)%Code(i), Res(1)%PRN(i),Res(1)%L(i), Res(1)%CS(i)
                call random_number(x)
                Res(1)%L(i)=Res(1)%L(i) !+x*0.03d0-0.015d0
                if ( CS_Array( Res(1)%PRN(i) )==1 ) Res(1)%CS(i)=1
            end do
            Res(1)%sow=sow
        end if

    end do
    Res(1)%N=i
    
end subroutine