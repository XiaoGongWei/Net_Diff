!==================== Get_ESC ==================
!
! PURPOSE:
!      This function is used to get Equivalent Satllite Clock
!
! INPUTS:
!       ObsWeek           
!       ObsSow
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!==================== End of Header ==============

subroutine Get_ESC(ObsWeek, ObsSow) 
use MOD_ESC
implicit none
    integer :: ObsWeek, Week
    real(8) :: ObsSow, Sow
    integer(1) :: PRN, error
    character(100) :: line, temp

    ESC=9999.d0
    do while(.true.)
        read(100,*, iostat=error) line, temp, week, sow
        if (error /=0) exit
        if ( (Week-ObsWeek)*604800.d0+(Sow-ObsSow)>0.d0) then
            backspace(100)
            exit
        end if
        do while (.true.)
            read(100,'(A)', iostat=error) line
            if (error /=0) exit
            if (index(line,'EPO')/=0) then
                backspace(100)
                exit
            end if
            if (index(line,"CLK  G") /= 0)   then
                read(line(12:13),'(I2)') PRN
                read(line(14:28),'(F15.8)') ESC(PRN)
            end if
        end do
    end do

end subroutine