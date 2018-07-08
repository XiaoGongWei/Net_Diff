! ==================RNXSMT=====================
! Function:
!          RINEX SMooThing, smooth the pseudorange observation
!    
! Reference: 
!          Bernase GPS Software Version 5.0 DRAFT, Page90-94
!
! Inputs:
!          Epoch         Epoch flag
!          PRN             PRN of the satellite
!          Range         LC combination pseurange
!          L1, L2          Phase observation, in meters
! 
! Called by:
!          Process.F90
!
! Written by: Yize Zhang
! ==================End of Header=====================

! ==========Ambigulty module========
module MOD_Amb
implicit none
    type type_Amb
        integer ::Epo(2,91)=0
        real(8) :: Amb(91)=0.d0
    end type
    type(type_Amb),save :: Amb
end module

subroutine RNXSMT(epoch, PRN, Range,L1,L2)
use MOD_FileID
use MOD_Amb
use MOD_constant
implicit none
    ! Intent
    integer :: epoch, PRN
    real(8) :: Range, L1, L2
    ! Local variables
    integer :: error
    real(8) :: tempamb, temp
    
!    if (Amb%Amb(PRN)==0.d0) then
!        read(unit=AmbID(PRN),fmt='(4X,2I6,F19.3)' ) Amb%Epo(1,PRN), Amb%Epo(2,PRN), Amb%Amb(PRN)
!    end if
    ! Get cycle compensation
    do while (epoch > Amb%Epo(2,PRN) )
        read(unit=AmbID(PRN),fmt="(2X)",iostat=error)  ! Check that if reach the end of the file
        backspace(AmbID(PRN))
        if (error ==0) then
            read( AmbID(PRN),'(4X,2I6,F19.3)' ) Amb%Epo(1,PRN),Amb%Epo(2,PRN), Amb%Amb(PRN)
        end if
        if (epoch <= Amb%Epo(2,PRN)) exit
    end do
    
    if (Amb%Epo(2,PRN)-Amb%Epo(1,PRN)<5) then    ! 
!        cycle        ! if the arc less than 5 epochs, not use
        return
!    elseif (epoch<Amb%Epo(1,PRN)) then
!!        cycle       ! if the epoch is not in this arc
!        return
    else
        tempamb=Amb%Amb(PRN)
    end if
                
    if ((L1 /=0.d0) .and. (L2/=0.d0)) then
!        temp=L1+1.54572778d0*(L1-L2)+tempamb   ! 平滑后的无电离层伪距
        temp=(f1*f1*L1-f2*f2*L2)/(f1+f2)/(f1-f2) +tempamb  ! 平滑后的无电离层距离
!        if (dabs(Range-temp)>2.d0) then
!            write(*,'(I4,I4,F12.3,A15)') epoch, PRN, Range-temp, 'Range-Phase'
!        end if
        Range=temp
    end if
end subroutine