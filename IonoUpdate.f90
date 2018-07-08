! ================ IonoUpdate ==============
! Function: 
!              This subroutine is to update the residuals
!  in double difference(DD) OMC, which mainly comes
!  from ionosohere difference in two stations
!
! INPUTS:
!      Week           ObsWeek
!       Sow             ObsSow
!      Epo_NEQ     Epo_NEQ struct
!
! WRITTEN BY: Yize Zhang, zhyize@163.com
! ================ End of Header ===============

subroutine IonoUpdate(Week, Sow, NEQ, Epo_NEQ)
use MOD_Epo_NEQ
use MOD_NEQ
use MOD_IonoDDRes
implicit none
    type(type_NEQ) :: NEQ
    type(type_Epo_NEQ) :: Epo_NEQ
    integer :: Week
    real(8) :: Sow
    integer(2) :: i, PRN

    ! 1. First update residuals
    if (any(Epo_NEQ%Vl1/=0.d0) .or. any(Epo_NEQ%Vl2/=0.d0)) then
        do i=1,Epo_NEQ%PRNS
            PRN=Epo_NEQ%PRN(i)
            IonoDDRes(PRN)%Sow=eoshift(IonoDDRes(PRN)%Sow,shift=1,boundary=Sow)
            if (Epo_NEQ%Vl1(i)/=0.d0) then
!                IonoDDRes(PRN)%Week1=eoshift(IonoDDRes(PRN)%Week1,shift=1,boundary=Week)
!                IonoDDRes(PRN)%Sow1=eoshift(IonoDDRes(PRN)%Sow1,shift=1,boundary=Sow)
                IonoDDRes(PRN)%L1=eoshift(IonoDDRes(PRN)%L1,shift=1,boundary=Epo_NEQ%Vl1(i)+IonoDDRes(PRN)%dL1)
                IonoDDRes(PRN)%amb1=eoshift(IonoDDRes(PRN)%amb1,shift=1,boundary=NEQ%amb_L1(PRN))
            end if
            if (Epo_NEQ%Vl2(i)/=0.d0) then
!                IonoDDRes(PRN)%Week2=eoshift(IonoDDRes(PRN)%Week2,shift=1,boundary=Week)
!                IonoDDRes(PRN)%Sow2=eoshift(IonoDDRes(PRN)%Sow2,shift=1,boundary=Sow)
                IonoDDRes(PRN)%L2=eoshift(IonoDDRes(PRN)%L2,shift=1,boundary=Epo_NEQ%Vl2(i)+IonoDDRes(PRN)%dL2)
                IonoDDRes(PRN)%amb2=eoshift(IonoDDRes(PRN)%amb2,shift=1,boundary=NEQ%amb_L2(PRN))
            end if
        end do
    end if

    return
end subroutine
