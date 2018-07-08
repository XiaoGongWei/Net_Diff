! ==============================================
!             Statiion Merge
! 
! PURPOSE:
!          Station merge in net difference. Merge the 1:n-1 th
!  zero-difference to one zero-difference
!
! INPUTS:
!     ZD                     raw zero difference strcuture, n*1 array
!     ZD_M                Zero difference strcuture after merge, 2*1 array
!     n                        total station number
!
! WRITTEN BY: Yize Zhang
! ============================================

subroutine Station_Merge(ZD, ZD_M, n)
use MOD_ZD
implicit none
    type(type_ZD) :: ZD(n), ZD_M(2)
    integer :: n
    ! Local variables
    integer :: i, j, k, PRN
    real(8) :: dL, Sum_P1, Sum_P2, P, dis

    ZD_M(2)=ZD(n)  ! User zero difference
    ZD_M(1)%PRNs=ZD_M(2)%PRNs
    ZD_M(1)%week=ZD_M(2)%week
    ZD_M(1)%sow=ZD_M(2)%sow

    do i=1, ZD_M(2)%PRNS
        PRN=ZD_M(2)%PRN(i)
        ZD_M(1)%PRN(i)=ZD_M(2)%PRN(i)
        ZD_M(1)%Ele(i)=ZD_M(2)%Ele(i)
        ZD_M(1)%P(i)=ZD_M(2)%P(i)
        ZD_M(1)%P1(i)=0.d0
        ZD_M(1)%P2(i)=0.d0
        Sum_P1=0.d0
        Sum_P2=0.d0
        do j=1, n-1
            P=1.d0   !!!!! According to distance??
            do k=1, ZD(j)%PRNS
                if (ZD(j)%PRN(k)==PRN) then
                    if ( ZD(j)%P1(k)/=0.d0 ) then
                        Sum_P1=Sum_P1+P
                        ZD_M(1)%P1(i)=ZD_M(1)%P1(i)+ZD(j)%P1(k)
                    end if
                    if ( ZD(j)%P2(k)/=0.d0 ) then
                        Sum_P2=Sum_P2+P
                        ZD_M(1)%P2(i)=ZD_M(1)%P2(i)+ZD(j)%P2(k)
                    end if
                    exit
                end if
            end do   ! do k=1, ZD(j)%PRNS
        end do   ! do j=1, n-1
        if (Sum_P1/=0.d0) then
            ZD_M(1)%P1(i)=ZD_M(1)%P1(i)/Sum_P1
        end if
        if (Sum_P2/=0.d0) then
            ZD_M(1)%P2(i)=ZD_M(1)%P2(i)/Sum_P2
        end if
    end do   ! do i=1, ZD_M(2)%PRNS

    return
end subroutine
