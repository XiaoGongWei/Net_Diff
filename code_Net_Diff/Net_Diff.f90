!============ Net_Diff ==========
!
! This function is used for net difference
!
!       L                    n*1 matrix
!       PRNPRN        PRN, n*1 matrix
!       Code             P1 or P2  
!       N                   PRN numbers
!       M                  parameter numbers
!       s                the order of the station

! Written by: Yize Zhang
!===========End of Header==========

subroutine Net_Diff(A,L,PRNPRN,Code,N,M,Num,S,epoch)
use MOD_STA
use MOD_Res
use MOD_FileID
use MOD_VAR
implicit none
    ! Intent
    integer :: N,Num,M, S
    integer :: PRNPRN(N)
    real(8) :: L(N), A(N,M)
    character(2) :: Code(N)
    ! Local variables
    integer :: i, j,k
    real(8) :: BLH(3)
    integer :: PRN
    character(2) :: Co
    real(8) :: P, sum_P
    real(8) :: dL, xx
    integer  :: epoch

    
!    BLH=STA%STA(sta)%BLH(1:3)
    do i=1,N
        PRN=PRNPRN(i)
        Co=Code(i)
        dL=0.d0
        sum_P=0.d0
!        do j=1,Sta%FixNum  ! For each fixed station
             P=1.d0
             j=1
             do k=1,Res(j)%N   ! Find the same PRN
                if ((Res(j)%PRN(k)==PRN) .and. (Res(j)%Code(k)==Co)) then
!                if ((Res(j)%PRN(k)==PRN) .and. (Res(j)%Code(k)=='LC')) then
                    dL=dL+P*Res(j)%L(k)
                    sum_P=sum_P+P
                    exit
                end if
            end do  ! do k=1,Res(j)%N
!        end do  ! do j=1,Sta%FixNum
        if (sum_P/=0.d0) then
            dL=dL/sum_P
            
            if (proc_mod==3)then
                if (Co=='LC') L(i)=L(i)-dL  ! For zone correction in Process_Corr
            else
                L(i)=L(i)-dL  ! In DSPP/DPPP
            end if
        else    ! If no common satellites
            L(i)=0.d0
            A(i,:)=0.d0
            Num=Num-1
        end if
        
    end do  ! do i=1,N

    return
end subroutine