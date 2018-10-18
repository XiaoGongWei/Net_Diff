! ============= Station Difference ===============
! Purpose:
!     Form station difference equation for each station
!
! INPUT:
!         ZD                  zero difference structure
! OUTPUT:
!         SD                  station difference structure
!
! WRITTEN BY: Yize Zhang
! ==========================================

subroutine Station_Difference(ZD, SD)
use MOD_ZD
use MOD_SD
use MOD_FileID
use MOD_VAR
implicit none
    type(type_ZD)  :: ZD(2)
    type(type_SD)  :: SD
    ! Local variables
    integer :: i, j, N
    

    SD%PRN=0;       SD%A=0.d0
    SD%Ele=0.d0;   SD%Q=0.d0;   
    SD%P1=0.d0;    SD%P2=0.d0
    SD%P1CS=0.d0;    SD%P2CS=0.d0
    SD%L1=0.d0;    SD%L2=0.d0
    SD%WL=0.d0;  SD%W4=0.d0
    SD%EWL=0.d0; SD%WL_amb=99.d0

    ! ******** Start to do station difference *******
    N=0
    do i=1,ZD(1)%PRNS  ! in the user station
        do j=1,ZD(2)%PRNS   ! in the reference station
            if (ZD(1)%PRN(i)==ZD(2)%PRN(j)) then  ! find the same satellite
                N=N+1
                SD%Sys(N)   =    ZD(2)%Sys(j)
                SD%System(N)   =    ZD(2)%System(j)
                SD%PRN(N)   =    ZD(2)%PRN(j)
                SD%PRN_S(N)   =    ZD(2)%PRN_S(j)
                SD%Ele(N)     =    0.5d0*ZD(1)%Ele(i)+0.5d0*ZD(2)%Ele(j)
                SD%Q(N)       =    0.5d0/ZD(1)%P(i)+0.5d0/ZD(2)%P(j)
                SD%A(N,:)      =    ZD(2)%A(j,:)
                SD%Corr(N)   =    ZD(2)%Corr(j)-ZD(1)%Corr(1)
                SD%s(N)       =    ZD(2)%s(j)-ZD(1)%s(i)
                if ( (ZD(1)%P1(i)/=0.d0) .and. (ZD(2)%P1(j)/=0.d0) ) then
                    SD%P1(N)=ZD(2)%P1(j)-ZD(1)%P1(i)
                end if
                if ( (ZD(1)%P2(i)/=0.d0) .and. (ZD(2)%P2(j)/=0.d0) ) then
                    SD%P2(N)=ZD(2)%P2(j)-ZD(1)%P2(i)
                end if
                if ( (ZD(1)%P1CS(i)/=0.d0) .and. (ZD(2)%P1CS(j)/=0.d0) ) then
                    SD%P1CS(N)=ZD(2)%P1CS(j)-ZD(1)%P1CS(i)   ! for cycle slip
                end if
                if ( (ZD(1)%P2CS(i)/=0.d0) .and. (ZD(2)%P2CS(j)/=0.d0) ) then
                    SD%P2CS(N)=ZD(2)%P2CS(j)-ZD(1)%P2CS(i)
                end if
                if ( (ZD(1)%L1(i)/=0.d0) .and. (ZD(2)%L1(j)/=0.d0) ) then
                    SD%L1(N)=ZD(2)%L1(j)-ZD(1)%L1(i)
                end if
                if ( (ZD(1)%L2(i)/=0.d0) .and. (ZD(2)%L2(j)/=0.d0) ) then
                    SD%L2(N)=ZD(2)%L2(j)-ZD(1)%L2(i)
                end if
                if ( (ZD(1)%WL(i)/=0.d0) .and. (ZD(2)%WL(j)/=0.d0) ) then
                    SD%WL(N)=ZD(2)%WL(j)-ZD(1)%WL(i)
                end if
                if ( (ZD(1)%WL_amb(SD%PRN(N))/=99.d0) .and. (ZD(2)%WL_amb(SD%PRN(N))/=99.d0) ) then  ! Just for test, not very good, because of the wrong rounding integer
!                    SD%WL_amb(N)=ZD(2)%WL_amb(j)-ZD(1)%WL_amb(i)  ! Wide Lane ambiguity, in cycle
                    SD%WL_amb(SD%PRN(N))=ZD(2)%WL_amb(SD%PRN(N))-ZD(1)%WL_amb(SD%PRN(N))  ! Wide Lane ambiguity, in cycle
                end if
                if ( (ZD(1)%W4(i)/=0.d0) .and. (ZD(2)%W4(j)/=0.d0) ) then
                    SD%W4(N)=ZD(2)%W4(j)-ZD(1)%W4(i)
                end if
                if ( (ZD(1)%EWL(i)/=0.d0) .and. (ZD(2)%EWL(j)/=0.d0) ) then
                    SD%EWL(N)=ZD(2)%EWL(j)-ZD(1)%EWL(i)
                    SD%EWL_amb(N)=ZD(2)%EWL_amb(j)-ZD(1)%EWL_amb(i)
                end if
                write(LogID,'(A6,1X,A1,I2,2F8.3,4F13.3)') '##SD', SD%System(N), SD%PRN_S(N),SD%P1(N),SD%P2(N), &
                                SD%L1(N),SD%L2(N),SD%WL(N),SD%W4(N)
                exit
            end if  ! if (ZD(1)%PRN(i)==ZD(2)%PRN(j)) then
        end do 
    end do ! do i=1,ZD(1)%PRNS

    SD%PRNS =  N
    SD%week =  ZD(2)%week
    SD%sow   =  ZD(2)%sow
    return
end subroutine