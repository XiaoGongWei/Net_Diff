!============ Sat_Diff ==========
!
! This function is used for satellite difference
!
!       L                    n*1 matrix
!       PRNPRN        PRN, n*1 matrix
!       Code             P1 or P2  
!       N                   PRN numbers
!       M                  parameter numbers
!       refSat            Reference Satellite

! Written by: Yize Zhang
!===========End of Header==========

subroutine Sat_Diff(A, L, PRNPRN, EleEle, N,M, RefSat)   
implicit none
    ! Intent
    integer :: N, M, RefSat
    real(8) :: A(N,M),L(N), PRNPRN(N), EleEle(N)
    ! Local variales
    integer :: i , RefSat_temp
    real(8) :: Maxele

    MaxEle=0.d0
    do i=1,N
        if ( (PRNPRN(i)==RefSat) .and. (EleEle(i)>=30.d0)) then
            goto 200
        elseif (EleEle(i)>MaxEle) then 
            RefSat_temp=PRNPRN(i)
            MaxEle=EleEle(i)
        end if
    end do
    RefSat=RefSat_temp

    200 do i=1,N
        

    end do
end 