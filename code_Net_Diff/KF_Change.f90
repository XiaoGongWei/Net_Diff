! ========================  KF_Change  ======================
!
! PURPOSE:
!           Change Inverse Nbb(InvN) in the NEQ(Normal EQuation) in Kalman Fileter
! 
! REFERENCE:
!          http://www.navipedia.net/index.php/Kalman_Filter
!          http://www.navipedia.net/index.php/Matrix_Definitions_Phi_and_Q
!          http://www.navipedia.net/index.php/Parameters_adjustment_for_PPP
! 
! INPUTS:
!           InvN               The normal equation  x=InvN*U, InvN is a low triangular matrix
!           X                    The normal equation  x=InvN*U
!           N                    The size of the normal equation
!           t                      Eliminate the parameter in the order of t
!
! OUTPUTS:
!              Nbb, U         The new normal equation.
! 
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! =================  End of header  =================

subroutine KF_Change(InvN,X, N,t,flag)
use MOD_VAR
implicit none
    ! Intent in and out
    real(8) :: InvN(N,N), X(N)
    integer :: N, t
    character(3) :: flag

    ! see Parameters adjustment for PPP in Navipedia
    ! http://www.navipedia.net/index_php/Parameters_adjustment_for_PPP
    if (index(flag,'pos') /= 0) then   ! Q 
        InvN(t,1:t)=0.d0   ! Eliminate InvN in row t and line t
        InvN(t:n,t)=0.d0
        InvN(t,t)=60.d0**2 ! 1.d8
        X(t)=0.d0
    elseif (flag(1:3)=='zwd') then
        if (InvN(t,t)==0.d0) InvN(t,t)= 0.05d0**2    ! P0  variance 0.5m
!        InvN(t,t)=InvN(t,t)+1.d-4*Interval/3600.d0 ! added in Zero_Difference.f90
    elseif (flag(1:3)=='clk') then
        InvN(t,1:t)=0.d0   ! Eliminate InvN in row t and line t
        InvN(t:n,t)=0.d0
        InvN(t,t)=9.d10   ! P0 variance 3*10d5m
        X(t)=0.d0
    elseif (flag(1:3)=='amb') then
        InvN(t,1:t)=0.d0   ! Eliminate InvN in row t and line t
        InvN(t:n,t)=0.d0
        InvN(t,t)=4.d2    ! P0 variance 20m
        X(t)=0.d0
    elseif (flag(1:3)=='ddp') then  ! double difference position
        InvN(t,1:t)=0.d0
        InvN(t:n,t)=0.d0
        InvN(t,t)=10.d0**2*10000.d0   ! precision is 10m
        X(t)=0.d0
    elseif (flag(1:3)=='dda') then  ! double difference ambiguity
        InvN(t,1:t)=0.d0
        InvN(t:n,t)=0.d0
        X(t)=0.d0
    elseif (flag(1:3)=='ddi') then  ! double difference ionosphere
        InvN(t,1:t)=0.d0
        InvN(t:n,t)=0.d0
        InvN(t,t)=0.2d0**2*10000.d0   ! precision is 0.2m
        X(t)=0.d0
    else 
        write(*,*) 'nargin error in KF_Change.f90'
    end if
    return
end subroutine
    