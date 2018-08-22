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
!           InvN               The normal equation  x=InvN*U
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
        InvN(t,:)=0.d0   ! Eliminate InvN in row t and line t
        InvN(:,t)=0.d0
        InvN(t,t)=60.d0**2   ! precision is 60m, in case of real kinematic
        X(t)=0.d0
    elseif (flag(1:3)=='zwd') then
        if (InvN(t,t)==0.d0) InvN(t,t)= 0.05d0**2    ! P0  variance 0.5m
!        InvN(t,t)=InvN(t,t)+1.d-4*Interval/3600.d0 ! added in Zero_Difference.f90
    elseif (flag(1:3)=='clk') then
        InvN(t,:)=0.d0   ! Eliminate InvN in row t and line t
        InvN(:,t)=0.d0
        InvN(t,t)=9.d10   ! P0 variance 3*10d5m
        X(t)=0.d0
    elseif (flag(1:3)=='amb') then
        InvN(t,:)=0.d0   ! Eliminate InvN in row t and line t
        InvN(:,t)=0.d0
        InvN(t,t)=4.d2    ! P0 variance 20m
        X(t)=0.d0
    elseif (flag(1:3)=='ddp') then  ! double difference position
        InvN(t,:)=0.d0
        InvN(:,t)=0.d0
        InvN(t,t)=60.d0**2   ! precision is 60m, in case of real kinematic
        X(t)=0.d0
    elseif (flag(1:3)=='dda') then  ! double difference ambiguity
        InvN(t,:)=0.d0
        InvN(:,t)=0.d0
        InvN(t,t)=1.d8
        X(t)=0.d0
    elseif (flag(1:3)=='ddi') then  ! double difference ionosphere
        InvN(t,:)=0.d0
        InvN(:,t)=0.d0
        InvN(t,t)=0.2d0**2   ! precision is 0.2m
        X(t)=0.d0
    elseif (flag(1:3)=='dds') then  ! unobserved satellite
        InvN(t,:)=0.d0
        InvN(:,t)=0.d0
        X(t)=0.d0
    else 
        write(*,*) 'nargin error in KF_Change.f90'
    end if
    return
end subroutine

! ========================  KF_Gain_one  ======================
!
! PURPOSE:
!           Kalman Filter for parameter constraint.

! INPUTS:
!           InvN               The normal equation, i.e. Pk  x=InvN*U
!           X                    The normal equation  x=InvN*U
!           N                    The size of the normal equation
!           t                      parameter position in InvN
!           y                     constaint value
!           R                    constaint squart variance
!
! OUTPUTS:
!           InvN, X            The new Kalman paramter and Co-Variance matrix
! 
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! =================  End of header  =================
subroutine KF_Gain_one(InvN, X, N, t, y, R)
implicit none
    ! Intent in and out
    real(8) :: InvN(N,N), X(N)
    integer :: N, t
    real(8) :: y, R
    ! Local variables
    real(8) :: Kk(N)
    integer :: i, j

    ! Kalman Filter
    ! Kk=Pk*Gt*[G*Pk*Gt+R]{-1}   X=X+Kk*(y-G*X)  Pk=(I-Kk*G)*Pk
    if (InvN(t,t)==0.d0) InvN(t,t)=1.d6
    Kk=InvN(t,:)/(InvN(t,t)+R**2)  ! Kalman Gain, N*1
    X=X+Kk*(y-X(t))   ! New X, X was 0.d0 here
    do i=1,N
        if (Kk(i)==0.d0) cycle
        do j=i, N
            if (InvN(t,j)==0.d0) cycle
            InvN(j,i)=InvN(j,i)-Kk(i)*InvN(t,j)  ! New InvN
            InvN(i,j)=InvN(j,i)
        end do
    end do

    return

end subroutine


! ========================  KF_Gain  ======================
!
! PURPOSE:
!           Traditional Kalman Filter.

! INPUTS:
!           InvN               The normal equation, i.e. Pk  x=InvN*U, N*N
!           X                    The normal equation  x=InvN*U, N*1
!           N                    The size of the normal equation
!           M                   The number of the error equations
!           G                    The design matrix of the error equations, M*N
!           y                     omc, M*1
!           R                    constaint  variance, M*M
!
! OUTPUTS:
!           InvN, X            The new Kalman paramter and Co-Variance matrix
! 
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! =================  End of header  =================
subroutine KF_Gain(InvN, X, N, M, G, y, R)
implicit none
    ! Intent in and out
    integer :: N, M
    real(8) :: InvN(N,N), X(N)
    real(8) :: G(M,N), y(M), R(M,M)
    ! Local variables
    real(8) :: Kk(N,M), InvNR(M,M)
    
    ! Kalman Filter
    ! Kk=Pk*Gt*[G*Pk*Gt+R]{-1}  X=X+Kk*(y-G*X)  Pk=(I-Kk*G)*Pk
    call InvSqrt(MATMUL(MATMUL(G,InvN),TRANSPOSE(G))+R, M, InvNR)  ! Inverse almost takes no time
    Kk=MATMUL(MATMUL(InvN,TRANSPOSE(G)), InvNR)    ! matrix multiply takes time
    X=X+MATMUL(Kk, y-MATMUL(G,X))
    InvN=InvN-MATMUL(MATMUL(Kk,G),InvN)

    return
end subroutine