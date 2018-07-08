!   ============== Polyfit =================
!
! PURPOSE:
!             This is the function for quadratic polyfit using least square
!             Use for doppler velocity prediction in RTK
! 
! INPUTS:
!           X          X vector(real)
!           Y          Y vector(real)
!           N         N points(integer)
!          Flag      RTK solution flag, 1: full fixed; 2: partial fixed; 3: float
!
! OUTPUTS:
!           Xp        interpolated point(real)
!           Yp        interpolated Y(real)
!           j           last epoch which has velocity
!           
! WRITTEN BY: Yize Zhang, zhyize@163.com
! 
! =============END OF HEADER============

subroutine Polyfit(X,Y, N, Flag, Xp, Yp, j)
implicit none
    ! Intents:
    real(8) :: X(N), Y(3,N)
    integer :: N
    integer(1) :: Flag(N)
    real(8) :: Xp, Yp(3)
    integer :: j
    ! Local variables
    integer :: i
    real(8) :: A(N,3), L(3,N), P, Nbb(3,3), InvN(3,3), U(3), dx(3)

    j=0
    do i=1,N
        if (all(Y(1:3,i)==0.d0) .or. X(i)-Xp<-6.d0) cycle  ! If no velocity or time difference > 5s
        j=j+1
        if (Flag(i)==1) then
            P=1.d0
        elseif (Flag(i)==2) then
            P=1.d0/2.d0**2
        elseif (Flag(i)==3) then
            P=1.d0/5.d0**2
        end if
        A(j,1:3)=(/1.d0, X(i), X(i)*X(i)/)*P
        L(:,j)=Y(:, i)*P
    end do

    if (j<=3) then
        Yp=0.d0  ! If less than 3 points
        return
    end if

    Nbb=matmul(transpose(A(1:j, :)), A(1:j, :))
    call InvSqrt(Nbb, 3, InvN)
    dx=matmul(InvN, matmul(transpose(A(1:j, :)), L(1,1:j)))
    Yp(1)=dx(1)+dx(2)*Xp+dx(3)*Xp**2
    dx=matmul(InvN, matmul(transpose(A(1:j, :)), L(2,1:j)))
    Yp(2)=dx(1)+dx(2)*Xp+dx(3)*Xp**2
    dx=matmul(InvN, matmul(transpose(A(1:j, :)), L(3,1:j)))
    Yp(3)=dx(1)+dx(2)*Xp+dx(3)*Xp**2

     return

end subroutine