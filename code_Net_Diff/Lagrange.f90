! ================  Lagrange  =============
!
! PURPOSE:
!     Lagrange interpolation
!
! INPUTS:
!     X      X, a line vector(N*1)£¬epoch time
!     Y      Y, a N*3 matrix, satellite coordinate, X,Y,Z
!     Xh     Interpolated point
!     N      The number of points
!
! OUTPUTS:
!     Yh     y value in the inperpolated point, 3*1
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!=============== End of Header ===============

subroutine Lagrange(X,Y,Xh,Yh, N)
implicit none
    real(8) :: X(N),Y(3,N), Xh, Yh(3)
    integer ::  N,i, j, k
    real(8) p
    
    Yh=0.0d0; k=0
    do i=1,N
        p=1.0d0
        if (all(Y(:,i)==0.d0) .or. all(Y(:,i)==9999.d0)) then
            cycle
        else
            k=k+1
        end if
        do j=1,N
            if (i /=j .and. all(Y(:,j)/=0.d0) .and. all(Y(:,j)/=9999.d0)) then
                p=p*(Xh-X(j))/(X(i)-X(j))
            end if
        end do
        Yh=Yh+Y(:,i)*p
    end do
    if (k<6) Yh=9999.d0
    return
    
end subroutine