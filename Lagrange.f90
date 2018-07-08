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
    integer ::  N,i, j
    real(8) p
    
    Yh=0.0d0
    do i=1,N
        p=1.0d0
        do j=1,N
            if (i /=j) then
                p=p*(Xh-X(j))/(X(i)-X(j))
            end if
        end do
        Yh=Yh+Y(:,i)*p
    end do
    return
    
end subroutine