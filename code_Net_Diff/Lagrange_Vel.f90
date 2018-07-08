! ================ Lagrange_Vel ==============
! 
! PURPOSE:
!     Lagrange velocity interpolation
!
! INPUTS:
!     X      X, a line vector(N*1)£¬epoch time
!     Y      Y, a N*3 matrix, satellite coordinate, X,Y,Z
!     Xh     Interpolated point
!     N      The number of points
!
! OUTPUTS:
!     Yh     velocity in the inperpolated point, 3*1
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!=============== End of Header ===============

subroutine Lagrange_Vel(X,Y,Xh,Yh, N)
implicit none
    real(8) :: X(N),Y(3,N), Xh, Yh(3)
    integer ::  N,i, j, k
    real(8) p, q, coe
    
    Yh=0.0d0
    do i=1,N
        p=0.0d0
        coe=1.d0
        do j=1,N
            if (i==j) cycle
            q=1.d0
            do k=1,N
                if (k==j .or. k==i) cycle
                q=q*(Xh-X(k))
                !p=p*(Xh-X(j))/(X(i)-X(j))
            end do
            p=p+q
            coe=coe*(X(i)-X(j))
        end do
        Yh=Yh+Y(:,i)*p/coe
    end do
    return
    
end subroutine