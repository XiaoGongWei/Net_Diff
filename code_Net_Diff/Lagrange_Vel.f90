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
    integer ::  N,i, j, k, m
    real(8) p, q, coe
    
    Yh=0.0d0; m=0
    do i=1,N
        p=0.0d0
        coe=1.d0
        if (all(Y(:,i)==0.d0) .or. all(Y(:,i)==9999.d0)) then
            cycle
        else
            m=m+1
        end if
        do j=1,N
            if (i==j .or. all(Y(:,j)==0.d0) .or. all(Y(:,j)==9999.d0)) cycle
            q=1.d0
            do k=1,N
                if (k==j .or. k==i .or. all(Y(:,k)==0.d0) .or. all(Y(:,k)==9999.d0)) cycle
                q=q*(Xh-X(k))
                !p=p*(Xh-X(j))/(X(i)-X(j))
            end do
            p=p+q
            coe=coe*(X(i)-X(j))
        end do
        Yh=Yh+Y(:,i)*p/coe
    end do
    if (m<6) Yh=9999.d0
    return
    
end subroutine