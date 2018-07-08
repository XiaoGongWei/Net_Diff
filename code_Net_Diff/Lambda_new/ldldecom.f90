subroutine ldldecom(log, n, Qahat, L, D)
!%LDLDECOM: Find LtDL-decompostion of Qahat-matrix
!%
!%           [L,D] = ldldecom(Qahat)
!%
!% This routine finds the LtDL decomposition of a given variance/
!% covariance matrix.
!%
!% Input arguments:
!%    Qahat: Symmetric n by n matrix to be factored
!%
!% Output arguments:
!%    L:     n by n factor matrix (strict lower triangular)
!%    D:     Diagonal n-vector

!% ----------------------------------------------------------------------
!% File.....: ldldecom
!% Date.....: 19-MAY-1999
!% Author...: Peter Joosten
!%            Mathematical Geodesy and Positioning
!%            Delft University of Technology
!% ----------------------------------------------------------------------

! Intents in:
integer :: log, n
real(8) :: Qahat(n,n)
! Intents out:
real(8) :: L(n,n), D(n)
! Local variables:
integer :: i, j
real(8) :: sum

do i = n, 1, -1

   D(i) = Qahat(i,i)
   L(i,1:i) = Qahat(i,1:i)/sqrt(Qahat(i,i))
   
   do j = 1, i-1
      Qahat(j,1:j) = Qahat(j,1:j)-L(i,1:j)*L(i,j)
   end do
   
   L(i,1:i) = L(i,1:i)/L(i,i)

end do

do i= 1, n
    sum=sum+D(i)
end do

if (sum < 1D-10) then

  write (log,'(4X,A)') 'Matrix on input is not positive definite!'

end if

end subroutine
!% ----------------------------------------------------------------------
!% End of routine: ldldecom
!% ----------------------------------------------------------------------
