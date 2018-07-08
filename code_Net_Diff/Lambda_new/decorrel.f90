subroutine decorrel(log, n, Qahat,ahat, Qzhat,Z,L,D,zhat,iZt)
!%DECORREL: Decorrelate a (co)variance matrix of ambiguities
!%
!%     [Qzhat,Z,L,D,zhat] = decorrel (Qahat,ahat)
!%
!% This routine creates a decorrelated Q-matrix, by finding the
!% Z-matrix and performing the corresponding transformation.
!%
!% The method is described in:
!% The routine is based on Fortran routines written by Paul de Jonge (TUD)
!% and on Matlab-routines written by Kai Borre.
!% The resulting Z-matrix can be used as follows:
!% zhat = Zt * ahat; \hat(z) = Z' * \hat(a);
!% Q_\hat(z) = Z' * Q_\hat(a) * Z
!%
!% Input arguments:
!%   Qahat: Variance-covariance matrix of ambiguities (original)
!%   ahat:  Original ambiguities (optional)
!%
!% Output arguments:
!%   Qzhat: Variance-covariance matrix of decorrelated ambiguities
!%   Z:     Z-transformation matrix
!%   L:     L matrix (from LtDL-decomposition of Qzhat)
!%   D:     D matrix (from LtDL-decomposition of Qzhat)
!%   zhat:  Transformed ambiguities (optional)
!%   iZt:   inv(Z')-transformation matrix
!%
!% ----------------------------------------------------------------------
!% Function.: decorrel
!% Date.....: 19-MAY-1999 / modified 12-APRIL-2012
!% Author...: Peter Joosten / Sandra Verhagen
!%            Mathematical Geodesy and Positioning
!%            Delft University of Technology
!% ----------------------------------------------------------------------

! Variables in:
integer :: log, n
real(8) :: Qahat(n,n), ahat(n)
! Variables out:
real(8) :: Qzhat(n,n), Z(n,n), L(n,n), D(n), zhat(n), iZt(n,n)
! Local variables:
integer :: i, j, i1
logical :: sw
integer(1) :: mu
real(8) :: delta, lambda, eta, help, v1(n), Qahat2(n,n), iZtt(n,n)

!%Tests on Inputs ahat and Qahat                           

!%Is the Q-matrix symmetric?
!if ~isequal(Qahat-Qahat'<1E-8,ones(size(Qahat)));
!  error ('Variance-covariance matrix is not symmetric!');
!end;

!%Is the Q-matrix positive-definite?
!if sum(eig(Qahat)>0) ~= size(Qahat,1);
!  error ('Variance-covariance matrix is not positive definite!');
!end;

!% -----------------------
!% --- Initialisations ---
!% -----------------------

do i=1,n
    do j=1,n
        iZt(i,j)=0.d0
    end do
    iZt(i,i)=1.d0
end do
i1   = n - 1
sw   = 1

!% --------------------------
!% --- LtDL decomposition ---
!% --------------------------
Qahat2 = Qahat
call ldldecom(log, n, Qahat2, L, D)

!% ------------------------------------------
!% --- The actual decorrelation procedure ---
!% ------------------------------------------

sw = .true.
do while (sw)

   i  = n   !%loop for column from n to 1
   sw = .false.

   do while (.not. sw  .and. i > 1)

      i = i - 1  !%the ith column
      if (i <= i1) then
      
         do j = i+1, n
            mu = nint(L(j,i))
            if (mu/=0) then !% if mu not equal to 0
               L(j:n,i) = L(j:n,i) - real(mu) * L(j:n,j)
               iZt(:,j) = iZt(:,j) + real(mu) * iZt(:,i)  !%iZt is inv(Zt) matrix 
            end if
         end do

      end if

      delta = D(i) + L(i+1,i)**2 * D(i+1)
      if (delta < D(i+1)) then

         lambda       = D(i+1) * L(i+1,i) / delta
         eta          = D(i) / delta
         D(i)         = eta * D(i+1)
         D(i+1)       = delta

         L(i:i+1,1:i-1) =matmul( reshape((/ -L(i+1,i),  eta, 1.d0,  lambda /),(/2,2/)), L(i:i+1,1:i-1) )
         L(i+1,i)     = lambda

!         % swap rows i and i+1
         do j = i+2, n
            help = L(j,i)              
            L(j,i) = L(j,i+1)              
            L(j,i+1) = help
         end do
        do j = 1, n
            help = iZt(j,i)              
            iZt(j,i) = iZt(j,i+1)              
            iZt(j,i+1) = help
        end do

         i1           = i
         sw           = .true.

      end if

   end do

end do

!% ---------------------------------------------------------------------
!% --- Return the transformed Q-matrix and the transformation-matrix ---
!% --- Return the decorrelated ambiguities, if they were supplied    ---
!% ---------------------------------------------------------------------

!call Inverse(n, iZt, iZtt)
call inv(iZt, iZtt, n)
Z = real(NINT(TRANSPOSE(iZtt)))
Qzhat = MATMUL(MATMUL(TRANSPOSE(Z), Qahat), Z)

!if (nargin == 2 .and. nargout >= 5) then
  zhat = MATMUL(TRANSPOSE(Z), ahat)
!end if

return

end subroutine

!% ----------------------------------------------------------------------
!% End of routine: decorrel
!% ----------------------------------------------------------------------
