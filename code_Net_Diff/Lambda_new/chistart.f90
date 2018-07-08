subroutine  chistart (D,L,ahat,ncands,factor, Chi2)
!%CHISTART: Computes the initial size of the search ellipsoid
!%
!%    Chi2 = chistart (D,L,ahat,ncands,factor)
!%
!% This routine computes or approximates the initial size of the search
!% ellipsoid. If the requested number of candidates is not more than the
!% dimension + 1, this is done by computing the squared distances of partially
!% conditionally rounded float vectors to the float vector in the metric of the
!% covariance matrix. Otherwise an approximation is used.
!%
!% Input arguments
!%    L,D   : LtDL-decomposition of the variance-covariance matrix of
!%            the float ambiguities (preferably decorrelated)
!%    ahat  : float ambiguites (preferably decorrelated)
!%    ncands: Requested number of candidates (default = 2)
!%    factor: Multiplication factor for the volume of the resulting
!%            search ellipsoid (default = 1.5)
!%
!% Output arguments:
!%    Chi2  : Size of the search ellipsoid
!
!% ----------------------------------------------------------------------
!% File.....: chistart.m
!% Date.....: 19-MAY-1999
!% Modified.: 05-MAR-2001, by P. Joosten
!% Author...: Peter Joosten
!%            Mathematical Geodesy and Positioning
!%            Delft University of Technology
!% ----------------------------------------------------------------------

!% ------------------
!% --- Initialize ---
!% ------------------

implicit none
! Intents in:
integer :: n
real(8) :: D(n), L(n,n), ahat(n), ncands,factor
! Intents out:
real(8) :: Chi2(ncands)
! Local variables:
integer :: i, j, k
real(8) :: iQ(n,n), afloat(n), afixed(n), dw, tmp

!if nargin < 4; ncands = 2  ; end;
!if nargin < 5; factor = 1.5; end;

!% ----------------------------------------------------------------------
!% --- Computation depends on the number of candidates to be computed ---
!% ----------------------------------------------------------------------

if (ncands <= n+1) then

!  % --------------------------------------------------------
!  % --- Computation based on the bootstrapping estimator ---
!  % --------------------------------------------------------

  Chi = []
  iQ  = (L \ diag(1./D)) / L'
  
  do k = n, 0, -1

    afloat = ahat
    afixed = ahat
  
    do i = n, 1, -1
        
        dw = 0
        do j = n, i, -1
            dw = dw + L(j,i) * (afloat(j) - afixed(j))
        end do
        
        afloat(i) = afloat(i) - dw
        if (i /= k) then
            afixed(i) = real(NINT(afloat(i)))
        else
            tmp   = real(NINT(afloat(i)))
            if (afloat(i)-tmp>0.d0) then
                afixed(i) = tmp + 1.d0
            elseif (afloat(i)-tmp<0.d0) then
                afixed(i) = tmp - 1.d0
            else
                afixed(i) = tmp
            end if   
        end if
        
    end do
    Chi = [Chi (ahat-afixed)' * iQ * (ahat-afixed)]

  end do

!  % ---------------------------------------------------------------
!  % --- Sort the results, and return the appropriate number     ---
!  % --- Add an "eps", to make sure there is no boundary problem ---
!  % ---------------------------------------------------------------

  Chi  = sort(Chi)
  Chi2 = Chi(ncands) + 1.d-6
   
else

!  % -----------------------------------------------------
!  % An approximation for the squared norm is computed ---
!  % ----------------------------------------------------- 
  
  Vn   = (2/n) * (pi ^ (n/2) / gamma(n/2))
  Chi2 = factor * (ncands / sqrt((prod(D)) * Vn)) ^ (2/n)

end if

!% ----------------------------------------------------------------------
!% End of routine: chistart
!% ----------------------------------------------------------------------

end subroutine
