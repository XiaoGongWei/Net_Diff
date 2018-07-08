subroutine lsearch (ahat,L,D,ncands, afixed,sqnorm)
!%LSEARCH: Integer ambiguity resolution, search
!%
!% This routine finds the integer vector which is closest to a given
!% float vector, in a least squares sence. This is the search-step in
!% integer ambiguity resolution. It is best to perform this search only
!% on ambiguities which have been decorrelated using LAMBDA.
!%
!% Input arguments:
!%    ahat   : Float ambiguities (should be decorrelated for computational
!%             efficiency)
!%    L,D    : LtDL-decomposition of the variance-covariance matrix of the 
!%           : float ambiguities ahat
!%    ncands : Number of requested candidates
!%
!% Output arguments:
!%    afixed : Estimated integers (n x ncands)
!%    sqnorm : Corresponding squared norms (vector, sorted)
!%
!% NOTE: the size of the search ellipsoid is determined with the routine
!% chistart
!%    Chi2 = chistart (D,L,ahat,ncands)
!
!% ----------------------------------------------------------------------
!% File.....: lsearch.m
!% Date.....: 19-MAY-1999
!% Date.....: 28-MAR-2012 modified by Sandra Verhagen
!% Author...: Peter Joosten
!%            Mathematical Geodesy and Positioning
!%            Delft University of Technology
!% ----------------------------------------------------------------------
!
!% ---------------------------------------------------------
!% --- Computes the initial size of the search ellipsoid ---
!% ---------------------------------------------------------
implicit none
! Intents in:
integer :: n, ncands
real(8) :: ahat(n),L(n,n),D(n)
! Intents out:
real(8) :: afixed_best(n),sqnorm(ncands)
! Local variables:
real(8) :: Chi2, dist(n)


call chistart (n, D,L,ahat,ncands,1.d0, Chi2)

!% -------------------------------
!% --- Initializing statements ---
!% -------------------------------

Linv      = inv(L)
Dinv      = 1./D

True      = 1
False     = 0

n         = max(size(ahat))

right     = [zeros(n,1) ; Chi2]
left      = [zeros(n+1,1)]
dq        = [Dinv(2:n)./Dinv(1:n-1) 1/Dinv(n)]

cand_n    = .false.
c_stop    = .false.
endsearch = .false.

ncan      = 0

i         = n + 1
iold      = i

afixed = zeros(n,ncands)
sqnorm = zeros(1,ncands)

!% ----------------------------------
!% --- Start the main search-loop ---
!% ----------------------------------

do while (.not. endsearch)

   i = i - 1

   if (iold <= i) then
      lef(i) = lef(i) + Linv(i+1,i)
   else
      lef(i) = 0.d0
      do j = i+1, n
         lef(i) = lef(i) + Linv(j,i)*distl(j,1)
      end do
   end if
   iold = i
   
   right(i)   = (right(i+1) - left(i+1)) * dq(i)
   reach      = sqrt(right(i))
   delta      = ahat(i) - reach - lef(i)
   distl(i,1) = ceil(delta) - ahat(i)
   
   if (distl(i,1) > reach - lef(i)) then

!%     ----------------------------------------------------
!%     --- There is nothing at this level, so backtrack ---
!%     ----------------------------------------------------

      cand_n = .false.
      c_stop = .false.
      
      do while ((.not. c_stop) .and. (i < n))
      
         i = i + 1;
         if (distl(i) < endd(i)) then
            distl(i) = distl(i) + 1
            left(i)  = (distl(i) + lef(i)) **2
            c_stop   = .true.
            if (i == n) cand_n = .true.
         end if
      
      end do
      
      if (i == n) && (~ cand_n); endsearch = True; end;
      
   else

!%     ----------------------------
!%     --- Set the right border ---
!%     ----------------------------

      endd(i) = reach - lef(i) - 1
      left(i) = (distl(i,1) + lef(i)) **2

   end if

   if (i == 1) then
   
!%     -------------------------------------------------------------------
!%     --- Collect the integer vectors and corresponding               ---
!%     --- squared distances, add to vectors "afixed" and "sqnorm" if: ---
!%     --- * Less then "ncands" candidates found so far                ---
!%     --- * The squared norm is smaller than one of the previous ones ---
!%     -------------------------------------------------------------------

      t       = Chi2 - (right(1)-left(1)) * Dinv(1)
      endd(1) = endd(1) + 1
      
      do while (distl(1) <= endd(1))

         if (ncan < ncands) then
         
            ncan             = ncan + 1
            afixed(:,ncan) = distl + ahat
            sqnorm(ncan)     = t

         else
         
            [maxnorm,ipos] = max(sqnorm)
            if (t < maxnorm) then
               afixed(:,ipos) = distl + ahat
               sqnorm(ipos)     = t
            end if
            
         end if

         t       = t + (2 * (distl(1) + lef(1)) + 1) * Dinv(1)
         distl(1) = distl(1) + 1

      end do
      

!%     -------------------------
!%     --- And backtrack ... ---
!%     -------------------------

      cand_n = .false.
      c_stop = .false.
      
      do while ((.not. c_stop) .and. (i < n))

         i = i + 1

         if (distl(i) < endd(i)) then
            distl(i) = distl(i) + 1
            left(i)  = (distl(i) + lef(i)) **2
            c_stop   = .true.
            if (i == n) cand_n = .true.
         end if
      
      end do

      if ((i == n) .and. (.not. cand_n)) endsearch = .true.
      
   end if

end do

!% ----------------------------------------------------------------------
!% --- Sort the resulting candidates, according to the norm
!% ----------------------------------------------------------------------

tmp    = sortrows ([sqnorm' afixed'])
sqnorm = tmp(:,1)'
afixed = tmp(:,2:n+1)'
afixed = round(afixed)

!% ------------------------
!% --- Check for errors ---
!% ------------------------

if (ncan < ncands) then
!    error ('Not enough candidates were found!!')
end if

!% ----------------------------------------------------------------------
!% End of routine: lsearch
!% ----------------------------------------------------------------------
end subroutine
