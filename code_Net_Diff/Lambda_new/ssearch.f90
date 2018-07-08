subroutine ssearch(n, ahat,L,D,ncands, afixed_best,sqnorm)
!%
!%        [afixed,sqnorm]=ssearch(ahat,L,D,ncands)
!%
!%   Integer ambiguity vector search by employing the search-and-shrink
!%   technique.
!%
!%INPUTS:
!%
!%    ahat   : Float ambiguities (should be decorrelated for computational
!%             efficiency)
!%    L,D    : LtDL-decomposition of the variance-covariance matrix of the 
!%           : float ambiguities ahat
!%    ncands : Number of requested candidates
!%
!%OUTPUTS:
!%
!% afixed: estimated integers (n x ncands )
!% sqnorm: corresponding squared norms (n-vector, ascendantly sorted )
!%
!%------------------------------------------------------------------|
!% DATE    : 02-SEPT-2010                                           |
!% Author  : Bofeng LI                                              |
!%           GNSS Research Center, Department of Spatial Sciences   |
!%           Curtin University of Technology                        |
! %e-mail  : bofeng.li@curtin.edu.au                                | 
!%------------------------------------------------------------------|
!%
!% REFERENCES:                                                      
!%  1. de Jonge P, Tiberius C(1996) The LAMBDA method of intger ambiguity 
!%     estimation:implementation aspects.
!%  2. Chang X ,Yang X, Zhou T(2005) MLAMBDA: a modified LAMBDA method for
!%     integer least-squares estimation
!%  3. Teunissen P(1993) Least-squares estimation of the integer GPS
!%     ambiguities. In: Invited lecture, section IV theory and methodology,
!%     IAG General Meeting, Beijing, China
!%  4. Teunissen P(1995) The least-squares ambiguity decorrelation
!%     adjustment: a method for fast GPS ambiguity estitmation. J Geod
!%     70:65?2
!
!%===========================START PROGRAM===============================%
implicit none
! Intents in:
integer :: n, ncands
real(8) :: ahat(n),L(n,n),D(n)
! Intents out:
real(8) :: afixed_best(n),sqnorm(ncands)
! Local variables:
real(8) :: Chi2, dist(n)
logical :: endsearch
integer :: count, imax, k
real(8) :: afixed(n,ncands)
real(8) :: acond(n), zcond(n), left, step(n), S(n,n), newdist, temp


!if size(ahat,2)~=1 || length(ahat)~=length(D)
!    error('Float ambiguity vector must be a column vector with same dimension as D');
!end

!%Initializing outputs
afixed = 0.d0
sqnorm = 0.d0

!%initializing the variables for searching
Chi2     = 1.0e+18        !%start search with an infinite chi^2
dist(n)  = 0              !%dist(k)=sum_{j=k+1}^{n}(a_j-acond_j)^2/d_j 
endsearch= .false.
count    = 0              !%the number of candidates

acond=0.d0
zcond=0.d0
step=0.d0
acond(n) = ahat(n)
zcond(n) = real(NINT(acond(n)))
left     = acond(n) - zcond(n)
if (left>0.d0) then
    step(n)  = 1.d0
elseif (left<0.d0) then
    step(n)  = -1.d0
else
    step(n)  = 0.d0
end if

!%-----------------------------------------------------------------------%
!%For a very occasional case when the value of float solution ahat(n)==0, we
!%compusively give a positive step to continue. This case can
!%actually never happen in reality, but only when the exact integer value
!%is specified for ahat. 
if (step(n)==0) then
    step(n) = 1
end if
!%------------------------------------------------------------------------%

imax     = ncands         !%initially, the maximum F(z) is at ncands

S(1:n, 1:n) = 0.d0           !%used to compute conditional ambiguities

k = n

!%Start the main search-loop
do while (.not. endsearch)
    !%newdist=sum_{j=k}^{n}(a_j-acond_j)^2/d_j=dist(k)+(a_k-acond_k)^2/d_k
    
    newdist = dist(k) + left**2/D(k)
    
    if (newdist < Chi2) then
        
        if (k/=1)  then       !%Case 1: move down
            k = k - 1
            dist(k)  = newdist
            S(k,1:k) = S(k+1,1:k) +(zcond(k+1)-acond(k+1))*L(k+1,1:k)
            
            acond(k) = ahat(k) + S(k, k)
            zcond(k) = real(NINT(acond(k)))
            left     = acond(k) - zcond(k)
            if (left>0.d0) then
                step(k) = 1.d0
            elseif (left<0.d0) then
                step(k) = -1.d0
            else
                step(k) = 0.d0
            end if
            
            !%-----------------------------------------------------------------------%
            !%For a very occasional case when the value of float solution ahat(n)==0, we
            !%compusively give a positive step to continue. This case can
            !%actually never happen in reality, but only when the exact integer value
            !%is specified for ahat. 
            if (step(k)==0)  step(k) = 1
            !%------------------------------------------------------------------------%
        else
            
            !%Case 2: store the found candidate and try next valid integer
            if (count < ncands - 1) then
                !%store the first ncands-1 initial points as candidates
                
                count = count + 1
                afixed(:, count) = zcond(1:n)
                sqnorm(count) = newdist          !%store F(zcond)
           
            else
                
                afixed(:,imax) = zcond(1:n)
                sqnorm(imax)   = newdist
                Chi2  = maxval(sqnorm)
                imax = maxloc(sqnorm,dim=1)
                
            end if
            
            zcond(1) = zcond(1) + step(1)     !%next valid integer
            left     = acond(1) - zcond(1)
            if (step(1)>0.d0) then
                step(1)  =-step(1)  - 1.d0
            else
                step(1)  =-step(1) + 1.d0
            end if
        end if
        
    else
        !%Case 3: exit or move up
        if (k == n) then
            endsearch = .true.
        else 
            k        = k + 1         !%move up
            zcond(k) = zcond(k) + step(k)  !%next valid integer
            left     = acond(k) - zcond(k)
            if (step(k)>0.d0) then
                step(k)  =-step(k)  - 1.d0
            else
                step(k)  =-step(k) + 1.d0
            end if
        end if
    end if
end do

![sqnorm, order]=sort(sqnorm)
!afixed = afixed(:,order)
! In fortran, we assume that ncands=2, and only choose the best fixed ambiguity
if (sqnorm(1)>sqnorm(2)) then
    afixed_best=afixed(:,2)
    temp=sqnorm(1)
    sqnorm(1)=sqnorm(2)
    sqnorm(2)=temp
else
    afixed_best=afixed(:,1)
end if

return

end subroutine
