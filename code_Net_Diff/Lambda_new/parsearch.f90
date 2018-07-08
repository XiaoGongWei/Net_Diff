subroutine parsearch(n, zhat,Qzhat,Z,L,D,P0,ncands, zpar,sqnorm,Qzpar,Zfixedpar,Ps,nfixed,zfixed)
!%
!%   [zpar,sqnorm,Qzpar,Zpar,Ps,nfixed,zfixed]=parsearch(zhat,Qzhat,Z,L,D,P0,ncands)
!%
!%This routine performs an integer bootstrapping procedure for partial 
!%ambiguity resolution (PAR) with user-defined success-rate P0
!%
!%INPUTS:
!%   zhat: decorrelated float ambiguities (zhat, must be a column!)
!%  Qzhat: variance-covariance matrix of decorrelated float ambiguities
!%      Z: Z-matrix from decorrel
!%    L,D: lower-triangular and diagonal matrix from LtDL-decomposition of Qzhat
!%     P0: Minimum required sucess rate [DEFAULT=0.995]
!% ncands: Number of requested integer candidate vectors [DEFAULT=2]
!%
!%OUTPUTS:
!%    zpar: subset of fixed ambiguities (nfixed x ncands) 
!%  sqnorm: squared norms corresponding to fixed subsets
!%   Qzpar: variance-covariance matrix of float ambiguities for the
!%          subset that is fixed
!%    Zpar: Z-matrix corresponding to the fixed subset
!%      Ps: Bootstrapped sucess rate of partial ambiguity resolution
!%  nfixed: The number of  fixed ambiguities
!%  zfixed: [OPTIONAL]: also give the complete 'fixed' ambiguity vector
!%          where the remaining (non-fixed) ambiguities are adjusted 
!%          according to their correlation with the fixed subset
!%
!%
!% NOTE:
!% This PAR algorithm should be applied to decorrelated ambiguities, which
!% can be obtained from the original ahat and its variance-covariance matrix 
!% Qahat as:
!% 
!%     [Qzhat,Z,L,D,zhat] = decorrel(Qahat,ahat);
!%
!% The fixed baseline solution can be obtained with (see documentation):
!%
!%     s      = length(zhat) - nfixed + 1;
!%     Qbz    = Qba * Zpar ;   
!%     bfixed = bhat - Qbz/Qzpar * (zhat(s:end)-zpar(:,1) );
!%     Qbfixed= Qbhat - Qbz/Qzpar * Qbz';
!%
!% Hence, zfixed is not required. LAMBDA, however, does give the solution
!% in terms of the full ambiguity vector.
!%
!%   *************************************************************
!%   *                Author: Sandra Verhagen                    *
!%   *                  Date: 04/APRIL/2012                      *
!%   *                  GNSS Research Centre                     *
!%   *  Dept.of Spatial Science, Curtin University of Technology *
!%   *************************************************************
!
!%============================START PROGRAM==========================%

implicit none
! Intents in:
integer :: n, ncands
real(8) :: zhat(n), Qzhat(n,n),Z(n,n),L(n,n),D(n),P0
! Intents out:
real(8) :: zpar(n),sqnorm(ncands),Qzpar(n,n),Zfixedpar(n,n),Ps,zfixed(n)
integer :: nfixed
! Local variables:
integer :: k, i
real(8) :: QP(n,n), QP2(n,n), temp

!if(nargin <4 )
!    P0=0.995;
!    warning(['user-defined success rate is necessary for PAR method',...
!        'the default value 0.995 is used']);
!end
!if nargin<6
!    ncands = 2;
!end

!n      = size (Qzhat,1);

!% bootstrapped success rate if all ambiguities would be fixed
!Ps = prod ( 2.d0 * normcdf(1./(2.d0*sqrt(D))) -1.d0 )
! In fortran, normcdf(x)=0.5d0*erf(x/sqrt(2.d0))+0.5d0
Ps=1.d0
do i=1,n
    Ps=Ps*(erf(1.d0/(2.d0*sqrt(D(i)))/sqrt(2.d0)))
end do
k = 1
do while (Ps < P0 .and. k < n)   
    k = k + 1
!    % bootstrapped success rate if the last n-k+1 ambiguities would be fixed
!    Ps = prod ( 2 * normcdf(1./(2*sqrt(D(k:n)))) -1 )  
    Ps=1.d0
    do i=k,n
        Ps=Ps*(erf(1.d0/(2.d0*sqrt(D(i)))/sqrt(2.d0)))
    end do
end do

if (Ps > P0) then
    
    !% last n-k+1 ambiguities are fixed to integers with ILS
    call ssearch(n-k+1, zhat(k:n),L(k:n,k:n),D(k:n),ncands, zpar(k:n),sqnorm)

    Qzpar(1:n-k+1,1:n-k+1) = Qzhat(k:n,k:n)
    Zfixedpar(:, 1:n-k+1)  = Z(:,k:n)

!    if (nargout > 6) then

       !% first k-1 ambiguities are adjusted based on correlation with the fixed
       !% ambiguities
!       QP = Qzhat(1:k-1,k:n) / Qzhat(k:n,k:n)
       call inv(Qzhat(k:n,k:n), QP2(k:n,k:n), n-k+1)
       QP(1:k-1,1:n-k+1) =MATMUL(Qzhat(1:k-1,k:n),QP2(k:n,k:n))
       if (k==1) then
            zfixed = zpar
       else
            ! In fortran, we assume that ncands=2, and only choose the best fixed ambiguity
            zfixed(1:k-1) = zhat(1:k-1) - MATMUL(QP(1:k-1,1:n-k+1), (zhat(k:n)-zpar(k:n)))
            zfixed(k:n)=zpar(k:n)
       end if
!    end if
    
    nfixed = n-k+1
else
    zfixed = zhat
    nfixed = 0
end if

return

end subroutine