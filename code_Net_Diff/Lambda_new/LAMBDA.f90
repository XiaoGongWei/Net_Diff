
subroutine LAMBDA(log, n,ahat,Qahat,method,afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu,dz)
!%
!%   [afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu]=LAMBDA(ahat,Qahat,method,varargin)
!%
!% This is the main routin of the LAMBDA software package. By default the ILS
!% method will be used for integer estimation based on the provided float
!% ambiguity vector ahat and associated variance-covariance matrix Qahat.
!% However, the user may also select other methods: integer rounding,
!% bootstrapping or Partial Ambiguity Resolution (PAR). Furthermore, there is
!% the option to apply the Ratio Test to decide on acceptance of the fixed
!% solution.
!%
!% NOTE 1: LAMBDA always first applies a decorrelation before the integer
!% estimation (for ILS this is required to guarantee an efficient search, for
!% rounding and bootstrapping it is required in order to get higher success
!% rates).
!%
!% NOTE 2: Two different ILS search strategies are implemented (see
!% documentation). The solution will be the same with both methods, but the
!% search-and-shrink method is known to be faster.
!%
!% INPUTS:
!%
!%     ahat: Float ambiguities (must be a column!)
!%    Qahat: Variance/covariance matrix of ambiguities 
!%   method: 1: ILS method based on search-and-shrink [DEFAULT]
!%           2: ILS method based enumeration in search
!%           3: integer rounding method
!%           4: integer bootstrapping method
!%           5: PAR with the input P0 of user-defined success rate
!%           6: ILS method with Ratio Test (uses search-and shrink)
!% varargin: Optional input arguments, see below        
!%
!% OUTPUTS:
!%
!%   afixed: Array of size (n x ncands) with the estimated integer
!%           candidates, sorted according to the corresponding squared norms, 
!%           best candidate first. 
!%           For integer rounding and bootstrapping: ncands = 1
!%   sqnorm: Distance between integer candidate and float ambiguity vectors 
!%           in the metric of the variance-covariance matrix Qahat.
!%           Only available for ILS.
!%       Ps: Bootstrapped success rate. 
!%           If ILS/PAR is used, Ps is its lower bound;
!%           If rounding is used, Ps is its upper bound;
!%           If bootstrapping is used, Ps is the exact success rate.
!%    Qzhat: Variance-covariance matrix of decorrelated ambiguities
!%           (corresponding to fixed subset in case of PAR).
!%        Z: Transformation matrix with 
!%           - dimension (n x n) for methods 1-4, 6
!%           - dimension (n x nfixed) for method 5 (PAR).
!%   nfixed: Number of fixed ambiguities 
!%           - with methods 1 to 4: will always be equal to n
!%           - with method 5 (PAR): will be equal to the number of 
!%             fixed decorrelated ambiguities
!%           - with method 6 (ILS + Ratio test): will be equal to n if fixed
!%             solution is accepted, and 0 otherwise
!%       mu: Threshold value used for Ratio Test
!%
!%
!% [AFIXED,...] = LAMBDA(ahat,Qahat,method,'option',value,...) can be used 
!% to specify various options. Possible options are
!%
!%   'ncands',value   Number of requested integer candidate vectors 
!%                    (only used with ILS/PAR, DEFAULT = 2) 
!%   'P0',value       - with method 5 (PAR): Minimum required success rate
!%                      [DEFAULT=0.995]
!%                    - with method 6 (ILS + Ratio test): Fixed failure rate 
!%                      (available options: 0.01 or 0.001) [DEFAULT=0.001]
!%   'MU',value       Fixed threshold value for Ratio Test 
!%                    (value must be between 0 and 1)
!%
!% The Ratio Test will only be applied if method 6 is selected. By default
!% the Fixed Failure rate Ratio Test will be applied. The threshold value MU 
!% is then determined such that the failure rate will not exceed the 
!% specified P0 (DEFAULT 0.001).
!% If MU is specified, the value for P0 is ignored.
!% NOTE: the Ratio Test used here is
!%
!%         Accept afixed iff: sqnorm(1)/sqnorm(2) <= MU
!%
!% Hence, the squared norm of the best (ILS) integer solution is in the
!% numerator. In literature often the reciprocal is used; the corresponding
!% critical value is then c = 1/MU.
!%
!%--------------------------------------------------------------------------
!% Release date  : 1-SEPT-2012                                          
!% Authors       : Bofeng LI and Sandra VERHAGEN                                             
!%  
!% GNSS Research Centre, Curtin University
!% Mathematical Geodesy and Positioning, Delft University of Technology                              
!%--------------------------------------------------------------------------
!%
!% REFERENCES: 
!%  1. LAMBDA Software Package: Matlab implementation, Version 3.0.
!%     Documentation provided with this software package.
!%  2. Teunissen P (1993) Least-squares estimation of the integer GPS
!%     ambiguities. In: Invited lecture, section IV theory and methodology,
!%     IAG General Meeting, Beijing, China
!%  3. Teunissen P (1995) The least-squares ambiguity decorrelation
!%     adjustment: a method for fast GPS ambiguity estitmation. J Geod
!%     70:651-7
!%  4. De Jonge P, Tiberius C (1996) The LAMBDA method of intger ambiguity 
!%     estimation:implementation aspects.
!%  5. Chang X ,Yang X, Zhou T (2005) MLAMBDA: a modified LAMBDA method for
!%     integer least-squares estimation

! Changed to fortran source code by Yize Zhang @ 2018/03/09 TUMSAT


!%=============================START PROGRAM===============================%
! Intents in:
integer :: log
integer :: n, method
real(8) :: ahat(n),  Qahat(n,n)
! Intents out:
real(8) :: afixed(n), sqnorm(2), Ps, Qzhat(n,n), Z(n,n), mu, dz(n)
integer :: nfixed
! Local variables:
integer :: ncands
real(8) :: P0, FFRT
real(8) :: zfixed(n), zpar(n)
real(8) :: incr(n)
real(8) :: L(n,n), D(n), zhat(n), iZt(n,n)
real(8) :: Qzhat2(n,n), Z2(n,n)

!if(nargin<2) then
!    error(['Not enough inputs: float solution',  ...
!        'and its variance-covariance matrix must be specified']);
!end

!%By default use ILS with search-and-shrink
!if (nargin<3)   method = 1

!%If ILS/PAR method is not used, only a unique solution is available
if (method==3 .or. method==4) then
    ncands = 1
else
    ncands = 2
end if

!%For PAR method, the user-defined success rate is necessary. 
!%Default value is 0.995
if (method==5) then
    P0     = 0.995d0
end if

!%Default values for Ratio Test
if (method==6) then
    P0     = 0.001d0
    FFRT   = 1.d0
else
    mu     = 1.d0
end if

!% Set the values of the optional input arguments in varargin
!i = 0;
!while i < length(varargin)
!   switch upper(varargin{i+1})
!      case 'P0'
!         P0 = varargin{i+2};         
!         if(method == 6 && (~ismember(P0,[0.01,0.001])))
!            error('Fixed failure rate must be either 0.01 or 0.001.');
!         elseif(method == 5 && P0>1.0)
!            error('User-defined success rate P0 cannot be larger than 1');
!         end
!         i=i+2;
!      case {'MU'}
!         mu   = varargin{i+2};
!         if(method == 6 && mu>1)
!            error('MU must be between 0 and 1');
!         end
!         if ~isempty(mu), FFRT = 0; end
!         i=i+2;
!      case {'NCANDS'}
!         if ismember(method,[1,2,5])
!            ncands = varargin{i+2};
!         end
!         i=i+2;
!      otherwise 
!         error('Unrecognised input option.') 
!   end
!end


zfixed = 0.d0
nfixed = n
!sqnorm = []

!%-----------------------------------------------------------------
!%Tests on Inputs ahat and Qahat                           

!%Is the Q-matrix symmetric?
!if ~isequal(Qahat-Qahat'<1E-8,ones(size(Qahat)));
!  error ('Variance-covariance matrix is not symmetric!');
!end;

!%Is the Q-matrix positive-definite?
!if sum(eig(Qahat)>0) ~= size(Qahat,1);
!  error ('Variance-covariance matrix is not positive definite!');
!end;

!%Do the Q-matrix and ambiguity-vector have identical dimensions?
!if length(ahat) ~= size(Qahat,1);
!  error (['Variance-covariance matrix and vector of ambiguities do', ...
!      'not have identical dimensions!']);
!end;

!%Is the ambiguity vector a column?  
!if size(ahat,2) ~= 1;
!  error ('Ambiguity vector should be a column vector');
!end;


!%remove integer numbers from float solution, so that all values are between
!%-1 and 1 (for computational convenience only)
incr = ahat - mod(ahat,1.d0)
ahat = mod(ahat,1.d0)

!%Compute Z matrix based on the decomposition  Q=L^T*D*L; The transformed
!%float solution: \hat{a} = Z^T *ahat, Qzhat = Z^T * Qahat * Z

call decorrel(log, n,Qahat,ahat, Qzhat,Z,L,D,zhat,iZt)

!%Compute the bootstrapped success rate
!if  (method == 6) then
!    Ps = prod(2*normcdf(0.5./sqrt(D))-1)
!end if

select case (method)
    case (1)
        !%ILS with shrinking search
        call ssearch(n, zhat,L,D,ncands, zfixed,sqnorm)
    case (2)
        !%ILS with enumeration search
!        call lsearch(n, zhat,L,D,ncands, zfixed,sqnorm)
    case (3)
        !%Integer rounding
        zfixed = nint(zhat)
    case (4)
        !%Integer bootstraping
!        call bootstrap(zhat,L, zfixed)
    case (5)
        call parsearch(n,zhat,Qzhat,Z,L,D,P0,ncands,  zpar,sqnorm,Qzhat2,Z2,Ps,nfixed,zfixed)
        Qzhat=Qzhat2
        Z=Z2
        if (nfixed == 0) then
            ncands = 1
        end if
    case (6)
        !% ILS with Ratio Test        
!        call ssearch(zhat,L,D,ncands, zfixed,sqnorm)
        if (FFRT==1) then
            if (1.d0-Ps > P0) then
!                mu = ratioinv(P0,1.d0-Ps,n)
            else !% if ILS failure rate smaller than P0: always accept
                mu = 1.d0
            end if
        end if
        !% Perform Ratio Test
        if (sqnorm(1)/sqnorm(2) > mu) then  !% rejection: keep float solution
            zfixed = zhat
            nfixed = 0
            ncands = 1
        end if
end select

!%Perform the back-transformation and add the increments
afixed = matmul(iZt, zfixed)
afixed = afixed + incr
dz=zhat-zfixed

return

end subroutine