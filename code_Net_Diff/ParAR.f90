! ==================== ParAR ====================
! PURPOSE:
!             partial ambiguity resolution
!
! INPUTS:
!              Q                        covariance matrix
!              amb                   ambiguity

! OUTPUTS:
!             Q2
!             amb2
! 
! WRITTEN BY: Yize Zhang
! ================================================

subroutine ParAR(Q, amb, PRNS, PRNPRN, P, Last_PAR, flag_partial, Q2, amb2, iPOS, iPOS2)
use MOD_constant
use MOD_VAR
use MOD_CycleSlip
use MOD_FileID
implicit none
    ! Intents in
    integer :: n, PRNS
    integer :: PRNPRN(PRNS)
    integer(1) :: flag_partial
    integer(1) :: Last_PAR(2*MaxPRN)
    ! Intents out:
    ! Local variables
    integer :: i, j, k, m, l, PRN
    integer :: npar, npar2, nfixed
    real(8) :: Ps, Qzhat(2*MaxPRN, 2*maxPRN), Z(2*MaxPRN, 2*maxPRN), mu
    real(8) :: Q(2*MaxPRN, 2*maxPRN), Q2(2*MaxPRN, 2*maxPRN), P(MaxPRN, maxPRN)
    real(8) :: disall(2), ratio, ratio2, amb(2*MaxPRN), amb2(2*MaxPRN), amb3(2*MaxPRN), dz(2*MaxPRN)
    integer :: iPOS(2*MaxPRN), iPOS2(2*MaxPRN), iPOS3(2*MaxPRN), minL, minL2, minLL(MaxPRN),minLL2(MaxPRN)
    real(8) :: minP, maxQ


    Q2=Q; amb2=amb; npar2=npar; iPOS2=iPOS; amb3=amb; iPOS3=iPOS  !; P=NEQ%P
    flag_partial=0; ratio2=0.d0; k=0; m=1; l=1; minLL=0; minLL2=0

    100 if (npar>1) then
        call LAMBDA(lambdaID, npar, amb(1:npar),Q(1:npar, 1:npar)/9.d0,1,amb(1:npar),disall,Ps,Qzhat(1:npar, 1:npar),Z(1:npar, 1:npar),nfixed,mu,dz(1:npar))
        if (nfixed==0) then
            ratio=0.d0
        else
            ratio=dabs(disall(2)/disall(1))
        end if
    end if
    if (ratio>minratio) then  ! ambiguity fix success
        return
    elseif (partial_AR) then ! partial ambigulty fixing
        minP=100.d0; maxQ=0.d0
        minL=0; minL2=0
        amb3=amb2
        iPOS=iPOS2
        do i=1,npar2  ! Start from the maximum variance satellite
            PRN=iPOS3(i)
            if (PRN==0) cycle
            if (Last_PAR(PRN)==1) then ! First find the previous PRN of partial AR
                minL=i  
                do j=i+1,npar2  ! Find the same partial AR PRN in L2
                    if (PRN+MaxPRN==iPOS3(j) .and. Last_PAR(iPOS3(j))==1) then
                        minL2=j
                        exit
                    end if
                end do
                exit
            end if
            !!   If not, start from the minmum satellite elevation
            if (PRN>maxPRN) PRN=PRN-maxPRN
            do j=1,PRNS
                if (PRNPRN(j)==PRN) then
                    exit
                end if
            end do
            if ( (PRN/=0) .and. (P(j,j)<minP-1.d-5) ) then
                minP=P(j,j)
                minL=i
            elseif ( (PRN/=0) .and. abs(P(j,j)-minP)<1.d-5 ) then
                minL2=i   ! dual frequency
            end if
        end do
        if (minL/=0) then
            k=k+1
            minLL(k)=minL
            iPOS3(minL)=0
            amb3(minL)=0.d0
            if (minL2/=0) then
                iPOS3(minL2)=0
                amb3(minL2)=0.d0
                minLL2(k)=minL2
            end if
            call LAMBDA_Prepare(Q2(1:npar2, 1:npar2), amb3(1:npar2), npar2, Q(1:npar2, 1:npar2), amb(1:npar2), npar, iPOS(1:npar2))
            flag_partial=1
            if (npar>=3) goto 100
        elseif ( (parARnum==2) ) then
         ! Which means ratio not less than minratio for only one satellite partial AR, then we try two.
         ! However, this is only for that last epoch is partial AR
            amb3=amb2
            iPOS=iPOS2
            l=l+1
            if (l>k) then
                m=m+1  ! The location of one satellite in minLL and minLL2
                l=m+1  ! The location of the other satellite in minLL and minLL2
            end if
            if (m<k) then
                if (minLL(m)/=0) then
                    amb3(minLL(m))=0.d0
                end if
                if (minLL2(m)/=0) then
                    amb3(minLL2(m))=0.d0
                end if
                if (minLL(l)/=0) then
                    amb3(minLL(l))=0.d0
                end if
                if (minLL2(l)/=0) then
                    amb3(minLL2(l))=0.d0
                end if

                call LAMBDA_Prepare(Q2(1:npar2, 1:npar2), amb3(1:npar2), npar2, Q(1:npar2, 1:npar2), amb(1:npar2), npar, iPOS(1:npar2))
                flag_partial=1
                if (npar>=3) goto 100
            elseif (m==k) then  ! If two satellites still doesn't fix, use the cycle slip one together with last partial AR
                amb3=amb2
                iPOS=iPOS2
                do i=1,npar2
                    PRN=iPOS2(i)
                    if (par_PRN(PRN)==1) amb3(i)=0.d0
                    if (PRN>MaxPRN) then
                        if (CycleSlip(1)%Slip(PRN-MaxPRN)==1 .or. CycleSlip(2)%Slip(PRN-MaxPRN)==1) then
                            amb3(i)=0.d0
                        end if
                    elseif (PRN<=MaxPRN .and. (CycleSlip(1)%Slip(PRN)==1 .or. CycleSlip(2)%Slip(PRN)==1)) then
                        amb3(i)=0.d0
                    end if
                end do
                call LAMBDA_Prepare(Q2(1:npar2, 1:npar2), amb3(1:npar2), npar2, Q(1:npar2, 1:npar2), amb(1:npar2), npar, iPOS(1:npar2))
                flag_partial=1
                if (npar>=3) goto 100
            end if
        end if  !  if (minL/=0) then
        
    end if

end subroutine