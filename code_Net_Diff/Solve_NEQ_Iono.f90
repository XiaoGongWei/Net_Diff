! ============= Sovle_NEQ_Iono===================== 
! PURPOSE:
!            Solve the normal equation and get the coordinate.
!        Steps:
!            1) Solve the normal equation and get the wide lane
!     and W4 combination ambiguity (float);
!            2) Fix the WL and W4 ambiguity  using LAMBDA method;
!            3) Change WL ambiguity  to pseudo-range and add to EPO_NEQ;
!            4) Solve L1 and L2 ambiguity(float) and ionosphere delay£»
!            5) Fix L1 and L2 ambiguity using LAMBDA method;
!            6) Estimate the coordinate by fixing the ambiguity of 
!      L1 and L2.
!
! INPUTS:
!         NEQ                         normal equation structure, used in Step 1
!         Epo_NEQ                 current epoch normal equation structure, 
!                                              used in Step 4
!
! OUTPUT:
!         Coor                        Coordinate

! WRITTEN BY: Yize Zhang
! =============================================

subroutine Solve_NEQ_Iono(NEQ, Epo_NEQ, Coor, Flag_Sln)
use MOD_FileID
use MOD_NEQ
use MOD_Epo_NEQ
use MOD_constant
use MOD_VAR
use MOD_STA
use MOD_GLO_Fre
use MOD_NEQ_DP
implicit none
    type(type_NEQ) :: NEQ
    type(type_Epo_NEQ) :: Epo_NEQ
    integer :: i, j, k, m, PRN, N, l
    logical :: AD_flag, flag_par_PRN
    integer(1) :: Flag_Sln
    real(8) :: Coor(4)
    real(8) :: maxV
    integer :: maxL
    ! for lambda
    integer :: npar, npar2, nfixed
    real(8) :: Ps, Qzhat(2*MaxPRN, 2*maxPRN), Z(2*MaxPRN, 2*maxPRN), mu
    real(8) :: Q(2*MaxPRN, 2*maxPRN), Q2(2*MaxPRN, 2*maxPRN), P(MaxPRN, maxPRN)
    real(8) :: disall(2), ratio, ratio2, amb(2*MaxPRN), amb2(2*MaxPRN), amb3(2*MaxPRN), dz(2*MAxPRN)
    real(8) :: dx(2*MaxPRN), temp_Nbb(ParaNum+MaxPRN,ParaNum+MaxPRN), temp_InvN(ParaNum+MaxPRN,ParaNum+MaxPRN)
    real(8) :: temp_U(ParaNum+MaxPRN), temp_dx(ParaNum+MaxPRN), temp_Nbb2(2*MaxPRN,ParaNum+MaxPRN)
    real(8) :: temp_QR(MaxPRN+ParaNum,2*MaxPRN), temp_QQ(2*MaxPRN, 2*MaxPRN)
    integer :: iPOS(2*MaxPRN), iPOS2(2*MaxPRN), iPOS3(2*MaxPRN), minL, minL2, minLL(MaxPRN),minLL2(MaxPRN)
    real(8) :: minP, maxQ
    integer(1) :: flag_partial, flag_fixed, freq, sys, par_PRNS

    !   Velocity estimation of Doppler Data
    if (Combination(3)) then
        Ad_Flag=.true.
        N=NEQ_DP%PRNS   ! Obs number in A
        if (N>=6) then
            NEQ_DP%Nbb=matmul(transpose(NEQ_DP%A(1:N,:)),NEQ_DP%A(1:N,:))
            NEQ_DP%U=matmul(transpose(NEQ_DP%A(1:N,:)),NEQ_DP%L(1:N))
            do while (Ad_Flag)
                Ad_Flag=.false.
                call InvSqrt(NEQ_DP%Nbb, 4 ,NEQ_DP%InvN)
                NEQ_DP%dx=MATMUL(NEQ_DP%InvN, NEQ_DP%U)
                if (any(isnan(NEQ_DP%dx))) then
                    write(*,*)  "-------ERROR-------: NEQ_DP%dx=nan, please check"
                    write(LogID, '(5X,A40)')  "-------ERROR-------: NEQ_DP%dx=nan, please check"
                    stop
                end if
                NEQ_DP%V(1:N)=MATMUL(NEQ_DP%A(1:N,:), NEQ_DP%dx)-NEQ_DP%L(1:N)
                NEQ_DP%sigma0=dsqrt(DOT_PRODUCT(NEQ_DP%V(1:N), NEQ_DP%V(1:N))/(N-4) )
                
                NEQ_DP%maxV=maxval(dabs(NEQ_DP%V(1:N)))
                NEQ_DP%maxL=maxloc(dabs(NEQ_DP%V(1:N)),dim=1)
                maxL=NEQ_DP%maxL
                write(LogID,"(A7,F10.3,3F10.2, I4,F10.4)") 'DPsig',NEQ_DP%sigma0, NEQ_DP%InvN(1,1), NEQ_DP%InvN(2,2), NEQ_DP%InvN(3,3), NEQ_DP%PRN(maxL),NEQ_DP%maxV
    
!                outlier=5.d0
                if ( (dabs(NEQ_DP%maxV)>5.d0)  ) then  !  if  (dabs(maxV)>3*sigma0) then ! .and. (dabs(maxV)>3*sigma0)   0.25
                    write(LogID,"(A10,I4,F10.3)") "DPoutlier", NEQ_DP%PRN(maxL),NEQ_DP%maxV
                    call Change_NEQ(NEQ_DP%Nbb, NEQ_DP%U, 4, NEQ_DP%A(maxL,:), NEQ_DP%L(maxL), "cut")
                    NEQ_DP%A(maxL,:)=0.d0
                    NEQ_DP%L(maxL)=0.d0
                    N=N-1
                    Ad_Flag=.true.
                end if
                if (N<6 .or. NEQ_DP%PRNS-N>3) then
                    NEQ_DP%dx=0.d0
                    exit
                end if
            end do
        end if
        write(LogID,"(A7,3F10.3)") 'Vel', NEQ_DP%dx(1:3)
    end if

    ! Step1:
    !     Solve the normal equation and get the wide line
    !     and W4 combination ambiguity (float)
    Ad_Flag=.true.
    NEQ%maxV=0.d0
    NEQ%maxL=0
    N=NEQ%PRNS
    do while(AD_flag)
        Ad_Flag=.false.
        call Invsqrt(NEQ%Nbb, NEQ%N, NEQ%InvN)
        NEQ%dx=matmul(NEQ%InvN, NEQ%U)   ! In distance(meter)
        if (any(isnan(NEQ%dx))) then
            write(*,*)  "-------ERROR-------: NEQ%dx=nan, please check"
            write(LogID, '(5X,A40)')  "-------ERROR-------: NEQ%dx=nan, please check"
            stop
        end if

        ! =================== Outliers Detect =====================
        NEQ%Vp1(1:N)=matmul(NEQ%Ap1(1:N, :), NEQ%dx(1:ParaNum)) - NEQ%Lp1(1:N)
        NEQ%maxV(1:1)=maxval(dabs(NEQ%Vp1(1:N)))
        NEQ%maxL(1:1)=maxloc(dabs(NEQ%Vp1(1:N)))
        NEQ%Vp2(1:N)=matmul(NEQ%Ap2(1:N, :), NEQ%dx(1:ParaNum)) - NEQ%Lp2(1:N)
        NEQ%maxV(2:2)=maxval(dabs(NEQ%Vp2(1:N)))
        NEQ%maxL(2:2)=maxloc(dabs(NEQ%Vp2(1:N)))
        NEQ%Vwl(1:N)=matmul(NEQ%Awl(1:N, :), NEQ%dx) - NEQ%Lwl(1:N)
        NEQ%maxV(3:3)=maxval(dabs(NEQ%Vwl(1:N)))
        NEQ%maxL(3:3)=maxloc(dabs(NEQ%Vwl(1:N)))
        NEQ%Vw4(1:N)=matmul(NEQ%Aw4(1:N, :), NEQ%dx) - NEQ%Lw4(1:N)
        NEQ%maxV(4:4)=maxval(dabs(NEQ%Vw4(1:N)))
        NEQ%maxL(4:4)=maxloc(dabs(NEQ%Vw4(1:N)))
        NEQ%Vewl(1:N)=matmul(NEQ%Aewl(1:N, :), NEQ%dx) - NEQ%Lewl(1:N)
        NEQ%maxV(5:5)=maxval(dabs(NEQ%Vewl(1:N)))
        NEQ%maxL(5:5)=maxloc(dabs(NEQ%Vewl(1:N)))

        maxV=maxval(dabs(NEQ%maxV))
        maxL=maxloc(dabs(NEQ%maxV),dim=1)
        write(unit=LogID,fmt='(5X,A5)',advance='no') 'PRN'
        do i=1,N
            write(unit=LogID,fmt='(I7)',advance='no') NEQ%PRN(i)
        end do
        write(unit=LogID,fmt='(A)') ''
        write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vp1'
        do i=1,N
            write(unit=LogID,fmt='(F7.3)',advance='no') NEQ%Vp1(i)
        end do
        write(unit=LogID,fmt='(A)') ''
        write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vp2'
        do i=1,N
            write(unit=LogID,fmt='(F7.3)',advance='no') NEQ%Vp2(i)
        end do
        write(unit=LogID,fmt='(A)') ''
        write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vwl'
        do i=1,N
            write(unit=LogID,fmt='(F7.3)',advance='no') NEQ%Vwl(i)
        end do
        write(unit=LogID,fmt='(A)') ''
        write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vw4'
        do i=1,N
            write(unit=LogID,fmt='(F7.3)',advance='no') NEQ%Vw4(i)
        end do
        write(unit=LogID,fmt='(A)') ''
        write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vewl'
        do i=1,N
            write(unit=LogID,fmt='(F7.3)',advance='no') NEQ%Vewl(i)
        end do
        write(unit=LogID,fmt='(A)') ''
        
        if ( dabs(maxV)>.4d0*a1*154.d0/(a1*154.d0+a2*120.d0)  ) then
            if ((maxL)==1) then   ! maxV in P1
                call Minus_NEQ( NEQ%Nbb(1:ParaNum,1:ParaNum), NEQ%U(1:ParaNum), NEQ%Ap1(1:N,:), NEQ%Lp1(1:N), &
                       NEQ%P(1:N, 1:N), ParaNum,  N, NEQ%maxL(1), NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in P1', 'PRN=',NEQ%PRN(NEQ%maxL(1)),'maxV=',maxV
            elseif  ((maxL)==2) then   ! maxV in P2
                call Minus_NEQ( NEQ%Nbb(1:ParaNum,1:ParaNum), NEQ%U(1:ParaNum), NEQ%Ap2(1:N,:), NEQ%Lp2(1:N), &
                       NEQ%P(1:N, 1:N), ParaNum,  N, NEQ%maxL(2), NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in P2', 'PRN=',NEQ%PRN(NEQ%maxL(2)),'maxV=',maxV
            elseif  ((maxL)==3) then
                call Minus_NEQ( NEQ%Nbb, NEQ%U, NEQ%Awl(1:N,:), NEQ%Lwl(1:N), &
                       NEQ%P(1:N, 1:N), NEQ%N,  N, NEQ%maxL(3), NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in WL', 'PRN=',NEQ%PRN(NEQ%maxL(3)),'maxV=',maxV
                if (ar_mode==3) then ! If fixed and hold mode
                    NEQ%fixed_amb(NEQ%PRN(NEQ%maxL(3)))=0.99d0  ! If outliers, do not hold the amiguity
                    NEQ%fixed_amb_num(NEQ%PRN(NEQ%maxL(3)))=0
                end if
                if (Var_smooth=='y' .or. Var_smooth=='Y') then ! If smooth pseudorange
                    STA%STA(1)%Pre(NEQ%PRN(NEQ%maxL(3)))%n=1.d0
                    STA%STA(2)%Pre(NEQ%PRN(NEQ%maxL(3)))%n=1.d0 
                end if
            elseif  ((maxL)==4) then
                call Minus_NEQ( NEQ%Nbb, NEQ%U, NEQ%Aw4(1:N,:), NEQ%Lw4(1:N), &
                       NEQ%P(1:N, 1:N), NEQ%N,  N, NEQ%maxL(4), NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in W4', 'PRN=',NEQ%PRN(NEQ%maxL(4)),'maxV=',maxV
                if (ar_mode==3) then ! If fixed and hold mode
                    NEQ%fixed_amb(NEQ%PRN(NEQ%maxL(4))+MaxPRN)=0.99d0  ! If outliers, do not hold the amiguity
                    NEQ%fixed_amb_num(NEQ%PRN(NEQ%maxL(4))+MaxPRN)=0
                end if
                if (Var_smooth=='y' .or. Var_smooth=='Y') then ! If smooth pseudorange
                    STA%STA(1)%Pre(NEQ%PRN(NEQ%maxL(4)))%n=1.d0
                    STA%STA(2)%Pre(NEQ%PRN(NEQ%maxL(4)))%n=1.d0 
                end if
            elseif  ((maxL)==5) then   ! maxV in EWL
                call Minus_NEQ( NEQ%Nbb(1:ParaNum,1:ParaNum), NEQ%U(1:ParaNum), NEQ%Aewl(1:N,:), NEQ%Lewl(1:N), &
                       NEQ%P(1:N, 1:N), ParaNum,  N, NEQ%maxL(5), NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in EWL', 'PRN=',NEQ%PRN(NEQ%maxL(5)),'maxV=',maxV
                NEQ%amb_EWL(NEQ%PRN(NEQ%maxL(5)))=0.d0
            end if
!            if (NEQ%SumN<3+N*2) exit
            Ad_Flag=.true.
            write(unit=LogID,fmt='(5X,A5,3F10.3,A8)') '!!!dx',NEQ%dx(1:3),'outlier'
        end if 
        ! =================== End of Outliers Detect =====================
        write(unit=LogID,fmt='(A10,3F7.2)') 'dx_float',NEQ%dx(1:3)
    end do
    ! If LC combination
    if (  ( (a1==0.d0) .and. (a2==0.d0) .and. (mod(b1,1.d0)/=0.d0) ) .or. ( (b1==0.d0) .and. (b2==0.d0) .and. (mod(a1,1.d0)/=0.d0) )  ) then
        Coor(1:3)=NEQ%dx(1:3)
        return
    end if
    ! If float solution
    if (ar_mode==0) then
        Coor(1:3)=NEQ%dx(1:3)
        Flag_Sln=3
        return
    end if
    
    ! Get the wide line and W4 combination ambigulties
    NEQ%amb_WL=NEQ%dx(ParaNum+1:ParaNum+MaxPRN)  ! *(a1*f1+a2*f2)/c   ! In cycle
    NEQ%amb_W4=NEQ%dx(ParaNum+1+MaxPRN:ParaNum+MaxPRN*2) ! *(b1*f1+b2*f2)/c
    

    ! Step2:
    !     Fix the WL and W4 amiguity using LAMBDA
    amb2=NEQ%dx(ParaNum+1:ParaNum+MaxPRN*2)
    do i=1, 2*maxPRN
        iPOS(i:i)=i
        if (ar_mode/=2 .and. i<=maxPRN) then
            if (NEQ%Ele(i)<FixEle) amb2(i)=0.d0    ! Fix ambiguity only when satellite elevation>FixEle
        elseif (ar_mode/=2 .and. i>maxPRN) then
            if (NEQ%Ele(i-maxPRN)<FixEle) amb2(i)=0.d0
        end if
    end do
    call LAMBDA_Prepare(NEQ%InvN(ParaNum+1:ParaNum+MaxPRN*2,ParaNum+1:ParaNum+MaxPRN*2), amb2, &
                                        MaxPRN*2, Q, amb, npar, iPOS)
    Q2=Q; amb2=amb; npar2=npar; iPOS2=iPOS; iPOS3=iPOS; P=NEQ%P
    flag_partial=0; ratio2=0.d0; k=0; m=1; l=1; minLL=0; minLL2=0
    write(LogID,'(A10)',advance='no') 'amb_float'
    if ( (a1*f1+a2*f2/=0.d0) .and. (b1*f1+b2*f2/=0.d0) ) then ! Dual frequency
        do i=1,npar
            write(LogID,'(I6,F10.1)',advance='no') iPOS2(i), amb(i)
            if (i==npar/2) then
                write(LogID,'(A)') ''
                write(LogID,'(A10)',advance='no') ''
            end if
        end do
    elseif (a1*f1+a2*f2/=0.d0) then ! L1 frequency
        do i=1,npar
            write(LogID,'(I6,2F8.1)',advance='no') iPOS2(i), amb(i)
        end do
    elseif (b1*f1+b2*f2/=0.d0) then ! L2 frequency
        do i=1,npar
            write(LogID,'(I6,2F8.1)',advance='no') iPOS2(i)-maxPRN, amb(i)
        end do
    end if
    write(LogID,'(A)') ''
    100 if (npar>1) then
!        write(LambdaID,'(A)') 'InvN'
!        write(LambdaID,'(14F11.3)') Q(1:npar, 1:npar)
!        write(LambdaID,'(14F11.3)') amb(1:npar)
!        call LAMBDA_zhang(lambdaID, npar, Q(1:npar, 1:npar), amb(1:npar), disall)
        call LAMBDA(lambdaID, npar, amb(1:npar),Q(1:npar, 1:npar)/10000.d0,1,amb(1:npar),disall,Ps,Qzhat(1:npar, 1:npar),Z(1:npar, 1:npar),nfixed,mu,dz(1:npar))
        if (nfixed==0) then
            ratio=0.d0
        else
            ratio=dabs(disall(2)/disall(1))
        end if
    end if
    if (ratio>minratio) then  ! ambiguity fix success
        if (flag_partial==1) then ! If partia       l ambiguity fixed
            ! recover the order of amb and iPOS
            do i=npar2,1,-1
                PRN=iPOS2(i)
                flag_fixed=0
                do j=i,1,-1
                    if (iPOS(j)==PRN) then
                        amb(i)=amb(j)
                        flag_fixed=1
                        exit
                    end if
                end do
                iPOS(i)=PRN
                if (flag_fixed==0) then
                    amb(i)=amb2(i)  ! unfixed ambiguity
                    iPOS(i)=0
                end if
            end do
            npar=npar2
        end if  ! if (flag_partial==1) then
        write(LogID,'(A10)',advance='no') 'amb_fix'
        if ( (a1*f1+a2*f2/=0.d0) .and. (b1*f1+b2*f2/=0.d0) ) then ! Dual frequency
            do i=1,npar
                if (iPOS(i)==0) then
                    if (iPOS2(i)<=maxPRN) then
                        NEQ%amb_WL(iPOS2(i))=0.d0
                    else
                        NEQ%amb_W4(iPOS2(i)-maxPRN)=0.d0
                    end if
                    NEQ%fixed_amb_num(iPOS2(i))=0
                else
                    if (abs(NEQ%fixed_amb(iPOS2(i))-amb(i))>0.001d0) then
                        NEQ%fixed_amb_num(iPOS2(i))=0  ! If fixed ambiguity  change
                    end if
                    NEQ%fixed_amb(iPOS2(i))=amb(i)   ! Save the fixed ambiguity
                    if (iPOS2(i)<=maxPRN) then
                        NEQ%amb_WL(iPOS2(i))=amb(i)
                    else
                        NEQ%amb_W4(iPOS2(i)-maxPRN)=amb(i)
                    end if
                end if
                write(LogID,'(I6,F8.1)',advance='no') iPOS2(i), amb(i)
                if (i==npar/2) then
                    write(LogID,'(A)') ''
                    write(LogID,'(A10)',advance='no') ''
                end if
            end do
        elseif (a1*f1+a2*f2/=0.d0) then ! L1 frequency
            do i=1,npar
                if (iPOS(i)==0) then
                    NEQ%amb_WL(iPOS2(i))=0.d0
                    NEQ%fixed_amb_num(iPOS2(i))=0
                else
                    if (abs(NEQ%fixed_amb(iPOS2(i))-amb(i))>0.001d0) then
                        NEQ%fixed_amb_num(iPOS2(i))=0  ! If fixed ambiguity change
                    end if
                    NEQ%fixed_amb(iPOS2(i))=amb(i)   ! Save the fixed ambiguity
                    NEQ%amb_WL(iPOS2(i))=amb(i)
                end if
                write(LogID,'(I6,F8.1)',advance='no') iPOS2(i), amb(i)
            end do
        elseif (b1*f1+b2*f2/=0.d0) then ! L2 frequency
            do i=1,npar
                if (iPOS(i)==0) then
                    NEQ%amb_W4(iPOS2(i)-maxPRN)=0.d0
                    NEQ%fixed_amb_num(iPOS2(i))=0
                else
                    if (abs(NEQ%fixed_amb(iPOS2(i))-amb(i))>0.001d0) then
                        NEQ%fixed_amb_num(iPOS2(i))=0  ! If fixed ambiguity does not change
                    end if
                    NEQ%fixed_amb(iPOS2(i))=amb(i)   ! Save the fixed ambiguity
                    NEQ%amb_W4(iPOS2(i)-maxPRN)=amb(i)
                end if
                write(LogID,'(I6,F8.1)',advance='no') iPOS2(i)-maxPRN, amb(i)
            end do
        end if
        write(LogID,'(A)') ''
        amb_success=amb_success+1
    elseif (partial_AR) then ! partial ambigulty fixing
        minP=100.d0
        maxQ=0.d0
        minL=0; minL2=0
        amb3=amb2
        iPOS=iPOS2
        do i=1,npar2  ! Start from the minimum elevation satellite
            PRN=iPOS3(i)
            if (PRN>maxPRN) PRN=PRN-maxPRN
            do j=1,NEQ%PRNS
                if (NEQ%PRN(j)==PRN) then
                    exit
                end if
            end do
            if ( (PRN/=0) .and. abs(P(j,j)-minP)<1.d-11 ) then
                minL2=i   ! dual frequency
            elseif ( (PRN/=0) .and. (P(j,j)<minP) ) then !    Q2(i, i)>maxQ
                minP=P(j,j)
                maxQ=Q(i,i)
                minL=i
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
            goto 100
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
            if (m<=k) then
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
                goto 100

            end if
        end if  !  if (minL/=0) then

    end if
    Flag_Sln=flag_partial+1
    NEQ%ratio=ratio
    if (ratio<minratio) then
        Flag_Sln=3   ! float solution
        Coor(1:3)=NEQ%dx(1:3)
    end if
    write(LogID,'(A10,F7.2, I7, I6)') 'ratio',ratio, Flag_Sln, amb_success

!    Coor=NEQ%dx(1:3)
!    write(unit=LogID,fmt='(A10,3F7.2)') 'fixcoor',NEQ%dx(1:3)
!    return

    !    Update the coordinate by fixing the WL ambiguity
    !   



    
    ! Step3:
    !      Change WL obaservation to pseudo-range and add to EPO_NEQ;
    N=NEQ%PRNS
    !¡¡3.1  First change the constant part of error equation in WL and W4
    do i=1, N
        PRN=Epo_NEQ%PRN(i)
        if (Epo_NEQ%Lwl(i)/=0.d0) then
            flag_par_PRN=.false.
            do j=1,npar
                if (iPOS(j)==PRN) then
                    flag_par_PRN=.true.
                    exit
                end if
            end do
            if (flag_par_PRN) then  ! If fixed
                Epo_NEQ%Lwl(i)= Epo_NEQ%Lwl(i) -  NEQ%amb_WL(PRN)*Epo_NEQ%Awl(i,ParaNum+PRN)
!            if (Epo_NEQ%amb_WL(PRN)/=99.d0) then  ! Just for test, not very good, because of the wrong rounding integer
!                Epo_NEQ%Lwl(i)= Epo_NEQ%Lwl(i) -  Epo_NEQ%amb_WL(PRN)*Epo_NEQ%Awl(i,ParaNum+PRN)
                Epo_NEQ%Awl(i,ParaNum+PRN)= 0.d0
!                ! Add N1-N2 constraints instead of Wide lane observation
!                EPO_NEQ%Nbb(ParaNum+PRN,ParaNum+PRN)=EPO_NEQ%Nbb(ParaNum+PRN,ParaNum+PRN)+1.d0/16.d0
!                EPO_NEQ%Nbb(ParaNum+MaxPRN+PRN,ParaNum+MaxPRN+PRN)=EPO_NEQ%Nbb(ParaNum+MaxPRN+PRN,ParaNum+MaxPRN+PRN)+1.d0/16.d0
!                EPO_NEQ%Nbb(ParaNum+MaxPRN+PRN,ParaNum+PRN)=EPO_NEQ%Nbb(ParaNum+MaxPRN+PRN,ParaNum+PRN)-1.d0/16.d0
!                EPO_NEQ%Nbb(ParaNum+PRN,ParaNum+MaxPRN+PRN)=EPO_NEQ%Nbb(ParaNum+PRN,ParaNum+MaxPRN+PRN)-1.d0/16.d0
!                EPO_NEQ%U(ParaNum+PRN)=EPO_NEQ%U(ParaNum+PRN)+NEQ%amb_WL(PRN)/16.d0  ! In cycle
!                EPO_NEQ%U(ParaNum+MaxPRN+PRN)=EPO_NEQ%U(ParaNum+MaxPRN+PRN)-NEQ%amb_WL(PRN)/16.d0
                ! EWL-WL ambiguity   ! No need, because WL already constraint the N1-N2 ambiguity. If triple frequency, then this can be can be applied 
                if (NEQ%amb_EWL(PRN)/=99.d0 .and. Epo_NEQ%Lw4(i)/=0.d0) then  ! 99 is in case that integer amb_EWL is 0
                    NEQ%amb_W4(PRN)=NEQ%amb_EWL(PRN)+NEQ%amb_WL(PRN)
                    Epo_NEQ%Lw4(i)= Epo_NEQ%Lw4(i) -  NEQ%amb_W4(PRN)* Epo_NEQ%Aw4(i,ParaNum+MaxPRN+PRN)
                    Epo_NEQ%Aw4(i,ParaNum+MaxPRN+PRN)= 0.d0
!                    ! Add N1-N3 constraints instead of Wide lane observation, but first should set 3 sets of ambiguities, i.e. L1, L2, L3
!                    EPO_NEQ%Nbb(ParaNum+PRN,ParaNum+PRN)=EPO_NEQ%Nbb(ParaNum+PRN,ParaNum+PRN)+1.d0
!                    EPO_NEQ%Nbb(ParaNum+2*MaxPRN+PRN,ParaNum+2*MaxPRN+PRN)=EPO_NEQ%Nbb(ParaNum+2*MaxPRN+PRN,ParaNum+2*MaxPRN+PRN)+1.d0
!                    EPO_NEQ%Nbb(ParaNum+2*MaxPRN+PRN,ParaNum+PRN)=EPO_NEQ%Nbb(ParaNum+2*MaxPRN+PRN,ParaNum+PRN)-1.d0
!                    EPO_NEQ%Nbb(ParaNum+PRN,ParaNum+2*MaxPRN+PRN)=EPO_NEQ%Nbb(ParaNum+PRN,ParaNum+2*MaxPRN+PRN)-1.d0
!                    EPO_NEQ%U(ParaNum+PRN)=EPO_NEQ%U(ParaNum+PRN)+1.d0*NEQ%amb_W4(PRN)  ! In cycle
!                    EPO_NEQ%U(ParaNum+2*MaxPRN+PRN)=EPO_NEQ%U(ParaNum+2*MaxPRN+PRN)-1.d0*NEQ%amb_W4(PRN)
                else
                    Epo_NEQ%Lw4(i)=0.d0
                    Epo_NEQ%Aw4(i,:)= 0.d0
                end if
            else
                Epo_NEQ%Lwl(i)=0.d0
                Epo_NEQ%Awl(i,:)= 0.d0
                Epo_NEQ%Lw4(i)=0.d0
                Epo_NEQ%Aw4(i,:)= 0.d0
            end if
        end if
    end do

    !    3.2  Then add to EPO_NEQ
    Epo_NEQ%Nbb     =  Epo_NEQ%Nbb +  matmul(  matmul( transpose(Epo_NEQ%Awl(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Awl(1:N,:)  )
    Epo_NEQ%Nbb     =  Epo_NEQ%Nbb +  matmul(  matmul( transpose(Epo_NEQ%Aw4(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Aw4(1:N,:)  )
    Epo_NEQ%U         =Epo_NEQ%U+  matmul(  matmul( transpose(Epo_NEQ%Awl(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Lwl(1:N) )
    Epo_NEQ%U         =Epo_NEQ%U+  matmul(  matmul( transpose(Epo_NEQ%Aw4(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Lw4(1:N)  )
    ! Add Wide lane (N1-N2) ambiguity constraints: No need, because actually ¡°pseudo-range¡± WL observation is the N1-N2 ambiguity constraint

    write(unit=LogID,fmt='(5X,A5)',advance='no') 'PRN'
    do i=1,N
        write(unit=LogID,fmt='(I7)',advance='no') Epo_NEQ%PRN(i)
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Lwl'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Lwl(i)*2.d0 ! sqrt(a1**2+a2**2)
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Lw4'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Lw4(i)*2.d0 ! sqrt(2.d0)
    end do
    write(unit=LogID,fmt='(A)') ''
    

    ! Step4:
    !      Solve L1 and L2 ambiguity(float) and ionosphere delay
    Ad_Flag=.true.
    N=Epo_NEQ%PRNS
    do while(AD_flag)
        Ad_Flag=.false.
        call Invsqrt(Epo_NEQ%Nbb, Epo_NEQ%N, Epo_NEQ%InvN)   ! Epo_NEQ%N=ParaNum+3*IonoNum
        Epo_NEQ%dx = matmul(Epo_NEQ%InvN, Epo_NEQ%U)
        if (any(isnan(Epo_NEQ%dx))) then
            write(*,*)  "-------ERROR-------: Epo_NEQ%dx=nan, please check"
            write(LogID, '(5X,A40)')  "-------ERROR-------: Epo_NEQ%dx=nan, please check"
            stop
        end if
        
        ! =================== Outliers Detect =====================
        Epo_NEQ%Vp1(1:N)=matmul(Epo_NEQ%Ap1(1:N, :), Epo_NEQ%dx) - Epo_NEQ%Lp1(1:N)
        Epo_NEQ%maxV(1:1)=maxval(dabs(Epo_NEQ%Vp1(1:N)))
        Epo_NEQ%maxL(1:1)=maxloc(dabs(Epo_NEQ%Vp1(1:N)))
        Epo_NEQ%Vp2(1:N)=matmul(Epo_NEQ%Ap2(1:N, :), Epo_NEQ%dx) - Epo_NEQ%Lp2(1:N)
        Epo_NEQ%maxV(2:2)=maxval(dabs(Epo_NEQ%Vp2(1:N)))
        Epo_NEQ%maxL(2:2)=maxloc(dabs(Epo_NEQ%Vp2(1:N)))
        Epo_NEQ%Vl1(1:N)=matmul(Epo_NEQ%Al1(1:N, :), Epo_NEQ%dx) - Epo_NEQ%Ll1(1:N)
        Epo_NEQ%maxV(3:3)=maxval(dabs(Epo_NEQ%Vl1(1:N)))
        Epo_NEQ%maxL(3:3)=maxloc(dabs(Epo_NEQ%Vl1(1:N)))
        Epo_NEQ%Vl2(1:N)=matmul(Epo_NEQ%Al2(1:N, :), Epo_NEQ%dx) - Epo_NEQ%Ll2(1:N)
        Epo_NEQ%maxV(4:4)=maxval(dabs(Epo_NEQ%Vl2(1:N)))
        Epo_NEQ%maxL(4:4)=maxloc(dabs(Epo_NEQ%Vl2(1:N)))
        Epo_NEQ%Vwl(1:N)=matmul(Epo_NEQ%Awl(1:N, :), Epo_NEQ%dx) - Epo_NEQ%Lwl(1:N)
        Epo_NEQ%maxV(5:5)=maxval(dabs(Epo_NEQ%Vwl(1:N)))
        Epo_NEQ%maxL(5:5)=maxloc(dabs(Epo_NEQ%Vwl(1:N)))
        Epo_NEQ%Vw4(1:N)=matmul(Epo_NEQ%Aw4(1:N, :), Epo_NEQ%dx) - Epo_NEQ%Lw4(1:N)
        Epo_NEQ%maxV(6:6)=maxval(dabs(Epo_NEQ%Vw4(1:N)))
        Epo_NEQ%maxL(6:6)=maxloc(dabs(Epo_NEQ%Vw4(1:N)))

        maxV=maxval(dabs(Epo_NEQ%maxV))
        maxL=maxloc(dabs(Epo_NEQ%maxV),dim=1)
        
        if ( (dabs(maxV)>.08d0*a1*154.d0/(a1*154.d0+a2*120.d0))  ) then
            if ((maxL)==1) then   ! maxV in P1
                call Minus_NEQ( Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%Ap1(1:N,:), Epo_NEQ%Lp1(1:N), &
                       Epo_NEQ%P(1:N, 1:N), Epo_NEQ%N,  N, Epo_NEQ%maxL(1), Epo_NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in P1', 'PRN=',Epo_NEQ%PRN(Epo_NEQ%maxL(1)),'maxV=',maxV
            elseif  ((maxL)==2) then   ! maxV in P2
                call Minus_NEQ( Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%Ap2(1:N,:), Epo_NEQ%Lp2(1:N), &
                       Epo_NEQ%P(1:N, 1:N), Epo_NEQ%N,  N, Epo_NEQ%maxL(2), Epo_NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in P2', 'PRN=',Epo_NEQ%PRN(Epo_NEQ%maxL(2)),'maxV=',maxV
            elseif  ((maxL)==3) then
                call Minus_NEQ( Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%Al1(1:N,:), Epo_NEQ%Ll1(1:N), &
                       Epo_NEQ%P(1:N, 1:N), Epo_NEQ%N,  N, Epo_NEQ%maxL(3), Epo_NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in L1', 'PRN=',Epo_NEQ%PRN(Epo_NEQ%maxL(3)),'maxV=',maxV
                if (ar_mode==3) then ! If fixed and hold mode
                    Epo_NEQ%fixed_amb(Epo_NEQ%PRN(Epo_NEQ%maxL(3)))=0.99d0  ! If outliers, do not hold the amiguity
                    Epo_NEQ%fixed_amb_num(Epo_NEQ%PRN(Epo_NEQ%maxL(3)))=0
                end if
                if (Var_smooth=='y' .or. Var_smooth=='Y') then ! If smooth pseudorange
                    STA%STA(1)%Pre(Epo_NEQ%PRN(Epo_NEQ%maxL(3)))%n=1.d0
                    STA%STA(2)%Pre(Epo_NEQ%PRN(Epo_NEQ%maxL(3)))%n=1.d0 
                end if
            elseif  ((maxL)==4) then
                call Minus_NEQ( Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%Al2(1:N,:), Epo_NEQ%Ll2(1:N), &
                       Epo_NEQ%P(1:N, 1:N), Epo_NEQ%N,  N, Epo_NEQ%maxL(4), Epo_NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in L2', 'PRN=',Epo_NEQ%PRN(Epo_NEQ%maxL(4)),'maxV=',maxV
                if (ar_mode==3) then ! If fixed and hold mode
                    Epo_NEQ%fixed_amb(Epo_NEQ%PRN(Epo_NEQ%maxL(4))+MaxPRN)=0.99d0  ! If outliers, do not hold the amiguity
                    Epo_NEQ%fixed_amb_num(Epo_NEQ%PRN(Epo_NEQ%maxL(4))+MaxPRN)=0
                end if
                if (Var_smooth=='y' .or. Var_smooth=='Y') then ! If smooth pseudorange
                    STA%STA(1)%Pre(Epo_NEQ%PRN(Epo_NEQ%maxL(4)))%n=1.d0
                    STA%STA(2)%Pre(Epo_NEQ%PRN(Epo_NEQ%maxL(4)))%n=1.d0 
                end if
           elseif ((maxL)==5) then   ! maxV in WL
                call Minus_NEQ( Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%Awl(1:N,:), Epo_NEQ%Lwl(1:N), &
                       Epo_NEQ%P(1:N, 1:N), Epo_NEQ%N,  N, Epo_NEQ%maxL(5), Epo_NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in WL', 'PRN=',Epo_NEQ%PRN(Epo_NEQ%maxL(5)),'maxV=',maxV
            elseif ((maxL)==6) then   ! maxV in W4
                call Minus_NEQ( Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%Aw4(1:N,:), Epo_NEQ%Lw4(1:N), &
                       Epo_NEQ%P(1:N, 1:N), Epo_NEQ%N,  N, Epo_NEQ%maxL(6), Epo_NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in W4', 'PRN=',Epo_NEQ%PRN(Epo_NEQ%maxL(6)),'maxV=',maxV
            end if

            Ad_Flag=.true.
        end if
        ! =================== End of Outliers Detect =====================
        write(unit=LogID,fmt='(A10,4F7.2)') 'dx_L1_float',Epo_NEQ%dx(1:4)
    end do
    
    ! Write residuals
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'PRN'
    do i=1,N
        write(unit=LogID,fmt='(I7)',advance='no') Epo_NEQ%PRN(i)
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vp1'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Vp1(i)
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vp2'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Vp2(i)
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vl1'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Vl1(i)
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vl2'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Vl2(i)
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vwl'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Vwl(i)
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vw4'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Vw4(i)
    end do
    write(unit=LogID,fmt='(A)') ''
    
    ! Write ionosphere parameter
    write(LogID,'(A10)',advance='no') 'iono'
    do i=ParaNum+2*IonoNum+1, ParaNum+3*IonoNum
        if (Epo_NEQ%dx(i)/=0.d0) then
            write(LogID,'(I6,F8.2)',advance='no') i-ParaNum-2*IonoNum, Epo_NEQ%dx(i)
        end if
    end do
    write(LogID,'(A)') ''

    ! Get  L1 and L2 ambigulties
    Epo_NEQ%amb_L1=Epo_NEQ%dx(ParaNum+1:ParaNum+MaxPRN)   ! In cycle
    Epo_NEQ%amb_L2=Epo_NEQ%dx(ParaNum+1+MaxPRN:ParaNum+MaxPRN*2)
    
    
    ! Step5:
    !     Fix  L1 and L2 amiguity using LAMBDA
    ! call LAMBDA(NEQ%InvN,dx, fixed_dx)
    amb2=Epo_NEQ%dx(ParaNum+1:ParaNum+MaxPRN*2)
    do i=1, 2*maxPRN
        iPOS(i:i)=i
        if (ar_mode/=2 .and. i<=maxPRN) then
            if (Epo_NEQ%Ele(i)<FixEle) amb2(i)=0.d0    ! Fix ambiguity only when satellite elevation>FixEle
        elseif (ar_mode/=2 .and. i>maxPRN) then
            if (Epo_NEQ%Ele(i-maxPRN)<FixEle) amb2(i)=0.d0
        end if
    end do
    call LAMBDA_Prepare(Epo_NEQ%InvN(ParaNum+1:ParaNum+MaxPRN*2,ParaNum+1:ParaNum+MaxPRN*2), amb2, &
                                       MaxPRN*2, Q, amb, npar, iPOS)
    Q2=Q; amb2=amb; npar2=npar; iPOS2=iPOS; iPOS3=iPOS; P=Epo_NEQ%P
    flag_partial=0; ratio2=0.d0; k=0; m=1; l=1; minLL=0; minLL2=0
    write(LogID,'(A10)',advance='no') 'L1_amb_float'
    do i=1,npar  ! Dual frequency
        write(LogID,'(I6,F8.1)',advance='no') iPOS2(i), amb(i)
        if (i==npar/2) then
            write(LogID,'(A)') ''
            write(LogID,'(A10)',advance='no') ''
        end if
    end do
    write(LogID,'(A)') ''

    200 if (npar>1) then
        call LAMBDA(lambdaID, npar, amb(1:npar),Q(1:npar, 1:npar)/10000.d0,1,amb(1:npar),disall,Ps,Qzhat(1:npar, 1:npar),Z(1:npar, 1:npar),nfixed,mu,dz(1:npar))
        if (nfixed==0) then
            ratio=0.d0
        else
            ratio=dabs(disall(2)/disall(1))
        end if
    end if
    if (ratio>minratio) then  ! ambiguity fix success
        if (flag_partial==1) then ! If partia       l ambiguity fixed
            ! recover the order of amb and iPOS
            do i=npar2,1,-1
                PRN=iPOS2(i)
                flag_fixed=0
                do j=i,1,-1
                    if (iPOS(j)==PRN) then
                        amb(i)=amb(j)
                        flag_fixed=1
                        exit
                    end if
                end do
                iPOS(i)=PRN
                if (flag_fixed==0) then
                    amb(i)=amb2(i)  ! unfixed ambiguity
                    iPOS(i)=0
                end if
            end do
            npar=npar2
            
        end if  ! if (flag_partial==1) then
        write(LogID,'(A10)',advance='no') 'amb_fix'

        do i=1,npar  ! Dual frequency
            if (iPOS(i)==0) then
!                if (iPOS2(i)<=maxPRN) then
!                    Epo_NEQ%amb_L1(iPOS2(i))=0.d0
!                else
!                    Epo_NEQ%amb_L2(iPOS2(i)-maxPRN)=0.d0
!                end if
                Epo_NEQ%fixed_amb_num(iPOS2(i))=0
            else
                if (abs(Epo_NEQ%fixed_amb(iPOS2(i))-amb(i))>0.001d0) then
                    Epo_NEQ%fixed_amb_num(iPOS2(i))=0  ! If fixed ambiguity change
                end if
                Epo_NEQ%fixed_amb(iPOS2(i))=amb(i)   ! Save the fixed ambiguity
                if (iPOS2(i)<=maxPRN) then
                    Epo_NEQ%amb_L1(iPOS2(i))=amb(i)
                else
                    Epo_NEQ%amb_L2(iPOS2(i)-maxPRN)=amb(i)
                end if
            end if
            write(LogID,'(I6,F8.1)',advance='no') iPOS2(i), amb(i)
            if (i==npar/2) then
                write(LogID,'(A)') ''
                write(LogID,'(A10)',advance='no') ''
            end if
        end do
        write(LogID,'(A)') ''
        amb_success2=amb_success2+1
    elseif (partial_AR) then ! partial ambigulty fixing
        minP=100.d0
        maxQ=0.d0
        minL=0; minL2=0
        amb3=amb2
        iPOS=iPOS2
        do i=1,npar2  ! Start from the minimum elevation satellite
            PRN=iPOS3(i)
            if (PRN>maxPRN) PRN=PRN-maxPRN
            do j=1,EPO_NEQ%PRNS
                if (EPO_NEQ%PRN(j)==PRN) then
                    exit
                end if
            end do
            if ( (PRN/=0) .and. abs(P(j,j)-minP)<1.d-11 ) then
                minL2=i   ! dual frequency
            elseif ( (PRN/=0) .and. (P(j,j)<minP) ) then !  Q2(i, i)>maxQ
                minP=P(j,j)
                maxQ=Q(i,i)
                minL=i
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
            goto 200
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
            if (m<=k) then
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
                goto 200

            end if
        end if  !  if (minL/=0) then

    end if
    Flag_Sln=flag_partial+1
    Epo_NEQ%ratio=ratio
    if (ratio<minratio) then
!        if (NEQ%ratio>minratio) then
!            Flag_Sln=3
!            Coor=NEQ%dx(1:3)   ! Use the wide lane ambiguity
!            write(unit=LogID,fmt='(A10,3F7.2)') 'fixcoor WL',NEQ%dx(1:3)
!        else
            Flag_Sln=3 
            Coor=Epo_NEQ%dx(1:4)
!        end if
    end if
    write(LogID,'(A10,F7.2, I7,I6)') 'ratio2',ratio, Flag_Sln, amb_success2
    if (Flag_Sln==3) return  ! Return float solution

    ! Step6:
!        Update the coordinate by fixing the ambiguity of L1 and L2.
!    write(*,*) 'Update the coordinate by fixing the ambiguity of L1 and L2'
!    dx(1:MaxPRN)=Epo_NEQ%dx(ParaNum+1:ParaNum+MaxPRN)-Epo_NEQ%amb_L1 ! residual of float ambiguity to integer ambiguity
!    dx(MaxPRN+1:2*MaxPRN)=Epo_NEQ%dx(ParaNum+MaxPRN+1:ParaNum+2*MaxPRN)-Epo_NEQ%amb_L2
!    temp_QR(1:3, 1:2*MaxPRN)=EPO_NEQ%InvN(1:3,ParaNum+1:ParaNum+2*MaxPRN)
!    temp_QR(ParaNum+1:ParaNum+MaxPRN, 1:2*MaxPRN)=EPO_NEQ%InvN(ParaNum+2*MaxPRN+1:ParaNum+3*MaxPRN,ParaNum+1:ParaNum+2*MaxPRN)
!
!    temp_QQ=EPO_NEQ%Nbb(ParaNum+1:ParaNum+MaxPRN*2,ParaNum+1:ParaNum+MaxPRN*2)
!    temp_dx(1:3)=Epo_NEQ%dx(1:3)
!    temp_dx(4:3+MaxPRN)=Epo_NEQ%dx(ParaNum+2*MaxPRN+1:ParaNum+3*MaxPRN)
!    temp_dx=temp_dx-MATMUL(MATMUL(temp_QR, temp_QQ), dx)
!    do i=1, npar
!        PRN=iPOS2(i)
!        do j=1, ParaNum
!            U(j)=U(j) - amb(i)*EPO_NEQ%Al1(,j
!        end do
!
!    end do

    dx(1:MaxPRN)=Epo_NEQ%amb_L1  ! Not best, as amb_L1 includes some unfixed float solution, should abandon them, but maybe absorbed in ionosphere parameter?
    dx(MaxPRN+1:2*MaxPRN)=Epo_NEQ%amb_L2
    temp_Nbb(1:ParaNum,1:ParaNum)=Epo_NEQ%Nbb(1:ParaNum,1:ParaNum)
    temp_Nbb(ParaNum+1:ParaNum+MaxPRN,1:ParaNum)=Epo_NEQ%Nbb(ParaNum+MaxPRN*2+1:ParaNum+MaxPRN*3,1:ParaNum)
    temp_Nbb(ParaNum+1:ParaNum+MaxPRN,ParaNum+1:ParaNum+MaxPRN)=Epo_NEQ%Nbb(ParaNum+MaxPRN*2+1:ParaNum+MaxPRN*3,ParaNum+MaxPRN*2+1:ParaNum+MaxPRN*3)
    call Invsqrt(temp_Nbb, ParaNum+MaxPRN, temp_InvN)
    temp_U(1:ParaNum)=Epo_NEQ%U(1:ParaNum)
    temp_U(ParaNum+1:ParaNum+MaxPRN)=Epo_NEQ%U(ParaNum+MaxPRN*2+1:ParaNum+MaxPRN*3)
    temp_Nbb2(1:2*MaxPRN,1:ParaNum)=EPO_NEQ%Nbb(ParaNum+1:ParaNum+2*MaxPRN,1:ParaNum)
    temp_Nbb2(1:2*MaxPRN,ParaNum+1:ParaNum+MaxPRN)=transpose(EPO_NEQ%Nbb(ParaNum+2*MaxPRN+1:ParaNum+3*MaxPRN,ParaNum+1:ParaNum+2*MaxPRN))

    temp_dx=MATMUL(temp_InvN, temp_U-MATMUL(transpose(temp_Nbb2),dx))
!    Epo_NEQ%dx(1:ParaNum)=temp_dx(1:ParaNum)
!    Epo_NEQ%dx(ParaNum+2*MaxPRN+1:ParaNum+3*MaxPRN)=temp_dx(ParaNum+1:ParaNum+MaxPRN)

    
    Coor=temp_dx(1:4)
    write(unit=LogID,fmt='(A10,3F7.2)') 'fixcoor',temp_dx(1:3)
    return
    
end subroutine