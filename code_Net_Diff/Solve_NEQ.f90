! ============= Sovle_NEQ===================== 
! PURPOSE:
!            Solve the normal equation and get the coordinate.
!        Steps:
!            1) Solve the normal equation and get the wide line
!     and W4 combination ambiguity (float);
!            2) Fix the WL and W4 amiguity using LAMBDA method;
!            3.a) Use the fixed ambigity to get fixed coordinate, return;
!            3.b) Use the fixed ambigity to get the ambiguity 
!      of L1 and L2;
!            4) Estimate the coordinate by fixing the ambiguity of 
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

subroutine Solve_NEQ(NEQ, Epo_NEQ, Coor, Flag_Sln)
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
    real(8) :: Coor(3)
    real(8) :: maxV
    integer :: maxL
    real(8) :: NEQ_dx(ParaNum+2*MaxPRN), NEQ_InvN(ParaNum+2*MaxPRN,ParaNum+2*MaxPRN)
    real(8) :: AA(4*NEQ%PRNS,ParaNum+2*MaxPRN), LL(4*NEQ%PRNS), RR(4*NEQ%PRNS,4*NEQ%PRNS)
    ! for lambda
    integer :: npar, npar2, nfixed
    real(8) :: Ps, Qzhat(2*MaxPRN, 2*maxPRN), Z(2*MaxPRN, 2*maxPRN), mu
    real(8) :: Q(2*MaxPRN, 2*maxPRN), Q2(2*MaxPRN, 2*maxPRN), P(MaxPRN, maxPRN)
    real(8) :: disall(2), ratio, ratio2, amb(2*MaxPRN), amb2(2*MaxPRN), amb3(2*MaxPRN), dz(2*MaxPRN)
    integer :: iPOS(2*MaxPRN), iPOS2(2*MaxPRN), iPOS3(2*MaxPRN), minL, minL2, minLL(MaxPRN),minLL2(MaxPRN)
    real(8) :: minP, maxQ
    integer(1) :: flag_partial, flag_fixed, freq, sys, par_PRNS
    real(8) :: temp_Nbb(ParaNum,ParaNum), temp_InvN(2*MaxPRN, 2*MaxPRN),temp_InvN2(ParaNum, 2*MaxPRN), dx(2*MaxPRN), temp_dx(ParaNum) !

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
                    write(LogID, '(5X,A50)')  "-------ERROR-------: dx=nan, please check"
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
    if (ADmethod=='KF') then
        NEQ_InvN=NEQ%InvN
        NEQ_dx=NEQ%dx ! temp save Pk and X in in case of outliers
    end if
    do while(AD_flag)
        Ad_Flag=.false.
        if (ADmethod=='LS') then
            call Invsqrt(NEQ%Nbb, NEQ%N, NEQ%InvN)
            NEQ%dx=matmul(NEQ%InvN, NEQ%U)   ! In distance(meter)
        elseif (ADmethod=='KF') then
            ! Add Kk information for each type of observation and get new dx and InvN
            K=0; AA=0.d0; LL=0.d0; RR=0.d0
            if (any(NEQ%Lp1(1:N)/=0.d0)) then
                AA(K+1:K+N,1:ParaNum)=NEQ%Ap1(1:N,1:ParaNum)
                LL(K+1:K+N)=NEQ%Lp1(1:N)
                RR(K+1:K+N,K+1:K+N)=NEQ%R(1:N,1:N)
                K=K+N
            end if
            if (any(NEQ%Lp2(1:N)/=0.d0)) then
                AA(K+1:K+N,1:ParaNum)=NEQ%Ap2(1:N,1:ParaNum)
                LL(K+1:K+N)=NEQ%Lp2(1:N)
                RR(K+1:K+N,K+1:K+N)=NEQ%R(1:N,1:N)
                K=K+N
            end if
            if (any(NEQ%Lwl(1:N)/=0.d0)) then
                AA(K+1:K+N,:)=NEQ%Awl(1:N,:)
                LL(K+1:K+N)=NEQ%Lwl(1:N)
                RR(K+1:K+N,K+1:K+N)=NEQ%R(1:N,1:N)
                K=K+N
            end if
            if (any(NEQ%Lw4(1:N)/=0.d0)) then
                AA(K+1:K+N,:)=NEQ%Aw4(1:N,:)
                LL(K+1:K+N)=NEQ%Lw4(1:N)
                RR(K+1:K+N,K+1:K+N)=NEQ%R(1:N,1:N)
                K=K+N
            end if
            call KF_Gain(NEQ%InvN, NEQ%dx, NEQ%N, K, AA(1:K,1:ParaNum+2*MaxPRN), LL(1:K), RR(1:K,1:K))
        end if
        if (any(isnan(NEQ%dx))) then
            write(*,*)  "-------ERROR-------: NEQ%dx=nan, please check"
            write(LogID, '(5X,A50)')  "-------ERROR-------: NEQ%dx=nan, please check"
            stop
        end if

        ! =================== Outliers Detect =====================
        NEQ%Vp1(1:N)=(matmul(NEQ%Ap1(1:N, :), NEQ%dx(1:ParaNum)) - NEQ%Lp1(1:N))
        NEQ%Vp2(1:N)=(matmul(NEQ%Ap2(1:N, :), NEQ%dx(1:ParaNum)) - NEQ%Lp2(1:N))
        NEQ%Vwl(1:N)=(matmul(NEQ%Awl(1:N, :), NEQ%dx) - NEQ%Lwl(1:N))
        NEQ%Vw4(1:N)=(matmul(NEQ%Aw4(1:N, :), NEQ%dx) - NEQ%Lw4(1:N))
        do i=1,N
             NEQ%Vp1(i)= NEQ%Vp1(i)*sigLC*sqrt(2.d0*NEQ%P(i,i))   ! unify to  carrier phase magnitude and the same weight, max P(i,i) is 0.5
             NEQ%Vp2(i)= NEQ%Vp2(i)*sigLC*sqrt(2.d0*NEQ%P(i,i))
             NEQ%Vwl(i)= NEQ%Vwl(i)*sigLC*sqrt(2.d0*NEQ%P(i,i))
             NEQ%Vw4(i)= NEQ%Vw4(i)*sigLC*sqrt(2.d0*NEQ%P(i,i))
        end do
        NEQ%maxV(1:1)=maxval(dabs(NEQ%Vp1(1:N)))
        NEQ%maxL(1:1)=maxloc(dabs(NEQ%Vp1(1:N)))
        NEQ%maxV(2:2)=maxval(dabs(NEQ%Vp2(1:N)))
        NEQ%maxL(2:2)=maxloc(dabs(NEQ%Vp2(1:N)))
        NEQ%maxV(3:3)=maxval(dabs(NEQ%Vwl(1:N)))
        NEQ%maxL(3:3)=maxloc(dabs(NEQ%Vwl(1:N)))
        NEQ%maxV(4:4)=maxval(dabs(NEQ%Vw4(1:N)))
        NEQ%maxL(4:4)=maxloc(dabs(NEQ%Vw4(1:N)))

        maxV=maxval(dabs(NEQ%maxV))
        maxL=maxloc(dabs(NEQ%maxV),dim=1)
        
        if ( dabs(maxV)>.21d0*c/(a1*154.d0+a2*120.d0)/10.23d6  ) then
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
!                if (Var_smooth=='y' .or. Var_smooth=='Y') then ! If smooth pseudorange
!                    STA%STA(1)%Pre(NEQ%PRN(NEQ%maxL(3)))%n=1.d0
!                    STA%STA(2)%Pre(NEQ%PRN(NEQ%maxL(3)))%n=1.d0 
!                end if
            elseif  ((maxL)==4) then
                call Minus_NEQ( NEQ%Nbb, NEQ%U, NEQ%Aw4(1:N,:), NEQ%Lw4(1:N), &
                       NEQ%P(1:N, 1:N), NEQ%N,  N, NEQ%maxL(4), NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in W4', 'PRN=',NEQ%PRN(NEQ%maxL(4)),'maxV=',maxV
                if (ar_mode==3) then ! If fixed and hold mode
                    NEQ%fixed_amb(NEQ%PRN(NEQ%maxL(4))+MaxPRN)=0.99d0  ! If outliers, do not hold the amiguity
                    NEQ%fixed_amb_num(NEQ%PRN(NEQ%maxL(4))+MaxPRN)=0
                end if
!                if (Var_smooth=='y' .or. Var_smooth=='Y') then ! If smooth pseudorange
!                    STA%STA(1)%Pre(NEQ%PRN(NEQ%maxL(4)))%n=1.d0
!                    STA%STA(2)%Pre(NEQ%PRN(NEQ%maxL(4)))%n=1.d0 
!                end if
            end if

            if (ADmethod=='KF') then
                NEQ%InvN=NEQ_InvN
                NEQ%dx=NEQ_dx ! temp save Pk and X in in case of outliers
            end if
            Ad_Flag=.true.
            write(unit=LogID,fmt='(5X,A5,3F10.3,A8)') '!!!dx',NEQ%dx(1:3),'outlier'
        end if
        
        ! =================== End of Outliers Detect =====================
        write(unit=LogID,fmt='(A10,3F7.2)') 'dx_float',NEQ%dx(1:3)
    end do
    
    ! Write residuals
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'PRN'
    do i=1,N
        write(unit=LogID,fmt='(I7)',advance='no') NEQ%PRN(i)
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vp1'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') NEQ%Vp1(i)*sigPC/sigLC/sqrt(2.d0*NEQ%P(i,i))   ! back to raw residuals
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vp2'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') NEQ%Vp2(i)*sigPC/sigLC/sqrt(2.d0*NEQ%P(i,i))
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vl1'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') NEQ%Vwl(i)/sqrt(2.d0*NEQ%P(i,i))
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vl2'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') NEQ%Vw4(i)/sqrt(2.d0*NEQ%P(i,i))
    end do
    write(unit=LogID,fmt='(A)') ''

    ! If LC combination
    if (  ( (a1==0.d0) .and. (a2==0.d0) .and. (mod(b1,1.d0)/=0.d0) ) .or. ( (b1==0.d0) .and. (b2==0.d0) .and. (mod(a1,1.d0)/=0.d0) )  ) then
        Coor=NEQ%dx(1:3)
        return
    end if
    ! If float solution
    if (ar_mode==0) then
        Coor=NEQ%dx(1:3)
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
            if (NEQ%Ele(i)<FixEle) then
                amb2(i)=0.d0    ! Fix ambiguity only when satellite elevation>FixEle
            end if
        elseif (ar_mode/=2 .and. i>maxPRN) then
            if (NEQ%Ele(i-maxPRN)<FixEle) then
                amb2(i)=0.d0
            end if
        end if
    end do
    call LAMBDA_Prepare(NEQ%InvN(ParaNum+1:ParaNum+MaxPRN*2,ParaNum+1:ParaNum+MaxPRN*2), amb2, &
                                       MaxPRN*2, Q, amb, npar, iPOS)
    Q2=Q; amb2=amb; npar2=npar; iPOS2=iPOS; amb3=amb; iPOS3=iPOS; P=NEQ%P
    flag_partial=0; ratio2=0.d0; k=0; m=1; l=1; minLL=0; minLL2=0
    write(LogID,'(A10)',advance='no') 'amb_float'
    if ( (a1*f1+a2*f2/=0.d0) .and. (b1*f1+b2*f2/=0.d0) ) then ! Dual frequency
        do i=1,npar
            write(LogID,'(I6,F8.1)',advance='no') iPOS2(i), amb(i)
            if (i==npar/2) then
                write(LogID,'(A)') ''
                write(LogID,'(A10)',advance='no') ''
            end if
        end do
    elseif (a1*f1+a2*f2/=0.d0) then ! L1 frequency
        do i=1,npar
            write(LogID,'(I6,F8.1)',advance='no') iPOS2(i), amb(i)
        end do
    elseif (b1*f1+b2*f2/=0.d0) then ! L2 frequency
        do i=1,npar
            write(LogID,'(I6,F8.1)',advance='no') iPOS2(i)-maxPRN, amb(i)
        end do
    end if
    write(LogID,'(A)') ''
    100 if (npar>1) then
!        write(LambdaID,'(A)') 'InvN'
!        write(LambdaID,'(38F11.3)') Q(1:npar, 1:npar)
!        write(LambdaID,'(38F11.3)') amb(1:npar)
!        call LAMBDA_zhang(lambdaID, npar, Q(1:npar, 1:npar), amb(1:npar), disall)
        call LAMBDA(lambdaID, npar, amb(1:npar),Q(1:npar, 1:npar)/9.d0,1,amb(1:npar),disall,Ps,Qzhat(1:npar, 1:npar),Z(1:npar, 1:npar),nfixed,mu,dz(1:npar))
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
                        NEQ%fixed_amb_num(iPOS2(i))=0  ! If fixed ambiguity does not change
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
                        NEQ%fixed_amb_num(iPOS2(i))=0  ! If fixed ambiguity does not change
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
        minP=100.d0; maxQ=0.d0
        minL=0; minL2=0
        amb3=amb2
        iPOS=iPOS2
        do i=1,npar2  ! Start from the maximum variance satellite
            PRN=iPOS3(i)
            if (PRN==0) cycle
            if (PRN>maxPRN) PRN=PRN-maxPRN
            do j=1,NEQ%PRNS
                if (NEQ%PRN(j)==PRN) then
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
                if (npar>=3) goto 100

            end if
        end if  !  if (minL/=0) then

    end if
    write(LogID,'(A10,2F7.2, I6)') 'ratio',ratio, ratio2, amb_success
    if (ratio<minratio) then
        Flag_Sln=3   ! Return float solution
        Coor=NEQ%dx(1:3)
        return
    end if
    Flag_Sln=flag_partial+1

    ! Setp3(a):
!        Update the coordinate by fixing the ambiguity of L1 and L2.
!    ! If there is no L1&L2 mixed combination, Still some problem, don't use this.
    if ((a1*a2==0.d0 .and. b1*b2==0.d0) .or. (a1*f1+a2*f2==0.d0) .or. (b1*f1+b2*f2==0.d0)) then
        ! Need some convengence time of the unfixed ambiguity, as amb_WL includes some unfixed float solution, should abandon them in U(1:ParaNum)
!        dx=NEQ%dx(ParaNum+1:ParaNum+2*maxPRN)
!        do i=1,npar
!            PRN=iPOS2(i)
!            dx(PRN)=amb(i)
!        end do
!        call InvSqrt(NEQ%Nbb(1:ParaNum,1:ParaNum), ParaNum, temp_Nbb)
!        temp_dx=MATMUL(temp_Nbb, NEQ%U(1:ParaNum)-MATMUL(transpose(NEQ%Nbb(ParaNum+1:ParaNum+2*MaxPRN,1:ParaNum)),dx))
!
!        do i=1,npar  ! wrong, don't know why
!            PRN=iPOS2(i)
!            temp_InvN2(:,i)=NEQ%InvN(1:ParaNum, ParaNum+PRN)
!        end do
!        temp_InvN2(:,1:npar)=MATMUL(temp_InvN2(:,1:npar), Z(1:npar, 1:npar))  ! Qbz
!        call InvSqrt(Qzhat(1:npar, 1:npar), npar, temp_Nbb)   ! Qzz-1  ! See lambda Eq-2.7
!        temp_dx=NEQ%dx(1:ParaNum)-MATMUL(MATMUL(temp_InvN2(:,1:npar), temp_Nbb), dz(1:npar))

        dx(1:npar)=0.d0; temp_InvN2(:,1:npar)=0.d0; j=0
        do i=1,npar
            PRN=iPOS(i)
            if (PRN==0) cycle
            j=j+1
            dx(j)=NEQ%dx(ParaNum+PRN)-amb(i)
            temp_InvN2(:,j)=NEQ%InvN(1:ParaNum, ParaNum+PRN)
        end do
        call InvSqrt(Q(1:npar, 1:npar), npar, temp_InvN(1:npar,1:npar)) !  Qaa-1 ! See lambda Eq-2.27, only the fixed InvN
        temp_dx=NEQ%dx(1:ParaNum)-MATMUL(MATMUL(temp_InvN2(:,1:npar), temp_InvN(1:npar,1:npar)), dx(1:npar))

        Coor=temp_dx(1:3)
        write(unit=LogID,fmt='(A10,3F7.2)') 'fixcoor', Coor
        return
    end if
    
    ! Step3(b):
    !      Use the float or fixed ambigity to get the 
    !      ambiguity  of L1 and L2
    if ( (a1*f1+a2*f2/=0.d0) .and. (b1*f1+b2*f2/=0.d0) ) then ! Dual frequency
        NEQ%amb_L1 = (b2*NEQ%amb_WL - a2* NEQ%amb_W4)/(a1*b2-a2*b1)   ! In cycle
        NEQ%amb_L2 = (-b1*NEQ%amb_WL + a1*NEQ%amb_W4)/(a1*b2-a2*b1)
    elseif (a1*f1+a2*f2/=0.d0) then ! L1 frequency
        NEQ%amb_L1=NEQ%amb_WL
    elseif (b1*f1+b2*f2/=0.d0) then ! L2 frequency
        NEQ%amb_L2=NEQ%amb_W4
    end if

    ! Step4: 
    !      Estimate the coordinate by fixing the ambiguity of 
    !      L1 and L2.
    N=Epo_NEQ%PRNS
    par_PRNS=0
    !¡¡First change the constant part of error equation in L1 and L2
    do i=1, N
        PRN=Epo_NEQ%PRN(i)
        sys=Epo_NEQ%Sys(i)

        if (sys==1) then   ! GPS/QZSS
            f1=10.23d6*154d0
            f2=10.23d6*120d0
        elseif (sys==2) then   ! GLONASS
            freq=Fre_Chann(PRN-GNum)
            f1=(1602.0d0+freq*0.5625d0)*1.0D6   ! f1=(1602.0d0+K*0.5625d0)*1.0d6
            f2=(1246.0d0+freq*0.4375d0)*1.0D6
        elseif  (sys==3) then   ! COMPASS
            if (freq_comb=='L1L2') then
                f1=10.23d6*152.6d0
                f2=10.23d6*118.0d0
!                f3=10.23d6*124.0d0
            elseif (freq_comb=='L1L3') then
                f1=10.23d6*152.6d0
                f2=10.23d6*124.0d0
            elseif (freq_comb=='L2L3') then
                f1=10.23d6*118.0d0
                f2=10.23d6*124.0d0
            end if
        elseif  (sys==4) then ! GALILEO
            if (freq_comb=='L1L2') then   ! E1 E5a
                f1=10.23d6*154.d0
                f2=10.23d6*115.0d0
            elseif (freq_comb=='L1L3') then   ! E1 E5b
                f1=10.23d6*154.d0
                f2=10.23d6*118.d0
            elseif (freq_comb=='L2L3') then   ! E5a E5b
                f1=10.23d6*115.0d0
                f2=10.23d6*118.0d0
            end if
        elseif (sys==5) then   ! IRNSS
            f1=10.23d6*115.d0
            f2=10.23d6*243.6d0
        else
            cycle
        end if

        if (Epo_NEQ%Ll1(i)/=0.d0) then
            flag_par_PRN=.true.
            do j=1,npar
                if (iPOS(j)==PRN) then
                    flag_par_PRN=.false.
                    exit
                end if
            end do
!            if (flag_par_PRN) then
!                par_PRNS=par_PRNS+1  ! If partial fixed
!            end if
            if (flag_par_PRN) then ! If partial fixed
!                Epo_NEQ%Al1(i,ParaNum+par_PRNS)=c/f1
!                Epo_NEQ%Ll1(i)=Epo_NEQ%Ll1(i) -  NEQ%amb_L1(PRN)* c/ f1 ! no much difference
                Epo_NEQ%Al1(i,:)=0.d0
                Epo_NEQ%Ll1(i)=0.d0
            else
                if ( (a1*f1+a2*f2/=0.d0) .and. (b1*f1+b2*f2/=0.d0) ) then ! Dual frequency
                    Epo_NEQ%Ll1(i)=Epo_NEQ%Ll1(i) -  NEQ%amb_L1(PRN)* c/f1/sigLC
                else
                    Epo_NEQ%Ll1(i)=Epo_NEQ%Ll1(i) -  NEQ%amb_L1(PRN)* c/(a1*f1+a2*f2)/sigLC
                end if
            end if
            NEQ%Lwl(i)=NEQ%Lwl(i) -  NEQ%amb_WL(PRN)*NEQ%Awl(i,3+PRN)
        end if
        if (Epo_NEQ%Ll2(i)/=0.d0) then
            flag_par_PRN=.true.
            do j=1,npar
                if (iPOS(j)-MaxPRN==PRN) then
                    flag_par_PRN=.false.
                    exit
                end if
            end do
!            if (flag_par_PRN) then
!                par_PRNS=par_PRNS+1  ! If partial fixed
!            end if
            if (flag_par_PRN) then ! If partial fixed
!                Epo_NEQ%Al2(i,ParaNum+par_PRNS)=c/f2
!                Epo_NEQ%Ll2(i)=Epo_NEQ%Ll2(i) -  NEQ%amb_L2(PRN)* c/ f2
                Epo_NEQ%Al2(i,:)=0.d0
                Epo_NEQ%Ll2(i)=0.d0
            else
                if ( (a1*f1+a2*f2/=0.d0) .and. (b1*f1+b2*f2/=0.d0) ) then ! Dual frequency
                    Epo_NEQ%Ll2(i)=Epo_NEQ%Ll2(i) -  NEQ%amb_L2(PRN)* c/ f2/sigLC
                else
                    Epo_NEQ%Ll2(i)=Epo_NEQ%Ll2(i) -  NEQ%amb_L2(PRN)* c/ (b1*f1+b2*f2)/sigLC
                end if
!                write(LogID,'(A6,2X, I2, 2F10.3)') 'L1L2',PRN,Epo_NEQ%Ll1(i),Epo_NEQ%Ll2(i)
            end if
            NEQ%Lw4(i)=NEQ%Lw4(i) -  NEQ%amb_W4(PRN)* NEQ%Aw4(i,3+MaxPRN+PRN)
        end if
    end do
    
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'PRN'
    do i=1,N
        write(unit=LogID,fmt='(I7)',advance='no') Epo_NEQ%PRN(i)
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Ll1'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Ll1(i)*sigLC
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Ll2'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Ll2(i)*sigLC
    end do
    write(unit=LogID,fmt='(A)') ''
    
    Epo_NEQ%Nbb     =  Epo_NEQ%Nbb +  matmul(  matmul( transpose(Epo_NEQ%Al1(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Al1(1:N,:)  )
    Epo_NEQ%Nbb     =  Epo_NEQ%Nbb +  matmul(  matmul( transpose(Epo_NEQ%Al2(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Al2(1:N,:)  )
    Epo_NEQ%U=Epo_NEQ%U+  matmul(  matmul( transpose(Epo_NEQ%Al1(1:N,:)), &
               Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Ll1(1:N) )
    Epo_NEQ%U=Epo_NEQ%U+  matmul(  matmul( transpose(Epo_NEQ%Al2(1:N,:)), &
               Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Ll2(1:N)  )

    ! Start to estimate the coordinate and ambiguities in L1 and L2
    Ad_Flag=.true.
    N=Epo_NEQ%PRNS
    do while(AD_flag)
        Ad_Flag=.false.
        call Invsqrt(Epo_NEQ%Nbb, ParaNum, Epo_NEQ%InvN)
        Epo_NEQ%dx = matmul(Epo_NEQ%InvN, Epo_NEQ%U)
        if (any(isnan(Epo_NEQ%dx))) then
            write(*,*)  "-------Epo_NEQ-------: Epo_NEQ%dx=nan, please check"
            write(LogID, '(5X,A50)')  "-------ERROR-------: Epo_NEQ%dx=nan, please check"
            stop
        end if
        
        ! =================== Outliers Detect =====================
        Epo_NEQ%Vp1(1:N)=(matmul(Epo_NEQ%Ap1(1:N, :), Epo_NEQ%dx(1:3)) - Epo_NEQ%Lp1(1:N))
        Epo_NEQ%Vp2(1:N)=(matmul(Epo_NEQ%Ap2(1:N, :), Epo_NEQ%dx(1:3)) - Epo_NEQ%Lp2(1:N))
        Epo_NEQ%Vl1(1:N)=(matmul(Epo_NEQ%Al1(1:N, :), Epo_NEQ%dx) - Epo_NEQ%Ll1(1:N))
        Epo_NEQ%Vl2(1:N)=(matmul(Epo_NEQ%Al2(1:N, :), Epo_NEQ%dx) - Epo_NEQ%Ll2(1:N))
        do i=1,N
             Epo_NEQ%Vp1(i)= Epo_NEQ%Vp1(i)*sigLC*sqrt(2.d0*Epo_NEQ%P(i,i))   ! unify to  carrier phase magnitude and the same weight, max P(i,i) is 0.5
             Epo_NEQ%Vp2(i)= Epo_NEQ%Vp2(i)*sigLC*sqrt(2.d0*Epo_NEQ%P(i,i))
             Epo_NEQ%Vl1(i)= Epo_NEQ%Vl1(i)*sigLC*sqrt(2.d0*Epo_NEQ%P(i,i))
             Epo_NEQ%Vl2(i)= Epo_NEQ%Vl2(i)*sigLC*sqrt(2.d0*Epo_NEQ%P(i,i))
        end do
        Epo_NEQ%maxV(1:1)=maxval(dabs(Epo_NEQ%Vp1(1:N)))
        Epo_NEQ%maxL(1:1)=maxloc(dabs(Epo_NEQ%Vp1(1:N)))
        Epo_NEQ%maxV(2:2)=maxval(dabs(Epo_NEQ%Vp2(1:N)))
        Epo_NEQ%maxL(2:2)=maxloc(dabs(Epo_NEQ%Vp2(1:N)))
        Epo_NEQ%maxV(3:3)=maxval(dabs(Epo_NEQ%Vl1(1:N)))
        Epo_NEQ%maxL(3:3)=maxloc(dabs(Epo_NEQ%Vl1(1:N)))
        Epo_NEQ%maxV(4:4)=maxval(dabs(Epo_NEQ%Vl2(1:N)))
        Epo_NEQ%maxL(4:4)=maxloc(dabs(Epo_NEQ%Vl2(1:N)))

        maxV=maxval(dabs(Epo_NEQ%maxV))
        maxL=maxloc(dabs(Epo_NEQ%maxV),dim=1)
        
        if ( dabs(maxV)>0.04d0 ) then
            if ((maxL)==1) then   ! maxV in P1
                call Minus_NEQ( Epo_NEQ%Nbb(1:ParaNum,1:ParaNum), Epo_NEQ%U(1:ParaNum), Epo_NEQ%Ap1(1:N,:), Epo_NEQ%Lp1(1:N), &
                       Epo_NEQ%P(1:N, 1:N), ParaNum,  N, Epo_NEQ%maxL(1), Epo_NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in P1', 'PRN=',Epo_NEQ%PRN(Epo_NEQ%maxL(1)),'maxV=',maxV
            elseif  ((maxL)==2) then   ! maxV in P2
                call Minus_NEQ( Epo_NEQ%Nbb(1:ParaNum,1:ParaNum), Epo_NEQ%U(1:ParaNum), Epo_NEQ%Ap2(1:N,:), Epo_NEQ%Lp2(1:N), &
                       Epo_NEQ%P(1:N, 1:N), ParaNum,  N, Epo_NEQ%maxL(2), Epo_NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in P2', 'PRN=',Epo_NEQ%PRN(Epo_NEQ%maxL(2)),'maxV=',maxV
            elseif  ((maxL)==3) then
                call Minus_NEQ( Epo_NEQ%Nbb(1:ParaNum,1:ParaNum), Epo_NEQ%U(1:ParaNum), Epo_NEQ%Al1(1:N,:), Epo_NEQ%Ll1(1:N), &
                       Epo_NEQ%P(1:N, 1:N), ParaNum,  N, Epo_NEQ%maxL(3), Epo_NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in L1', 'PRN=',Epo_NEQ%PRN(Epo_NEQ%maxL(3)),'maxV=',maxV
                if (ar_mode==3) then ! If fixed and hold mode
                    NEQ%fixed_amb(Epo_NEQ%PRN(Epo_NEQ%maxL(3)))=0.99d0  ! If outliers, do not hold the amiguity
                    NEQ%fixed_amb_num(Epo_NEQ%PRN(Epo_NEQ%maxL(3)))=0
                end if
!                if (Var_smooth=='y' .or. Var_smooth=='Y') then ! If smooth pseudorange
!                    STA%STA(1)%Pre(Epo_NEQ%PRN(Epo_NEQ%maxL(3)))%n=1.d0
!                    STA%STA(2)%Pre(Epo_NEQ%PRN(Epo_NEQ%maxL(3)))%n=1.d0 
!                end if
            elseif  ((maxL)==4) then
                call Minus_NEQ( Epo_NEQ%Nbb(1:ParaNum,1:ParaNum), Epo_NEQ%U(1:ParaNum), Epo_NEQ%Al2(1:N,:), Epo_NEQ%Ll2(1:N), &
                       Epo_NEQ%P(1:N, 1:N), ParaNum,  N, Epo_NEQ%maxL(4), Epo_NEQ%SumN)
                write(LogID,'(10X,A14,A5,I3, A6, F10.3)') 'outlier in L2', 'PRN=',Epo_NEQ%PRN(Epo_NEQ%maxL(4)),'maxV=',maxV
                if (ar_mode==3) then ! If fixed and hold mode
                    NEQ%fixed_amb(Epo_NEQ%PRN(Epo_NEQ%maxL(4))+MaxPRN)=0.99d0  ! If outliers, do not hold the amiguity
                    NEQ%fixed_amb_num(Epo_NEQ%PRN(Epo_NEQ%maxL(4))+MaxPRN)=0
                end if
!                if (Var_smooth=='y' .or. Var_smooth=='Y') then ! If smooth pseudorange
!                    STA%STA(1)%Pre(Epo_NEQ%PRN(Epo_NEQ%maxL(4)))%n=1.d0
!                    STA%STA(2)%Pre(Epo_NEQ%PRN(Epo_NEQ%maxL(4)))%n=1.d0 
!                end if
            end if

            Ad_Flag=.true.
        end if
        ! =================== End of Outliers Detect =====================
    end do
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'PRN'
    do i=1,N
        write(unit=LogID,fmt='(I7)',advance='no') Epo_NEQ%PRN(i)
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vp1'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Vp1(i)*sigPC/sigLC/sqrt(2.d0*Epo_NEQ%P(i,i))    ! back to raw residuals
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vp2'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Vp2(i)*sigPC/sigLC/sqrt(2.d0*Epo_NEQ%P(i,i))
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vl1'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Vl1(i)/sqrt(2.d0*Epo_NEQ%P(i,i))
    end do
    write(unit=LogID,fmt='(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'Vl2'
    do i=1,N
        write(unit=LogID,fmt='(F7.3)',advance='no') Epo_NEQ%Vl2(i)/sqrt(2.d0*Epo_NEQ%P(i,i))
    end do
    write(unit=LogID,fmt='(A)') ''

    Coor=Epo_NEQ%dx(1:3)
    write(unit=LogID,fmt='(A10,3F7.2)') 'fixcoor',Epo_NEQ%dx(1:3)
    return
    
end subroutine