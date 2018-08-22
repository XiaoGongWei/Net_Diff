! ==================== Form NEQ ====================
! PURPOSE:
!             To form the normal equation information at this epoch
!
! INPUTS:
!              SD                        station difference structure
!              DD                       double difference structure

! OUTPUTS:
!              NEQ                     normal equation structure, used in 
!                                              WL and W4 ambiguity solution
!             Epo_NEQ              current epoch normal equation structure, 
!                                              used in coordinate estimation after fix 
!                                              the L1 and L2 ambiguity
! 
! WRITTEN BY: Yize Zhang
! ================================================

subroutine Form_NEQ(SD, DD, NEQ, Epo_NEQ)
use MOD_SD
use MOD_DD
use MOD_NEQ
use MOD_Epo_NEQ
use MOD_constant
use MOD_Var
use MOD_CycleSlip
use MOD_GLO_Fre
use MOD_FileID
use MOD_NEQ_DP
implicit none
    type(type_SD)   ::  SD
    type(type_DD)   ::  DD
    type(type_NEQ) ::  NEQ
    type(type_Epo_NEQ) :: Epo_NEQ
    ! Local variables
    integer :: N, i, j, k, PRN, sys, freq
    real(8)  :: Awl(MaxPRN, ParaNum+2*MaxPRN+IonoNum),Aw4(MaxPRN, ParaNum+2*MaxPRN+IonoNum)
    real(8) :: Aewl(MaxPRN,ParaNum)
    real(8)  :: Ap1(MaxPRN, ParaNum+3*IonoNum), Ap2(MaxPRN, ParaNum+3*IonoNum)
    real(8)  :: Al1(MaxPRN, ParaNum+3*IonoNum), Al2(MaxPRN, ParaNum+3*IonoNum)
    logical :: flag_del_PRN
    real(8) :: dT(2), factor
    real(8) :: Kk(ParaNum+3*IonoNum,1)  ! Kalman gain

    
    Awl=0.d0; Aw4=0.d0
    Aewl=0.d0
    Al1=0.d0;  Al2=0.d0
    Ap1=0.d0; Ap2=0.d0
    NEQ%Lp1=0.d0; NEQ%Lp2=0.d0
    NEQ%Lwl=0.d0; NEQ%Lw4=0.d0

    if (ar_mode==2) then ! Instantaneous AR
        NEQ%Nbb=0.d0
        NEQ%InvN=0.d0
        NEQ%U=0.d0
        NEQ%dx=0.d0
        if (If_Est_Iono .and. IonoNum>0) then
            Epo_NEQ%Nbb=0.d0
            Epo_NEQ%InvN=0.d0
            Epo_NEQ%U=0.d0
            Epo_NEQ%dx=0.d0
        end if
    else
        ! Eliminate the unobserved satellite at this epoch
        do PRN=1, MaxPRN !DD%PRNS
            flag_del_PRN=.true.
            do i=1,DD%PRNS
                if (DD%PRN(i)==PRN) then
                    flag_del_PRN=.false.   ! If this satellite exist in this epoch, do not delete
                    exit
                end if
            end do
            dT(1)=(CycleSlip(1)%CS(PRN)%WeekPrev(3)-DD%Week)*604800.0d0+CycleSlip(1)%CS(PRN)%SowPrev(3)-DD%Sow
            dT(2)=(CycleSlip(2)%CS(PRN)%WeekPrev(3)-DD%Week)*604800.0d0+CycleSlip(2)%CS(PRN)%SowPrev(3)-DD%Sow
            ! First is the ionosphere parameter
            if (If_Est_Iono .and. IonoNum>0) then
                if (flag_del_PRN .and. any(dT<-600.d0)) then ! If not observed more than 10 minutes
!                    if (ADmethod=='LS') then
!                            call Elimi_Para(Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%N, ParaNum+2*IonoNum+PRN) !  ionosphere parameter
!                    elseif (ADmethod=='KF') then  ! As ionosphere parameter is estimated as randon walk, so we use transformed Kalman Filter instead of Least Square
                            call KF_Change(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, ParaNum+2*IonoNum+PRN, 'dds')  ! ionosphere parameter
!                    end if
                end if
            end if
            ! Then comes to ambiguity parameter
            if (ADmethod=='LS' .and. abs(NEQ%Nbb(ParaNum+PRN,ParaNum+PRN))<1.d-11) then
                cycle
            elseif (ADmethod=='KF' .and. abs(NEQ%InvN(ParaNum+PRN,ParaNum+PRN))<1.d-11) then
                cycle
            end if
            if (flag_del_PRN .and. any(dT<-3.d0*Interval)) then   !  If satellite unobserved more than 3 epoches. result shows not very good, so don't use.
                if (ADmethod=='LS') then
                    call Elimi_Para(NEQ%Nbb, NEQ%U, NEQ%N, ParaNum+PRN)
                    call Elimi_Para(NEQ%Nbb, NEQ%U, NEQ%N, ParaNum+MaxPRN+PRN)
                    NEQ%dx(ParaNum+PRN)=0.d0
                    NEQ%dx(ParaNum+MaxPRN+PRN)=0.d0
                elseif (ADmethod=='KF') then
                    call KF_Change(NEQ%InvN, NEQ%dx,NEQ%N, ParaNum+PRN, 'dds')
                    call KF_Change(NEQ%InvN, NEQ%dx,NEQ%N, ParaNum+MaxPRN+PRN, 'dds')
                end if
                if (If_Est_Iono .and. IonoNum>0) then
!                    if (ADmethod=='LS') then
!                        call Elimi_Para(Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%N, ParaNum+PRN)  ! L1 ambiguity
!                        call Elimi_Para(Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%N, ParaNum+IonoNum+PRN) ! L2 ambiguity
!                        Epo_NEQ%dx(ParaNum+PRN)=0.d0
!                        Epo_NEQ%dx(ParaNum+IonoNum+PRN)=0.d0
!                    elseif (ADmethod=='KF') then ! As ionosphere parameter is estimated as randon walk, so we use transformed Kalman Filter instead of Least Square
                        call KF_Change(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, ParaNum+PRN, 'dds')
                        call KF_Change(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, ParaNum+IonoNum+PRN, 'dds')
!                    end if
                end if
                if (ar_mode==3) then ! If fixed and hold mode
                    NEQ%fixed_amb(PRN)=0.99d0   ! Reinitialize the fixed ambiguity
                    NEQ%fixed_amb(PRN+MaxPRN)=0.99d0
                    NEQ%fixed_amb_num(PRN)=0
                    NEQ%fixed_amb_num(PRN+MaxPRN)=0
                    if (If_Est_Iono .and. IonoNum>0) then
                        Epo_NEQ%fixed_amb(PRN)=0.99d0   ! Reinitialize the fixed ambiguity
                        Epo_NEQ%fixed_amb(PRN+MaxPRN)=0.99d0
                        Epo_NEQ%fixed_amb_num(PRN)=0
                        Epo_NEQ%fixed_amb_num(PRN+MaxPRN)=0
                    end if
                end if
            end if
        end do
    end if
    if (Pos_State=="K") then    ! on the fly, kinematic
        if (ADmethod=='LS') then
            call Elimi_Para(NEQ%Nbb, NEQ%U, NEQ%N, 1)
            call Elimi_Para(NEQ%Nbb, NEQ%U, NEQ%N, 2)
            call Elimi_Para(NEQ%Nbb, NEQ%U, NEQ%N, 3)
        elseif (ADmethod=='KF') then
            call KF_Change(NEQ%InvN, NEQ%dx,NEQ%N, 1, 'ddp')
            call KF_Change(NEQ%InvN, NEQ%dx,NEQ%N, 2, 'ddp')
            call KF_Change(NEQ%InvN, NEQ%dx,NEQ%N, 3, 'ddp')
        end if
        if (If_Est_Iono .and. IonoNum>0) then  ! When estimate ionosphere parameter
            call KF_Change(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, 1, 'ddp')
            call KF_Change(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, 2, 'ddp')
            call KF_Change(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, 3, 'ddp')
        else ! If not estimate ionosphere, Least square is enough for Epo_NEQ as it has only coordinate paramters
            call Elimi_Para(Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%N, 1)
            call Elimi_Para(Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%N, 2)
            call Elimi_Para(Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%N, 3)
        end if
    elseif (Pos_State=="S") then
        if (ADmethod=='KF') then
            if (NEQ%InvN(1,1)==0.d0) then ! If first epoch
                call KF_Change(NEQ%InvN, NEQ%dx,NEQ%N, 1, 'ddp')
                call KF_Change(NEQ%InvN, NEQ%dx,NEQ%N, 2, 'ddp')
                call KF_Change(NEQ%InvN, NEQ%dx,NEQ%N, 3, 'ddp')
            end if
        end if
        if (Epo_NEQ%InvN(1,1)==0.d0 .and. If_Est_Iono .and. IonoNum>0) then  ! When estimate ionosphere parameter
            call KF_Change(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, 1, 'ddp')
            call KF_Change(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, 2, 'ddp')
            call KF_Change(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, 3, 'ddp')
        end if
    elseif (Pos_State=="F") then
        DD%A(:,1:3)=0.d0
    end if
    Epo_NEQ%Al1=0.d0
    Epo_NEQ%Al2=0.d0
    Epo_NEQ%Ll1=0.d0
    Epo_NEQ%Ll2=0.d0
    Epo_NEQ%Vl1=0.d0
    Epo_NEQ%Vl2=0.d0

    N=DD%PRNS
    NEQ%SumN=0
    Epo_NEQ%SumN=0

    do i=1, DD%PRNS
        PRN=DD%PRN(i)
        sys=DD%Sys(i)
        
        if (sys==1) then   ! GPS/QZSS
            f1=10.23d6*154d0
            f2=10.23d6*120d0
            f3=10.23d6*115d0
        elseif (sys==2) then   ! GLONASS
            freq=Fre_Chann(PRN-GNum)
            f1=(1602.0d0+freq*0.5625d0)*1.0D6   ! f1=(1602.0d0+K*0.5625d0)*1.0d6
            f2=(1246.0d0+freq*0.4375d0)*1.0D6
        elseif  (sys==3) then   ! COMPASS
            if (freq_comb=='L1L2') then
                f1=10.23d6*152.6d0
                f2=10.23d6*118.0d0
                f3=10.23d6*124.0d0
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
                f3=10.23d6*118.0d0
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

        if (ar_mode/=2) then ! If not instantaneous AR
            if ( (CycleSlip(1)%Slip(PRN)==1) .or. (CycleSlip(2)%Slip(PRN)==1) ) then ! Cycle Slip in 
                 write(LogID,'(A10,I3,A11)') 'PRN', PRN,'cycle slip'  
                 if (ADmethod=='LS') then
                     call Elimi_Para(NEQ%Nbb, NEQ%U, NEQ%N, ParaNum+PRN)
                     call Elimi_Para(NEQ%Nbb, NEQ%U, NEQ%N, ParaNum+MaxPRN+PRN)
                 elseif (ADmethod=='KF') then
                    call KF_Change(NEQ%InvN, NEQ%dx,NEQ%N, ParaNum+PRN, 'dda')  ! L1 ambiguity
                    call KF_Change(NEQ%InvN, NEQ%dx,NEQ%N, ParaNum+MaxPRN+PRN, 'dda')  ! L2 ambiguity
                 end if
                 if (If_Est_Iono .and. IonoNum>0) then
!                    if (ADmethod=='LS') then
!                         call Elimi_Para(Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%N, ParaNum+PRN)  ! L1 ambiguity
!                         call Elimi_Para(Epo_NEQ%Nbb, Epo_NEQ%U, Epo_NEQ%N, ParaNum+IonoNum+PRN) ! L2 ambiguity
!                    elseif (ADmethod=='KF') then ! As ionosphere parameter is estimated as randon walk, so we use transformed Kalman Filter instead of Least Square
                        call KF_Change(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, ParaNum+PRN, 'dda')  ! L1 ambiguity
                        call KF_Change(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, ParaNum+IonoNum+PRN, 'dda')  ! L2 ambiguity
!                    end if
                 end if
                 if (ar_mode==3) then ! If fixed and hold mode
                    NEQ%fixed_amb(PRN)=0.99d0   ! Reinitialize the fixed ambiguity if cycle slip occurs
                    NEQ%fixed_amb(PRN+MaxPRN)=0.99d0
                    NEQ%fixed_amb_num(PRN)=0
                    NEQ%fixed_amb_num(PRN+MaxPRN)=0
                    if (If_Est_Iono .and. IonoNum>0) then
                        Epo_NEQ%fixed_amb(PRN)=0.99d0   ! Reinitialize the fixed ambiguity if cycle slip occurs
                        Epo_NEQ%fixed_amb(PRN+MaxPRN)=0.99d0
                        Epo_NEQ%fixed_amb_num(PRN)=0
                        Epo_NEQ%fixed_amb_num(PRN+MaxPRN)=0
                    end if
                 end if
            end if
        end if

        ! 将信息加入到法方程中
         ! For P1
         if ( DD%P1(i)/=0.d0 ) then
             Ap1(i,1:ParaNum)=DD%A(i,1:ParaNum)
             if (If_Est_Iono .and. IonoNum>0) then 
                Ap1(i,ParaNum+IonoNum*2+PRN)=1.d0  ! Ionosphere parameter
                if (If_IonoCompensate .and. Epo_NEQ%ratio>minratio .and. EPO_NEQ%fixed_amb_num(PRN)>10 .and. DD%Ele(PRN)>FixEle) then ! if EPO_NEQ fix and hold ambiguity
                    NEQ%Lp1(i)=(DD%P1(i) - Epo_NEQ%dx(ParaNum+IonoNum*2+PRN))/sigPC
                else
                    NEQ%Lp1(i)=DD%P1(i)/sigPC
                end if
             else  ! Normal double difference
                NEQ%Lp1(i)=DD%P1(i)/sigPC
             end if
             NEQ%SumN=NEQ%SumN+1
             Epo_NEQ%SumN=Epo_NEQ%SumN+1
         end if

         ! For P2
         if ( DD%P2(i)/=0.d0 ) then
             Ap2(i,1:ParaNum)=DD%A(i,1:ParaNum)
             if (If_Est_Iono .and. IonoNum>0) then 
                Ap2(i,ParaNum+IonoNum*2+PRN)=f1**2/f2**2  ! Ionosphere parameter
                if (If_IonoCompensate .and. Epo_NEQ%ratio>minratio .and. EPO_NEQ%fixed_amb_num(PRN)>5 .and. DD%Ele(PRN)>FixEle) then ! if EPO_NEQ fix and hold ambiguity
                    NEQ%Lp2(i)=(DD%P2(i) - f1**2/f2**2*Epo_NEQ%dx(ParaNum+IonoNum*2+PRN))/sigPC
                else
                    NEQ%Lp2(i)=DD%P2(i)/sigPC
                end if
             else  ! Normal double difference
                NEQ%Lp2(i)=DD%P2(i)/sigPC
             end if
             NEQ%SumN=NEQ%SumN+1
             Epo_NEQ%SumN=Epo_NEQ%SumN+1
         end if

         ! For WL
         if ( DD%WL(i)/=0.d0 ) then
             if (DD%P1(i)/=0.d0 .or. DD%P2(i)/=0.d0 .or. NEQ%InvN(ParaNum+PRN, ParaNum+PRN)/=0.d0) then
                 Awl(i,1:ParaNum)=DD%A(i,1:ParaNum)   ! For first epoch, if no psedo-range, NEQ can't be solved
                 Awl(i,ParaNum+PRN)=c/(a1*f1+a2*f2)  ! 1.d0
                 NEQ%SumN=NEQ%SumN+1
                 if (If_Est_Iono .and. IonoNum>0) then
                    Awl(i,ParaNum+IonoNum*2+PRN)= f1/f2 ! Ionosphere parameter
                    if (If_IonoCompensate .and. Epo_NEQ%ratio>minratio .and. EPO_NEQ%fixed_amb_num(PRN)>5 .and. DD%Ele(PRN)>FixEle) then ! if EPO_NEQ fix and hold ambiguity
                        NEQ%Lwl(i)=(DD%WL(i) - f1/f2*Epo_NEQ%dx(ParaNum+IonoNum*2+PRN))/sigLC/sqrt(a1**2+a2**2)  ! should be sqrt((a1*f1)**2+(a2*f2)**2)/(a1*f1-a2*f2)
                    else
                        NEQ%Lwl(i)=DD%WL(i)/sigLC/sqrt(a1**2+a2**2)  ! should be sqrt((a1*f1)**2+(a2*f2)**2)/(a1*f1-a2*f2)
                    end if
                    Epo_NEQ%SumN=Epo_NEQ%SumN+1
                 else  ! Normal double difference
                    NEQ%Lwl(i)=DD%WL(i)/sigLC/sqrt(a1**2+a2**2)  ! should be sqrt((a1*f1)**2+(a2*f2)**2)/(a1*f1-a2*f2)
                 end if
             else
                 DD%WL(i)=0.d0
             end if
         end if

         ! For W4
         if ( DD%W4(i)/=0.d0 ) then
             if (DD%P1(i)/=0.d0 .or. DD%P2(i)/=0.d0 .or. NEQ%InvN(ParaNum+MaxPRN+PRN, ParaNum+MaxPRN+PRN)/=0.d0) then
                 Aw4(i,1:ParaNum)=DD%A(i,1:ParaNum)   ! For first epoch, if no psedo-range, NEQ can't be solved
                 Aw4(i,ParaNum+MaxPRN+PRN)=c/(b1*f1+b2*f2)  !1.d0
                 NEQ%SumN=NEQ%SumN+1
                 if (If_Est_Iono .and. IonoNum>0) then
                    Aw4(i,ParaNum+MaxPRN+PRN)=c/(f1-f3)    ! L1-L3 wide lane ambiguity parameter
                    Aw4(i,ParaNum+IonoNum*2+PRN)= f1/f3 ! Ionosphere parameter
                    Epo_NEQ%SumN=Epo_NEQ%SumN+1
                    NEQ%SumN=NEQ%SumN-1  ! W4 not used in NEQ for long baseline
                 end if
             else
                DD%W4(i)=0.d0
             end if
         end if

         ! For EWL
         if ( DD%EWL(i)/=0.d0 ) then
             Aewl(i,1:ParaNum)=DD%A(i,1:ParaNum)
             NEQ%SumN=NEQ%SumN+1
         end if

         ! For L1
         if ( DD%L1(i)/=0.d0 ) then
             Al1(i,1:ParaNum)=DD%A(i,1:ParaNum)
             if (If_Est_Iono .and. IonoNum>0) then     ! Danger: For first epoch, if no psedo-range or WL help, NEQ can't be solved
                Al1(i,ParaNum+PRN)= c/f1   ! ambiguity parameter 
                Al1(i,ParaNum+2*IonoNum+PRN)= -1.d0 ! Ionosphere parameter
             end if
             Epo_NEQ%SumN=Epo_NEQ%SumN+1
         end if

         ! For L2
         if ( DD%L2(i)/=0.d0 ) then
             Al2(i,1:ParaNum)=DD%A(i,1:ParaNum)
             if (If_Est_Iono .and. IonoNum>0) then 
                Al2(i,ParaNum+MaxPRN+PRN)= c/f2 ! ambiguity parameter
                Al2(i,ParaNum+2*IonoNum+PRN)=  -f1**2/f2**2  ! Ionosphere parameter
             end if
             Epo_NEQ%SumN=Epo_NEQ%SumN+1
         end if
    end do

    NEQ%PRNS=DD%PRNS
    NEQ%PRN=DD%PRN
    NEQ%Ele=DD%Ele
    NEQ%Sys=DD%Sys
    NEQ%System=DD%System
    NEQ%P(1:N,1:N)=DD%P(1:N,1:N)
    NEQ%R(1:N,1:N)=DD%Q(1:N,1:N)
    NEQ%Ap1(1:N,:)=Ap1(1:N,1:ParaNum)/sigPC      ! P=P*0.01d0
    NEQ%Ap2(1:N,:)=Ap2(1:N,1:ParaNum)/sigPC     ! P=P*0.01d0
!    NEQ%Lp1(1:N)=DD%P1(1:N)/sigPC         ! Already done in former loop
!    NEQ%Lp2(1:N)=DD%P2(1:N)/sigPC         ! Already done in former loop
    NEQ%Aewl(1:N,:)=Aewl(1:N,:)/sigLC/30.d0     ! Extra Wide Lane P=P/(sqrt((a1*f1)**2+(a2*f2)**2)/(a1*f1-a2*f2)), for GPS 33.2, for BeiDou 28.5, for Galileo 54.9
    NEQ%Lewl(1:N)=DD%EWL(1:N)/sigLC/30.d0 
    NEQ%amb_EWL=DD%EWL_amb
    if ( (a1/=0.d0) .or. (a2/=0.d0) ) then
        NEQ%Awl(1:N,:)=Awl(1:N,1:ParaNum+MaxPRN*2)/sigLC/sqrt(a1**2+a2**2)
!        NEQ%Lwl(1:N)=DD%WL(1:N)/sigLC/sqrt(a1**2+a2**2)         ! Already done in former loop
    end if
    if ( (b1/=0.d0) .or. (b2/=0.d0) ) then
        NEQ%Aw4(1:N,:)=Aw4(1:N,1:ParaNum+MaxPRN*2)/sigLC/sqrt(b1**2+b2**2)
        NEQ%Lw4(1:N)=DD%W4(1:N)/sigLC/sqrt(b1**2+b2**2)
    end if
    if (ParaNum==4) then    ! Don't estimate tropsphere parameter in NEQ Only in EPO_NEQ for long baseline need to estimate it.
        NEQ%Ap1(1:N,4)=0.d0
        NEQ%Ap2(1:N,4)=0.d0
        NEQ%Aewl(1:N,4)=0.d0
        NEQ%Awl(1:N,4)=0.d0
        NEQ%Aw4(1:N,4)=0.d0
    end if
    
    if (GloParaNum>0) then  ! If GLONASS
        do i=1,GloParaNum
            if (NEQ%InvN(ParaNum-GloParaNum+i,ParaNum-GloParaNum+i)==0.d0) then
                NEQ%InvN(ParaNum-GloParaNum+i,ParaNum-GloParaNum+i)=0.5d0**2  ! 1m
            end if
        end do
    end if
    ! Fix and Hold mode, add constraints to fixed ambiguity, See RTKLIB rtkpos.c  holdamb()
    if (ar_mode==3) then !
        do PRN=1, 2*maxPRN
            if (NEQ%fixed_amb(PRN)/=0.99d0) then
                NEQ%fixed_amb_num(PRN)=NEQ%fixed_amb_num(PRN)+1
                if (PRN<=maxPRN) then
                    if (NEQ%fixed_amb_num(PRN)>10 .and. NEQ%Ele(PRN)>HoldEle) then  ! If the same ambiguity >10 epoches
                        if (ADmethod=='LS') then
                            NEQ%Nbb(PRN+ParaNum,PRN+ParaNum)=NEQ%Nbb(PRN+ParaNum,PRN+ParaNum)+1.d0*(1.d0/0.01**2) ! 1cm
                            NEQ%U(PRN+ParaNum)=NEQ%U(PRN+ParaNum)+real(NEQ%fixed_amb(PRN))*(1.d0/0.01**2) ! 1cm
                        elseif (ADmethod=='KF') then
                            call KF_Gain_one(NEQ%InvN, NEQ%dx,NEQ%N, PRN+ParaNum, 1.d0*(NEQ%fixed_amb(PRN)), 0.01d0)
                        end if
                    end if
                elseif (PRN>maxPRN) then
                    if (NEQ%fixed_amb_num(PRN)>10 .and. NEQ%Ele(PRN-maxPRN)>HoldEle) then  ! If the same ambiguity >10 epoches
!                        if (ADmethod=='LS') then
                            NEQ%Nbb(PRN+ParaNum,PRN+ParaNum)=NEQ%Nbb(PRN+ParaNum,PRN+ParaNum)+1.d0*(1/0.01**2) ! 1cm
                            NEQ%U(PRN+ParaNum)=NEQ%U(PRN+ParaNum)+real(NEQ%fixed_amb(PRN))*(1.d0/0.01**2) ! 1cm
!                        elseif (ADmethod=='KF') then
!                            call KF_Gain_one(NEQ%InvN, NEQ%dx,NEQ%N, PRN+ParaNum, 1.d0*(NEQ%fixed_amb(PRN)), 0.01d0)
!                        end if
                    end if
                end if
            end if
        end do
    end if

    ! If Doppler velocity used, add constraints to fixed coordinate
    if (Combination(3) .and. Vel_Used==1) then
        do i=1,3
            if (NEQ_DP%Flag_Sln(5)==1) then ! If full fixed
                if (ADmethod=='LS') then
                    NEQ%Nbb(i,i)=NEQ%Nbb(i,i)+1.d0/0.01d0/NEQ_DP%dt  ! This is wrong for Nbb
                elseif (ADmethod=='KF') then
                    NEQ%InvN(i,i)=NEQ%InvN(i,i)+0.01d0*NEQ_DP%dt    ! Precision of estimated position using doppler velocity is 0.1m, (1/0.1**2) ! 0.1m
                end if
            elseif (NEQ_DP%Flag_Sln(5)==2) then ! If partial fixed
                if (ADmethod=='LS') then
                    NEQ%Nbb(i,i)=NEQ%Nbb(i,i)+1.d0/0.25d0*NEQ_DP%dt  ! This is wrong for Nbb
                elseif (ADmethod=='KF') then
                    NEQ%InvN(i,i)=NEQ%InvN(i,i)+0.25d0*NEQ_DP%dt     ! (1/0.5**2) ! 0.5m
                end if
            elseif (NEQ_DP%Flag_Sln(5)==3) then ! If not fixed
                if (ADmethod=='LS') then
                    NEQ%Nbb(i,i)=NEQ%Nbb(i,i)+1.d0/4.d0*NEQ_DP%dt   ! This is wrong for Nbb
                elseif (ADmethod=='KF') then
                    NEQ%InvN(i,i)=NEQ%InvN(i,i)+4.d0*NEQ_DP%dt    ! (1/2.d0**2) ! 2m
                end if
            end if
        end do
    end if
    
    if (ADmethod=='LS') then     ! When Least square, add Nbb
        NEQ%Nbb(1:ParaNum,1:ParaNum) = NEQ%Nbb(1:ParaNum,1:ParaNum) +  matmul(  matmul( transpose(NEQ%Ap1(1:N,:)), NEQ%P(1:N,1:N) ), NEQ%Ap1(1:N,:)  )
        NEQ%Nbb(1:ParaNum,1:ParaNum) = NEQ%Nbb(1:ParaNum,1:ParaNum) +  matmul(  matmul( transpose(NEQ%Ap2(1:N,:)), NEQ%P(1:N,1:N) ), NEQ%Ap2(1:N,:)  )
        NEQ%Nbb(1:ParaNum,1:ParaNum) = NEQ%Nbb(1:ParaNum,1:ParaNum) +  matmul(  matmul( transpose(NEQ%Aewl(1:N,:)), NEQ%P(1:N,1:N) ), NEQ%Aewl(1:N,:)  )  ! EWL
        NEQ%Nbb             =  NEQ%Nbb               +  matmul(  matmul( transpose(NEQ%Awl(1:N,:)), NEQ%P(1:N,1:N) ), NEQ%Awl(1:N,:)  )
        NEQ%Nbb             =  NEQ%Nbb               +  matmul(  matmul( transpose(NEQ%Aw4(1:N,:)), NEQ%P(1:N,1:N) ), NEQ%Aw4(1:N,:)  )

        NEQ%U(1:ParaNum)          =  NEQ%U(1:ParaNum)           +  matmul(  matmul( transpose(NEQ%Ap1(1:N,:)), NEQ%P(1:N,1:N) ), NEQ%Lp1(1:N)  )
        NEQ%U(1:ParaNum)          =  NEQ%U(1:ParaNum)           +  matmul(  matmul( transpose(NEQ%Ap2(1:N,:)), NEQ%P(1:N,1:N) ), NEQ%Lp2(1:N)  )
        NEQ%U(1:ParaNum)          =  NEQ%U(1:ParaNum)           +  matmul(  matmul( transpose(NEQ%Aewl(1:N,:)), NEQ%P(1:N,1:N) ), NEQ%Lewl(1:N)  )  ! EWL
        NEQ%U                 =  NEQ%U                   +  matmul(  matmul( transpose(NEQ%Awl(1:N,:)), NEQ%P(1:N,1:N) ), NEQ%Lwl(1:N)  )
        NEQ%U                 =  NEQ%U                   +  matmul(  matmul( transpose(NEQ%Aw4(1:N,:)), NEQ%P(1:N,1:N) ), NEQ%Lw4(1:N)  )
    elseif (ADmethod=='KF') then
        ! Kalmen Filter will done in Solve_NEQ_Iono, or Solve_NEQ,  KF_Gain function
    end if

    Epo_NEQ%PRNS=DD%PRNS
    Epo_NEQ%PRN=DD%PRN
    Epo_NEQ%Ele=DD%Ele
    Epo_NEQ%Sys=DD%Sys
    Epo_NEQ%System=DD%System
    Epo_NEQ%P(1:N,1:N)=DD%P(1:N,1:N)
    Epo_NEQ%R(1:N,1:N)=DD%Q(1:N,1:N)
    Epo_NEQ%Ap1(1:N,:)=Ap1(1:N,:)/sigPC
    Epo_NEQ%Lp1(1:N)=DD%P1(1:N)/sigPC
    Epo_NEQ%Ap2(1:N,:)=Ap2(1:N,:)/sigPC
    Epo_NEQ%Lp2(1:N)=DD%P2(1:N)/sigPC
    Epo_NEQ%Al1(1:N,:)=Al1(1:N,:)/sigLC
    Epo_NEQ%Ll1(1:N)=DD%L1(1:N)/sigLC
    Epo_NEQ%Al2(1:N,:)=Al2(1:N,:)/sigLC
    Epo_NEQ%Ll2(1:N)=DD%L2(1:N)/sigLC
    Epo_NEQ%amb_WL=DD%WL_amb  ! Just for test, not very good, because of the wrong rounding integer due to code multipath
    if (GloParaNum>0) then  ! If GLONASS
        do i=1,GloParaNum
            if (Epo_NEQ%InvN(ParaNum-GloParaNum+i,ParaNum-GloParaNum+i)==0.d0) then
                Epo_NEQ%InvN(ParaNum-GloParaNum+i,ParaNum-GloParaNum+i)=0.5d0**2  ! 1m
            end if
        end do
    end if
    if (If_Est_Iono .and. IonoNum>0) then 
        if ( (a1/=0.d0) .or. (a2/=0.d0) ) then
            Epo_NEQ%Awl(1:N,:)=Awl(1:N,:)/sigLC/2.d0 ! sqrt(a1**2+a2**2)   ! The observation noise of WL may be greater 
            Epo_NEQ%Lwl(1:N)=DD%WL(1:N)/sigLC/2.d0 ! sqrt(a1**2+a2**2)   ! should be sqrt((a1*f1)**2+(a2*f2)**2)/(a1*f1-a2*f2), for GPS, it is 5.74, for BeiDou, it is 5.57, for Galileo, it is 4.93
        end if
!        if ( (b1/=0.d0) .or. (b2/=0.d0) ) then  ! only for triple frequency
            Epo_NEQ%Aw4(1:N,:)=Aw4(1:N,:)/sigLC/2.d0 ! sqrt(2.d0)  ! 6.9d0 ! should be sqrt((b1*f1)**2+(b2*f2)**2)/(b1*f1-b2*f2)
            Epo_NEQ%Lw4(1:N)=DD%W4(1:N)/sigLC/2.d0 ! sqrt(2.d0)
!        end if
        
        ! Initial precision and random walk of ionosphere parameter
        do i=2*IonoNum+ParaNum+1, 3*IonoNum+ParaNum
!            if (ADmethod=='LS') then
!                if (Epo_NEQ%Nbb(i,i)==0.d0 .and. any(Epo_NEQ%Al1(:,i)/=0.d0)) then
!                    Epo_NEQ%Nbb(i,i)=Epo_NEQ%Nbb(i,i)+(1.d0/0.2d0**2)    !    (1/0.2**2) ! 0.2m
!                end if
!            elseif (ADmethod=='KF') then ! As ionosphere parameter is estimated as randon walk, so we use transformed Kalman Filter instead of Least Square
                if (Epo_NEQ%InvN(i,i)==0.d0 .and. any(Epo_NEQ%Al1(:,i)/=0.d0)) then
!                    Epo_NEQ%InvN(i,i)=Epo_NEQ%InvN(i,i)+0.01d0     !    (0.1**2) ! 0.1m
                elseif (Epo_NEQ%InvN(i,i)>0.d0) then
                    ! Random walk of ionosphere delay
                    Epo_NEQ%InvN(i,i) = Epo_NEQ%InvN(i,i)+0.5**2/3600.d0*Interval   ! (4m2/3600*Interval)
                end if
!            end if
        end do
        ! Initial precision and random walk of troposphere parameter
!        if (TropLen/=0.d0 .and. ADmethod=='LS') then
!            Epo_NEQ%Nbb(4,4)=Epo_NEQ%Nbb(4,4)+(1.d0/0.02d0**2)    !    (1/0.02**2) ! 0.02m
        if (TropLen/=0.d0) then ! As ionosphere parameter is estimated as randon walk, so we use transformed Kalman Filter instead of Least Square
            if (Epo_NEQ%InvN(4,4)==0.d0) then
!                Epo_NEQ%InvN(4,4)=Epo_NEQ%InvN(4,4)+0.01d0**2   !    (0.01**2) ! 0.01m
            else
                Epo_NEQ%InvN(4,4)=Epo_NEQ%InvN(4,4)+(0.01**2/3600.d0*Interval)    !    (0.01^2/3600*Interval) ! 0.01cm
            end if
        end if

        ! *****************************************************
        if (ADmethod=='LS') then  
        ! This is transformed Kalman Filter, similar as Least Square.
        ! For long baseline, we use transformed KF instead of LS
        ! Traditional Kalman Filter formulation is more efficient when parameter is very huge
            call InvSqrt(Epo_NEQ%InvN, Epo_NEQ%N, Epo_NEQ%Nbb)     ! Inverse is a bit slow
            Epo_NEQ%U=MATMUL(Epo_NEQ%Nbb, Epo_NEQ%dx)   ! Update Nbb and U
        end if
        ! Now we can start traditional Least Square
        ! ****************************************************

        ! Add constraints to ionosphere
        factor=(Baseline*6.d-6)*exp((90.d0-Min_Lat)/50.d0-1.d0)  ! distance related and latitude related
        if (Baseline>200.d3) then
            factor=1.2d0*exp((90.d0-Min_Lat)/50.d0-1.d0)
        end if
        do i=2*IonoNum+ParaNum+1, 3*IonoNum+ParaNum
            if (any(Epo_NEQ%Al1(:,i)/=0.d0)) then
                    if (ADmethod=='LS') then
                        Epo_NEQ%Nbb(i,i)=Epo_NEQ%Nbb(i,i)+(1.d0/factor**2) ! 0.4m
                    elseif (ADmethod=='KF') then
                        call KF_Gain_one(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, i, 0.d0, factor)
                    end if
            end if
        end do
        ! Add constraints to troposphere
        factor=log(1.d0+Baseline/5.d3)*0.01d0+Diff_Hgt*1.d-5  ! distance related and height related
        if (TropLen/=0.d0) then
            if (ADmethod=='LS') then
                Epo_NEQ%Nbb(4,4)=Epo_NEQ%Nbb(4,4)+(1.d0/factor**2) ! 0.05m
            elseif (ADmethod=='KF') then
                call KF_Gain_one(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, 4, 0.d0, factor)
            end if
        end if
        
        ! Fix and Hold mode, add constraints to fixed ambiguity, See RTKLIB rtkpos.c  holdamb()
        if (ar_mode==3) then !
            do PRN=1, 2*maxPRN
                if (If_Est_Iono .and. IonoNum>0 .and. Epo_NEQ%fixed_amb(PRN)/=0.99d0) then
                    Epo_NEQ%fixed_amb_num(PRN)=Epo_NEQ%fixed_amb_num(PRN)+1
                    if (PRN<=maxPRN) then
                        if (Epo_NEQ%fixed_amb_num(PRN)>50 .and. Epo_NEQ%Ele(PRN)>HoldEle) then  ! If the same ambiguity > 50 epoches
                            if (ADmethod=='LS') then
                                Epo_NEQ%Nbb(PRN+ParaNum,PRN+ParaNum)=Epo_NEQ%Nbb(PRN+ParaNum,PRN+ParaNum)+(1.d0/0.01d0**2) !  ! 1cm
                                Epo_NEQ%U(PRN+ParaNum)=Epo_NEQ%U(PRN+ParaNum)+1.d0*(Epo_NEQ%fixed_amb(PRN))*(1.d0/0.01d0**2)
                            elseif (ADmethod=='KF') then
                                call KF_Gain_one(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, PRN+ParaNum, 1.d0*(Epo_NEQ%fixed_amb(PRN)), 0.01d0)
                            end if
                        end if
                    elseif (PRN>maxPRN) then
                        if (Epo_NEQ%fixed_amb_num(PRN)>50 .and. Epo_NEQ%Ele(PRN-maxPRN)>HoldEle) then  ! If the same ambiguity >50 epoches
                            if (ADmethod=='LS') then
                                Epo_NEQ%Nbb(PRN+ParaNum,PRN+ParaNum)=Epo_NEQ%Nbb(PRN+ParaNum,PRN+ParaNum)+(1.d0/0.01d0**2) ! (1/0.01**2) ! 1cm
                                Epo_NEQ%U(PRN+ParaNum)=Epo_NEQ%U(PRN+ParaNum)+real(Epo_NEQ%fixed_amb(PRN))*(1.d0/0.01d0**2)
                             elseif (ADmethod=='KF') then
                                call KF_Gain_one(Epo_NEQ%InvN, Epo_NEQ%dx,Epo_NEQ%N, PRN+ParaNum, 1.d0*(Epo_NEQ%fixed_amb(PRN)), 0.01d0)
                            end if
                        end if
                    end if
                end if
            end do
        end if
        
        if (ADmethod=='LS') then     ! When Least square, add Nbb
            Epo_NEQ%Nbb     =  Epo_NEQ%Nbb + matmul(  matmul( transpose(Epo_NEQ%Ap1(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Ap1(1:N,:)  )
            Epo_NEQ%Nbb     =  Epo_NEQ%Nbb +  matmul(  matmul( transpose(Epo_NEQ%Ap2(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Ap2(1:N,:)  )
            Epo_NEQ%Nbb     =  Epo_NEQ%Nbb +  matmul(  matmul( transpose(Epo_NEQ%Al1(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Al1(1:N,:)  )
            Epo_NEQ%Nbb     =  Epo_NEQ%Nbb +  matmul(  matmul( transpose(Epo_NEQ%Al2(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Al2(1:N,:)  )

            Epo_NEQ%U         =  Epo_NEQ%U           +  matmul(  matmul( transpose(Epo_NEQ%Ap1(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Lp1(1:N)  )
            Epo_NEQ%U         =  Epo_NEQ%U           +  matmul(  matmul( transpose(Epo_NEQ%Ap2(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Lp2(1:N)  )
            Epo_NEQ%U         =  Epo_NEQ%U           +  matmul(  matmul( transpose(Epo_NEQ%Al1(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Ll1(1:N)  )
            Epo_NEQ%U         =  Epo_NEQ%U           +  matmul(  matmul( transpose(Epo_NEQ%Al2(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Ll2(1:N)  )
        elseif (ADmethod=='KF') then
            ! Kalmen Filter will done in Solve_NEQ_Iono, KF_Gain function
        end if
    else 
        ! Traditional Least Square
        Epo_NEQ%Nbb(1:ParaNum,1:ParaNum)     =  Epo_NEQ%Nbb(1:ParaNum,1:ParaNum) + matmul(  matmul( transpose(Epo_NEQ%Ap1(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Ap1(1:N,:)  )
        Epo_NEQ%Nbb(1:ParaNum,1:ParaNum)     =  Epo_NEQ%Nbb(1:ParaNum,1:ParaNum) +  matmul(  matmul( transpose(Epo_NEQ%Ap2(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Ap2(1:N,:)  )    
        Epo_NEQ%U(1:ParaNum)         =  Epo_NEQ%U(1:ParaNum)           +  matmul(  matmul( transpose(Epo_NEQ%Ap1(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Lp1(1:N)  )
        Epo_NEQ%U(1:ParaNum)         =  Epo_NEQ%U(1:ParaNum)           +  matmul(  matmul( transpose(Epo_NEQ%Ap2(1:N,:)), Epo_NEQ%P(1:N,1:N) ), Epo_NEQ%Lp2(1:N)  )    
    end if

    return
end subroutine