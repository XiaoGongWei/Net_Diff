! ================= Reset_Amb ===================== 
! PURPOSE:
!            This subroutine is to reset ambiguity if reference satellite
!   change in tight combination.
!

! INPUTS:
!         SD                             zero difference structure
!         DD                            double difference structure
!         RefSat                       reference satellite
!         NEQ                         normal equation structure
!         Epo_NEQ                 current epoch normal equation structure
! OUTPUT:
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! =================End of header=================

subroutine Reset_Amb(ZD, DD, RefSat, NEQ, Epo_NEQ)
use MOD_ZD
use MOD_DD
use MOD_NEQ
use MOD_Epo_NEQ
use MOD_GLO_Fre
use MOD_Var
implicit none
    type(type_ZD) :: ZD(2)
    type(type_DD) :: DD
    integer :: RefSat
    type(type_NEQ) :: NEQ
    type(type_Epo_NEQ) :: Epo_NEQ
    ! lLocal variables
    integer :: i, PRN, sys, freq
    real(8) :: damb_rov(3), damb_ref(3), damb(3)
    real(8) :: ref_f1, ref_f2, ref_f3

    do i=1, DD%PRNS
        PRN=DD%PRN(i)
        sys=DD%Sys(i)
        damb_rov=ZD(2)%amb1(PRN,1:3)-ZD(1)%amb1(PRN,1:3)-ZD(2)%amb0(PRN,1:3)+ZD(1)%amb0(PRN,1:3)
        damb_ref=ZD(2)%amb1(RefSat,1:3)-ZD(1)%amb1(RefSat,1:3)-ZD(2)%amb0(RefSat,1:3)+ZD(1)%amb0(RefSat,1:3)

        damb=damb_rov-damb_ref   ! In cycle

!                if (any(damb(1:2)/=0.d0)) then  ! If DD intial ambiguity change, then change it
        NEQ%dx(ParaNum+PRN)=NEQ%dx(ParaNum+PRN)-damb(1)
        NEQ%dx(ParaNum+PRN+SatNum)=NEQ%dx(ParaNum+PRN+SatNum)-damb(2)
        if (If_Est_Iono .and. IonoNum>0) then
            Epo_NEQ%dx(ParaNum+PRN)=Epo_NEQ%dx(ParaNum+PRN)-damb(1)
            Epo_NEQ%dx(ParaNum+PRN+SatNum)=Epo_NEQ%dx(ParaNum+PRN+SatNum)-damb(2)
        end if
        ZD(1)%amb0(PRN,1:3)=ZD(1)%amb1(PRN,1:3)     ! Reset ZD ambiguity, in cycle
        ZD(2)%amb0(PRN,1:3)=ZD(2)%amb1(PRN,1:3)
            
        if (sys==1) then   ! GPS/QZSS
            f1=f_L1
            f2=f_L2
            f3=f_L5
        elseif (sys==2) then   ! GLONASS
            freq=Fre_Chann(PRN-GNum)
            f1=(1602.0d0+freq*0.5625d0)*1.0D6   ! f1=(1602.0d0+K*0.5625d0)*1.0d6
            f2=(1246.0d0+freq*0.4375d0)*1.0D6
        elseif  (sys==3) then   ! COMPASS
            if (freq_comb=='L1L2') then
                f1=f_B1
                f2=f_B2
                f3=f_B3
            elseif (freq_comb=='L1L3') then
                f1=f_B1
                f2=f_B3
            elseif (freq_comb=='L2L3') then
                f1=f_B2
                f2=f_B3
            end if
        elseif  (sys==4) then ! GALILEO
            if (freq_comb=='L1L2') then   ! E1 E5a
                f1=f_E1
                f2=f_E5
                f3=f_E5b
            elseif (freq_comb=='L1L3') then   ! E1 E5b
                f1=f_E1
                f2=f_E5b
            elseif (freq_comb=='L2L3') then   ! E5a E5b
                f1=f_E5a
                f2=f_E5b
            end if
        elseif (sys==5) then   ! IRNSS
            f1=f_L1
            f2=f_S
        end if
        if (RefSat<=GNum) then
            ref_f1=f_L1
            ref_f2=f_L2
            ref_f3=f_L5
        elseif (RefSat<=GNum+RNum) then
            freq=Fre_Chann(PRN-GNum)
            ref_f1=(1602.0d0+freq*0.5625d0)*1.0D6   ! f1=(1602.0d0+K*0.5625d0)*1.0d6
            ref_f2=(1246.0d0+freq*0.4375d0)*1.0D6
        elseif  (RefSat<=GNum+RNum+CNum) then   ! COMPASS
            if (freq_comb=='L1L2') then
                ref_f1=f_B1
                ref_f2=f_B2
                ref_f3=f_B3
            elseif (freq_comb=='L1L3') then
                ref_f1=f_B1
                ref_f2=f_B3
            elseif (freq_comb=='L2L3') then
                ref_f1=f_B2
                ref_f2=f_B3
            end if
        elseif  (RefSat<=GNum+RNum+CNum+NumE) then ! GALILEO
            if (freq_comb=='L1L2') then   ! E1 E5a
                ref_f1=f_E1
                ref_f2=f_E5
                ref_f3=f_E5b
            elseif (freq_comb=='L1L3') then   ! E1 E5b
                ref_f1=f_E1
                ref_f2=f_E5b
            elseif (freq_comb=='L2L3') then   ! E5a E5b
                ref_f1=f_E5a
                ref_f2=f_E5b
            end if
        elseif (RefSat<=GNum+RNum+CNum+NumE+JNum) then  ! QZSS
            ref_f1=f_L1
            ref_f2=f_L2
            ref_f3=f_L5
        elseif (RefSat<=GNum+RNum+CNum+NumE+JNum+INum) then   ! IRNSS
            ref_f1=f_L1
            ref_f2=f_S
        end if
        DD%L1(i)=DD%L1(i)-damb_rov(1)*c/f1 + damb_ref(1)*c/ref_f1   ! Reset DD ambiguity, in meter
        DD%L2(i)=DD%L2(i)-damb_rov(2)*c/f2 + damb_ref(2)*c/ref_f2
        DD%WL(i)=DD%WL(i)-damb_rov(1)*c/f1 + damb_ref(1)*c/ref_f1
        DD%W4(i)=DD%W4(i)-damb_rov(2)*c/f2 + damb_ref(2)*c/ref_f2

!        if (PRN==DD%RefSat(1)) then  ! reset DISB due to frequency difference
!            NEQ%dx(4+sys*4-1)=NEQ%dx(4+sys*4-1) + damb_rov(1)*(c/f_L1-c/f_B1)  ! old reference ambiguity cycle change
!            NEQ%dx(4+sys*4)=NEQ%dx(4+sys*4) + damb_rov(2)*(c/f_L2-c/f_B2)
!        end if
!                end if
    end do
    ZD(1)%amb0(RefSat,1:3)=ZD(1)%amb1(RefSat,1:3)     ! Reset intial ambiguity of reference satellite
    ZD(2)%amb0(RefSat,1:3)=ZD(2)%amb1(RefSat,1:3)

    return
end subroutine