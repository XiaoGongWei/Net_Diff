! ========================================
!                         CS_Detect
!  PURPOSE
!     Cyccle-slip detect in station difference
!
! Inputs:
!           SD                         Station difference strcuture
!           PRN                      PRN
!           Ele                        Satellite elevation
!
! Output:
!            Slip                       Cycle slip flag
!                                         ¡¯0¡® for no, '1' for cycle slip
! 
! Written by : Yize Zhang
! =======================================

subroutine  CS_Detect(SD, PRN, Ele, Slip) 
use MOD_PreSD
use MOD_SD
use MOD_constant
use MOD_VAR
use MOD_FileID
implicit none
    type(type_SD)  :: SD
    integer             ::  PRN
    real(8)              :: Ele
    integer             ::  Slip
    ! Local variables
    integer             ::  i, j
    real(8)              ::  L1, L2, P1, P2
    real(8)              :: GF, MW
    real(8)              :: dT, dGF, dMW
    real(8)              :: TT(3), GFEst, sigma2, meanGFPrev, meanTT, coe, sigma
    real(8)              :: GFThreshold=0.d0, MWThreshold=0.d0
    real(8)              :: P, L, L1P1, dL1P1, L1P1Threshold=0.d0
    
    Slip=0

    do i=1, SD%PRNS
        if (SD%PRN(i)==PRN) exit
    end do
    L1=SD%L1(i)      ! In distance(meter)
    L2=SD%L2(i)

!    ! ======== GF combination ======
!!    if ( (index(CSmethod,"GF")/=0) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
!        if ((L1==0.d0) .or. (L2==0.d0) ) then
!            return  ! need to change if single frequency data
!        end if
!        GF=L1*f1/c-L2*f1/c         ! In cycle
!        dT=(SD%Week-PreSD(PRN)%WeekGF)*604800.d0+(SD%sow-PreSD(PRN)%sowGF)
!        if (dT>300.d0) then
!            Slip=1
!        elseif  (PreSD(PRN)%GF/=0.d0)  then
!            dGF=GF - PreSD(PRN)%GF
!            dGF=dGF*dsind(Ele)
!            if ( abs(dGF)>=0.4d0 ) Slip=1
!        end if
!
!        PreSD(PRN)%WeekGF=SD%Week
!        PreSD(PRN)%sowGF=SD%sow
!        PreSD(PRN)%GF=GF      ! In cycle
!!    end if
!
!    ! =====M-W combination======
!    P1=SD%P1CS(i)
!    P2=SD%P2CS(i)
!!    if ( (index(CSmethod,"MW")/=0) .and. (P1/=0.d0) .and. (P2/=0.d0) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
!    if ((P1 ==0.0d0) .or. (P2 ==0.0d0)) return  ! If no P1 or P2 data of this satellite
!        MW=(f1*L1-f2*L2)/(f1-f2) - (f1*P1 +f2*P2)/(f1+f2)  ! In distance(meter)
!        MW=MW/c*(f1-f2)      ! In cycle 
!        if  (PreSD(PRN)%MW/=0.d0)  then
!            dMW=MW - PreSD(PRN)%MW
!            if ( abs(dMW)>=3.d0 ) Slip=1
!        end if
!
!        PreSD(PRN)%WeekMW=SD%Week
!        PreSD(PRN)%sowMW=SD%sow
!        PreSD(PRN)%MW=MW      ! In cycle
!!    end if

     ! ======== GF combination ======
    if ( (index(CSmethod,"GF")/=0) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
        GF=L1*f1/c-L2*f1/c         ! In cycle
        TT=(PreSD(PRN)%WeekPrev-SD%Week)*604800.0d0+PreSD(PRN)%SowPrev-SD%Sow
        if ( TT(3)<-Interval*10.d0) then
            Slip=1  ! If data missing too long or new satellite ascending
        elseif (PreSD(PRN)%arcLengthGF<3)then
            Slip=2  ! If not enough epoches,accumulate data
        else
            ! Quadratic polyfit
            !call Lagrange(TT,PreSD(PRN)%GFPrev,0.d0,GFEst,2) 
            GFEst=PreSD(PRN)%GFPrev(1)*TT(2)*TT(3)/(TT(1)-TT(2))/(TT(1)-TT(3))+ &
                        PreSD(PRN)%GFPrev(2)*TT(1)*TT(3)/(TT(2)-TT(1))/(TT(2)-TT(3))+ &
                        PreSD(PRN)%GFPrev(3)*TT(1)*TT(2)/(TT(3)-TT(1))/(TT(3)-TT(2))
            ! Linear polyfit
            meanGFPrev=sum(PreSD(PRN)%GFPrev(1:3))/3
            meanTT=sum(TT(1:3))/3
            coe=DOT_PRODUCT(TT(1:3)-meanTT,PreSD(PRN)%GFPrev(1:3))/ &
                DOT_PRODUCT(TT(1:3)-meanTT,TT(1:3)-meanTT)
            GFEst=meanGFPrev-meanTT*coe  ! a0
            sigma=sqrt(DOT_PRODUCT(PreSD(PRN)%GFPrev(1:3)-GFEst-a1*TT(1:3),PreSD(PRN)%GFPrev(1:3)-GFEst-a1*TT(1:3)))
            GFThreshold=csGFmax+(csGFmin-csGFmax)*exp(TT(3)/csGFdt)
            if (abs(GF-GFEst)*c/f1>GFThreshold*1.d0) Slip=1  ! in diatance
        end if
    end if
    
    ! ======== M-W combination ======
    ! Lw=(f1*L1-f2*L2)/(f1-f2), Pw=(f1*P1 +f2*P2)/(f1+f2)
    !  Warn: Above L1 and L2 is in meter
    !  rantaMW*(N1-N2)=Lw-Pw=(f1*L1-f2*L2)/(f1-f2) - (f1*P1 +f2*P2)/(f1+f2), rantaMW=c/(f1-f2)
    !  MW=N1-N2
    !  dMW=MW2-MW1
    
    P1=SD%P1CS(i)
    P2=SD%P2CS(i)
    if ( (index(CSmethod,"MW")/=0) .and. (P1/=0.d0) .and. (P2/=0.d0) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
!    if ((P1 ==0.0d0) .or. (P2 ==0.0d0)) return  ! If no P1 or P2 data of this satellite, impossible happen
    MW=(f1*L1-f2*L2)/(f1-f2) - (f1*P1 +f2*P2)/(f1+f2)  ! In distance(meter)
    MW=MW/c*(f1-f2)                                              ! In cycle 
!    MW=L1 - L2 - (P1*f1 + P2*f2)/(f1+f2)*(f1-f2)/c ! L1, L2 here is in cycle
    if  (PreSD(PRN)%arcLengthMW>=3)  then
        dMW=MW-PreSD(PRN)%nMWmean
        sigma2=PreSD(PRN)%nMWmean2 - PreSD(PRN)%nMWmean**2
        MWThreshold=min(csMWmax, max(csMWslope*sqrt(sigma2),csMWmin))
        if (abs(dMW)>=MWThreshold) Slip=1
    end if 
    end if

    ! initial cycle slip information if cycle slip happens
    if ( (Slip==1) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
        PreSD(PRN)%WeekPrev=0.d0
        PreSD(PRN)%SowPrev=0.d0
        PreSD(PRN)%GFPrev=0.d0
        PreSD(PRN)%arcLengthGF=0
        PreSD(PRN)%arcLengthMW=0
        write(CSID,"(I3,F6.2,I6,5F10.4,I3)")  PRN, ele, int(-TT(3)), abs(GF-GFEst)*c/f1, GFThreshold,3.d0*sigma*c/f1, dMW, MWThreshold, Slip
    elseif (Slip==2) then
       Slip=1
    end if
        
    ! Update GF information
    if ( (index(CSmethod,"GF")/=0) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
        PreSD(PRN)%arcLengthGF=PreSD(PRN)%arcLengthGF+1
        if (PreSD(PRN)%arcLengthGF>3) PreSD(PRN)%arcLengthGF=3
        PreSD(PRN)%WeekPrev=eoshift(PreSD(PRN)%WeekPrev,shift=1,boundary=SD%Week)
        PreSD(PRN)%SowPrev=eoshift(PreSD(PRN)%SowPrev,shift=1,boundary=SD%Sow)
        PreSD(PRN)%GFPrev=eoshift(PreSD(PRN)%GFPrev,shift=1,boundary=GF)
    end if

    ! Update MW information
    if ( (index(CSmethod,"MW")/=0) .and. (P1/=0.d0) .and. (P2/=0.d0) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
    PreSD(PRN)%arcLengthMW=PreSD(PRN)%arcLengthMW+1
    if (PreSD(PRN)%arcLengthMW>10) PreSD(PRN)%arcLengthMW=10
    PreSD(PRN)%nMWmean=PreSD(PRN)%nMWmean+(MW-PreSD(PRN)%nMWmean)/PreSD(PRN)%arcLengthMW
    PreSD(PRN)%nMWmean2=PreSD(PRN)%nMWmean2+(MW**2-PreSD(PRN)%nMWmean2)/PreSD(PRN)%arcLengthMW
    end if

    ! ========= LLI ============
    ! Loss of Lock Indicator (LLI) 
    ! From RINEX 3.02 and gLAB preprocessing.c
	!		//   Range: 0-7
	!		//      0 or blank: OK or not known
 	!		//      1, 3, 5 and 7: Lost lock between previous and current observation: cycle slip possible
    if (index(CSmethod,"LLI")/=0) then
        if (index(ObsCombine,"PC")/=0) then
            if ( (mod(PreSD(PRN)%LLI1,2)==1) .or. (mod(PreSD(PRN)%LLI2,2)==1) ) then
    !        if ( (LLI1>0) .or. (LLI2>0) ) then
                Slip=1
                write(CSID,"(I3,F6.2,2I5,I3)")  PRN, ele, PreSD(PRN)%LLI1, PreSD(PRN)%LLI2, Slip
            end if
        elseif ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0)) then
            if (mod(PreSD(PRN)%LLI1,2)==1)  then
                Slip=1
                write(CSID,"(I3,F6.2,2I5,I3)")  PRN, ele, PreSD(PRN)%LLI1, PreSD(PRN)%LLI2, Slip
            end if
        elseif ( ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) .and. (freq_comb=='L2L3')) then
            if (mod(PreSD(PRN)%LLI1,2)==1)  then
                Slip=1
                write(CSID,"(I3,F6.2,2I5,I3)")  PRN, ele, PreSD(PRN)%LLI1, PreSD(PRN)%LLI2, Slip
            end if
        else
            if (mod(PreSD(PRN)%LLI2,2)==1)  then
                Slip=1
                write(CSID,"(I3,F6.2,2I5,I3)") PRN, ele, PreSD(PRN)%LLI1, PreSD(PRN)%LLI2, Slip
            end if        
        end if
    end if

        
    ! ========= LLI ============
    ! Loss of Lock Indicator (LLI) 
    ! From RINEX 3.02 and gLAB preprocessing.c
	!		//   Range: 0-7
	!		//      0 or blank: OK or not known
 	!		//      1, 3, 5 and 7: Lost lock between previous and current observation: cycle slip possible
!    if (index(CSmethod,"LLI")/=0) then
!        if ( (mod(LLI1,2)==1) .or. (mod(LLI2,2)==1) ) then
!            Slip=1
!            write(CSID,"(I3,F6.2,2I5,I3)")  PRN, ele, LLI1, LLI2, Slip
!        end if
!    end if
    
end subroutine