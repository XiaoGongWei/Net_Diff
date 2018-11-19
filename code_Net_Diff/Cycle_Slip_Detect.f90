! ==========Cycle Slip Detection=======
! Inputs:
!       ObsData    Obsdata
!       epoch        Judge that if it is epoch 1
!       j                 The order of the PRN in ObsData
!       Ele              Elevation of the satellite
!
! Output:
!       Slip             Cycle Slip Flag, 1 means that a cycle slip happens
!
! Written by: Yize Zhang
! =============End of Header============

subroutine Cycle_Slip_Detect(sta, week, sow, P1, P2, L1, L2, LLI1, LLI2, PRN, Ele, Slip) 
use MOD_constant
use MOD_CycleSlip
use MOD_VAR
use MOD_FileID
implicit none
    integer            ::  sta, week
    real(8)              :: sow 
    integer(1)        :: LLI1, LLI2
    real(8)              ::  L1, L2, P1, P2
    integer             ::  PRN
    real(8)              :: Ele
    ! Intent out
    integer(1) :: Slip
    
    real(8)              :: GF, MW
    real(8)              :: dT, dGF, dMW
    real(8)              :: TT(3), GFEst, sigma2, meanGFPrev, meanTT, coe, sigma
    real(8)              :: GFThreshold=0.d0, MWThreshold=0.d0
    real(8)              :: P, L, L1P1, dL1P1, L1P1Threshold=0.d0
    
!    Slip=0
!    CycleSlip(sta)%CScount=99  ! do not use satellite difference
    if (CycleSlip(sta)%CScount==99) Slip=0
        
!    CycleSlip(sta)%GF(1,PRN)=CycleSlip(sta)%GF(2,PRN)
!    CycleSlip(sta)%GF(2,PRN)=9999.d0
!    CycleSlip(sta)%MW(1,PRN)=CycleSlip(sta)%MW(2,PRN)
!    CycleSlip(sta)%MW(2,PRN)=9999.d0
!        
!    ! L4 combination, ionosphere only combination, unit in cycle
!    ! GF=(L1-L2)m=ranta1*Phase1 - ranta2*Phase2
!    ! dL=GF2-GF1
!    ! According to Han(Han, Quality Control Issuses Relating to Ambiguity Resolution for Real-Time  GPS Kinematic Positioning[J]. Journal of Geodesy,1997,71(6):351-361
!    ! If dL<5.0cm, it means that a cycle slip happens. It equals 0.205 dGF here.
!    CycleSlip(sta)%GF(2,PRN)=L1-L2*f1/f2    ! L1, L2 unit in cycle
!    if ( (CycleSlip(sta)%dT(PRN)>=180.d0) .and. (CycleSlip(sta)%CScount==99) ) then
!        Slip=1
!    elseif (dabs(CycleSlip(sta)%GF(1,PRN)- 9999.0d0)>0.1d0) then
!        dGF=CycleSlip(sta)%GF(2,PRN)-CycleSlip(sta)%GF(1,PRN)
!        dGF=dGF/dsqrt(CycleSlip(sta)%dT(PRN))*dsind(Ele)   ! dsqrt(CycleSlip(sta)%dT(PRN))
!        if ( (abs(dGF)>=0.07d0) .and. (CycleSlip(sta)%CScount==99) ) then   ! 0.8 cycle
!            Slip=1
!        end if
!    else
!        Slip=1
!    end if
!        
!    CycleSlip(sta)%dT(PRN)=0.d0
!        
!    if ((P1 ==0.0d0) .or. (P2 ==0.0d0)) return  ! If no P1 or P2 data of this satellite
!        
!    !! =====M-W combination, unit in cycle======
!    ! Lw=(f1*L1-f2*L2)/(f1-f2), Pw=(f1*P1 +f2*P2)/(f1+f2)
!    !  Warn: Above L1 and L2 is in meter
!    !  rantaMW*(N1-N2)=Lw-Pw=(f1*L1-f2*L2)/(f1-f2) - (f1*P1 +f2*P2)/(f1+f2), rantaMW=c/(f1-f2)
!    !  MW=N1-N2
!    !  dMW=MW2-MW1
!    CycleSlip(sta)%MW(2,PRN)=L1 - L2 - (P1*f1 + P2*f2)/(f1+f2)*(f1-f2)/c ! L1, L2 here is in cycle
!    if (dabs(CycleSlip(sta)%MW(1,PRN)-9999.0d0)>0.1d0) then
!        dMW=CycleSlip(sta)%MW(2,PRN)-CycleSlip(sta)%MW(1,PRN)
!        if ( (abs(dMW)>=3.d0) .and. (CycleSlip(sta)%CScount==99) ) then  ! 5 cycle
!            Slip=1
!        end if
!    end if
!
!    if (Slip==1)then
!        write(CSID,"(A5,I5,3F8.2,I3)")  "CS",PRN, ele,dGF, dMW, Slip
!    end if
    
    ! ======== GF combination ======
    if ( (index(CSmethod,"GF")/=0) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
        GF=L1-L2*f1/f2        ! In cycle
        TT=(CycleSlip(sta)%CS(PRN)%WeekPrev-Week)*604800.0d0+CycleSlip(sta)%CS(PRN)%SowPrev-Sow
        if ( (TT(3)<-Interval*10.d0) .and. (CycleSlip(sta)%CScount==99) ) then
            Slip=1  ! If data missing too long or new satellite ascending
        elseif ( (CycleSlip(sta)%CS(PRN)%arcLengthGF<3) .and. (CycleSlip(sta)%CScount==99) )then
!            Slip=2  ! If not enough epoches,accumulate data
            dGF=GF-CycleSlip(sta)%CS(PRN)%GFPrev(3)
!            dGF=dGF*dsind(Ele)   ! dsqrt(CycleSlip(sta)%dT(PRN))
            if ( abs(dGF)*c/f1>=csGFmin )  Slip=1   ! Be strict for the first 3 epoches
        else
            ! Quadratic polyfit
            !call Lagrange(TT,CycleSlip(sta)%CS(PRN)%GFPrev,0.d0,GFEst,2) 
            GFEst=CycleSlip(sta)%CS(PRN)%GFPrev(1)*TT(2)*TT(3)/(TT(1)-TT(2))/(TT(1)-TT(3))+ &
                        CycleSlip(sta)%CS(PRN)%GFPrev(2)*TT(1)*TT(3)/(TT(2)-TT(1))/(TT(2)-TT(3))+ &
                        CycleSlip(sta)%CS(PRN)%GFPrev(3)*TT(1)*TT(2)/(TT(3)-TT(1))/(TT(3)-TT(2))
            ! Linear polyfit
            meanGFPrev=sum(CycleSlip(sta)%CS(PRN)%GFPrev(1:3))/3.d0
            meanTT=sum(TT(1:3))/3.d0
            coe=DOT_PRODUCT(TT(1:3)-meanTT,CycleSlip(sta)%CS(PRN)%GFPrev(1:3))/ &
                DOT_PRODUCT(TT(1:3)-meanTT,TT(1:3)-meanTT)
            GFEst=meanGFPrev-meanTT*coe  ! a0
            sigma=sqrt(DOT_PRODUCT(CycleSlip(sta)%CS(PRN)%GFPrev(1:3)-GFEst-a1*TT(1:3),CycleSlip(sta)%CS(PRN)%GFPrev(1:3)-GFEst-a1*TT(1:3)))
            GFThreshold=(csGFmax+(csGFmin-csGFmax)*exp(TT(3)/csGFdt))*sigLC/0.01d0  ! Default sigLC is 0.01m
            if (abs(GF-GFEst)*c/f1>GFThreshold*0.8d0 .and. CycleSlip(sta)%CScount==99 ) Slip=1  ! in diatance, in case of 1 cycle jump in L1 and L2 in RTK
        end if
    end if
    
    ! ======== M-W combination ======
    ! Lw=(f1*L1-f2*L2)/(f1-f2), Pw=(f1*P1 +f2*P2)/(f1+f2)
    !  Warn: Above L1 and L2 is in meter
    !  rantaMW*(N1-N2)=Lw-Pw=(f1*L1-f2*L2)/(f1-f2) - (f1*P1 +f2*P2)/(f1+f2), rantaMW=c/(f1-f2)
    !  MW=N1-N2
    !  dMW=MW2-MW1
    if ( (index(CSmethod,"MW")/=0) .and. (P1/=0.d0) .and. (P2/=0.d0) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
    !    if ((P1 ==0.0d0) .or. (P2 ==0.0d0)) return  ! If no P1 or P2 data of this satellite, impossible happen
    !    MW=(f1*L1-f2*L2)/(f1-f2) - (f1*P1 +f2*P2)/(f1+f2)  ! In distance(meter)
    !    MW=MW/c*(f1-f2)                                              ! In cycle 
        MW=L1 - L2 - (P1*f1 + P2*f2)/(f1+f2)*(f1-f2)/c ! L1, L2 here is in cycle
        if  (CycleSlip(sta)%CS(PRN)%arcLengthMW>=3)  then
            dMW=MW-CycleSlip(sta)%CS(PRN)%nMWmean
            sigma2=CycleSlip(sta)%CS(PRN)%nMWmean2 - CycleSlip(sta)%CS(PRN)%nMWmean**2
            MWThreshold=min(csMWmax, max(csMWslope*sqrt(sigma2),csMWmin))
            MWThreshold=MWThreshold*sigLC/0.01d0  ! Default sigLC is 0.01m
            if ( (abs(dMW)>=MWThreshold)  .and. (CycleSlip(sta)%CScount==99) ) Slip=1
        end if 
    end if

    ! initial cycle slip information if cycle slip happens
    if ( (Slip==1) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
        CycleSlip(sta)%CS(PRN)%WeekPrev=0.d0
        CycleSlip(sta)%CS(PRN)%SowPrev=0.d0
        CycleSlip(sta)%CS(PRN)%GFPrev=0.d0
        CycleSlip(sta)%CS(PRN)%arcLengthGF=0
        CycleSlip(sta)%CS(PRN)%arcLengthMW=0
        write(CSID,"(2I3,F6.2,I6,5F10.4,I3)")  sta, PRN, ele, int(-TT(3)), abs(GF-GFEst)*c/f1, GFThreshold,3.d0*sigma*c/f1, dMW, MWThreshold, Slip
    elseif (Slip==2) then
       Slip=1
    end if
        
    ! Update GF information
    if ( (index(CSmethod,"GF")/=0) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
        CycleSlip(sta)%CS(PRN)%arcLengthGF=CycleSlip(sta)%CS(PRN)%arcLengthGF+1
        if (CycleSlip(sta)%CS(PRN)%arcLengthGF>3) CycleSlip(sta)%CS(PRN)%arcLengthGF=3
        CycleSlip(sta)%CS(PRN)%WeekPrev=eoshift(CycleSlip(sta)%CS(PRN)%WeekPrev,shift=1,boundary=Week)
        CycleSlip(sta)%CS(PRN)%SowPrev=eoshift(CycleSlip(sta)%CS(PRN)%SowPrev,shift=1,boundary=Sow)
        CycleSlip(sta)%CS(PRN)%GFPrev=eoshift(CycleSlip(sta)%CS(PRN)%GFPrev,shift=1,boundary=GF)
    end if

    ! Update MW information
    if ( (index(CSmethod,"MW")/=0) .and. (P1/=0.d0) .and. (P2/=0.d0) .and. (L1/=0.d0) .and. (L2/=0.d0) ) then
        CycleSlip(sta)%CS(PRN)%arcLengthMW=CycleSlip(sta)%CS(PRN)%arcLengthMW+1
        if (CycleSlip(sta)%CS(PRN)%arcLengthMW>120) CycleSlip(sta)%CS(PRN)%arcLengthMW=120
        CycleSlip(sta)%CS(PRN)%nMWmean=CycleSlip(sta)%CS(PRN)%nMWmean+(MW-CycleSlip(sta)%CS(PRN)%nMWmean)/CycleSlip(sta)%CS(PRN)%arcLengthMW
        CycleSlip(sta)%CS(PRN)%nMWmean2=CycleSlip(sta)%CS(PRN)%nMWmean2+(MW**2-CycleSlip(sta)%CS(PRN)%nMWmean2)/CycleSlip(sta)%CS(PRN)%arcLengthMW
    end if

    ! ========= LLI ============
    ! Loss of Lock Indicator (LLI) 
    ! From RINEX 3.02 and gLAB preprocessing.c
	!		//   Range: 0-7
	!		//      0 or blank: OK or not known
 	!		//      1, 3, 5 and 7: Lost lock between previous and current observation: cycle slip possible
    if (index(CSmethod,"LLI")/=0) then
        if (index(ObsCombine,"PC")/=0) then
            if ( (mod(LLI1,2)==1) .or. (mod(LLI2,2)==1) ) then
                Slip=1
            end if
        end if
        if ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0)) then
            if (mod(LLI1,2)==1)  then
                Slip=1
            end if
        end if
        if ( ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) .and. (freq_comb=='L2L3')) then
            if (mod(LLI1,2)==1)  then
                Slip=1
            end if
        end if
        if ( ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) .and. (freq_comb=='L1L2')) then
            if (mod(LLI2,2)==1)  then
                Slip=1
            end if        
        end if
        if ((index(ObsCombine,"P3")/=0) .or. (index(ObsCombine,"G3")/=0)) then
            if (mod(LLI2,2)==1)  then
                Slip=1
            end if
        end if
        if (Slip==1) then
            write(CSID,"(2I3,F6.2,2I5,I3)")  sta, PRN, ele, LLI1, LLI2, Slip
        end if
    end if

!    else  ! If single frequency, use L1P1 combination
    ! ======== P1L1 combination ==========
    if ( (index(CSmethod,"P1L1")/=0) .and. (index(ObsCombine,"PC")==0) ) then  ! only valid when single frequency
        if ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0)) then
            P=P1
            L=L1*c/f1
        elseif ( ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) .and. (freq_comb=='L2L3')) then
            P=P1
            L=L1*c/f1
        else
            P=P2
            L=L2*c/f2
        end if
        if (L==0.d0) return  ! If no phase data
        L1P1=L-P
        dL1P1=L1P1-CycleSlip(sta)%CS(PRN)%L1P1mean
        L1P1Threshold=0.d0
        if ((CycleSlip(sta)%CS(PRN)%LastWeek-Week)*604800.0d0+CycleSlip(sta)%CS(PRN)%LastSow-Sow<-Interval*10.d0) then
            Slip=2    ! If data missing too long or new satellite ascending
        elseif  (CycleSlip(sta)%CS(PRN)%arcLengthL1P1<3)  then
            Slip=1    ! If not enough epoches,accumulate data
        elseif (CycleSlip(sta)%CS(PRN)%arcLengthL1P1>=3)  then
            sigma2=CycleSlip(sta)%CS(PRN)%L1P1mean2-CycleSlip(sta)%CS(PRN)%L1P1mean**2
            ! Weight current sigma with an initial value
            sigma2=sigma2+(csL1P1init-sigma2)/CycleSlip(sta)%CS(PRN)%arcLengthL1P1
            L1P1Threshold=min(csL1P1slope*sqrt(sigma2),csL1P1max)
            if (abs(dL1P1)>L1P1Threshold) Slip=2
        end if

        !  initial cycle slip information if cycle slip happens
        if (Slip==2) then
            CycleSlip(sta)%CS(PRN)%arcLengthL1P1=0
            write(CSID,"(2I3,F6.2,2F10.4,I3)")  sta, PRN, ele, dL1P1, L1P1Threshold, Slip
            Slip=1
        end if

        ! Update L1P1 information
        CycleSlip(sta)%CS(PRN)%arcLengthL1P1=CycleSlip(sta)%CS(PRN)%arcLengthL1P1+1
        if (CycleSlip(sta)%CS(PRN)%arcLengthL1P1>int(csL1P1maxLength/Interval)) CycleSlip(sta)%CS(PRN)%arcLengthL1P1=int(csL1P1maxLength/Interval)
        CycleSlip(sta)%CS(PRN)%L1P1mean=CycleSlip(sta)%CS(PRN)%L1P1mean+(L1P1-CycleSlip(sta)%CS(PRN)%L1P1mean)/CycleSlip(sta)%CS(PRN)%arcLengthL1P1
        CycleSlip(sta)%CS(PRN)%L1P1mean2=CycleSlip(sta)%CS(PRN)%L1P1mean2+(L1P1**2-CycleSlip(sta)%CS(PRN)%L1P1mean2)/CycleSlip(sta)%CS(PRN)%arcLengthL1P1
        
    end if
    CycleSlip(sta)%CS(PRN)%LastWeek=Week   ! Also used in RTK, to delect whether to delete satellite, see Form_NEQ.f90
    CycleSlip(sta)%CS(PRN)%LastSow=Sow
    if (Slip==1) then
        CycleSlip(sta)%CS(PRN)%arcLengthSlip=0   ! Record this for reference satellite selection
    elseif (Slip==0) then
        CycleSlip(sta)%CS(PRN)%arcLengthSlip=CycleSlip(sta)%CS(PRN)%arcLengthSlip+1
        if (CycleSlip(sta)%CS(PRN)%arcLengthSlip>10) CycleSlip(sta)%CS(PRN)%arcLengthSlip=10
    end if

end subroutine