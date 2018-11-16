! ================   Cal Sat Clk  ===============
!
! PURPOSE:
!        Calculate the GNSS clock correction in the specific time,
!     using the precise clock by linear interpolation.
!
! INPUTS:
!       GPSweek    GPS week
!       GPSsec       GPS second
!       PRN            Satellite Number
!
! OUTPUT:
!       Sat_clk        Satellite clock correction, unit in second
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ===============  End of Header  ==============

subroutine Cal_Sat_Clk(GPSweek, GPSsec, PRN, Sat_clk)
use MOD_ClkData
use MOD_VAR
implicit none
    ! Intent in
    integer :: GPSweek
    real(8) :: GPSsec
    integer :: PRN
    ! Intent out
    real(8) :: Sat_Clk
    ! tgd
    real(8) :: TGD(2)
    
    real(8) :: dGPST(2)
    real(8) :: Clk1, Clk2,ClkVel
    
    dGPST=(ClkData%GPSweek(:)-GPSweek)*604800.0d0+ClkData%GPSsec(:)-GPSsec
    
     ! Linear Interpolation
        Clk1=ClkData%AS(PRN)%Clk(1)
        Clk2=ClkData%AS(PRN)%Clk(2)
        ClkVel=ClkData%AS(PRN)%ClkVel
        if ( (dabs(Clk1-9999.0d0)>0.01d0) .and. (dabs(Clk2-9999.0d0)>0.01d0) ) then ! Normal
            Sat_Clk=Clk1-dGPST(1)/(dGPST(2)-dGPST(1))*(Clk2-Clk1)
        else if ( (abs(Clk1-9999.0d0)>0.01d0) .and. (ClkVel/=9999.d0) ) then  ! If only Clk1 normal
            Sat_Clk=Clk1-dGPST(1)*ClkVel
        else if  ( (abs(Clk2-9999.0d0)>0.01d0) .and. (ClkVel/=9999.d0) ) then ! If only Clk2 normal
            Sat_Clk=Clk2-dGPST(2)*ClkVel
        else  ! If both Clk1 and Clk2 nan
            Sat_Clk=9999.0d0
        end if
        
        if ( ((PRN<=GNum+RNum) .or. (PRN>GNum+RNum+CNum)) .and. (Sat_Clk /=9999.d0)) then  ! For GPS/GLONASS/Galileo/QZSS
            if (index(ObsCombine,"PC")/=0) then
                    if ((PRN<=GNum) .or. (PRN>GNum+RNum+CNum+NumE)) then ! GPS/QZSS
                        if (freq_comb=='L1L3') then
                            Sat_Clk=Sat_Clk + 1.5457*DCBBSX(1,PRN)- 1.2606*DCBBSX(2,PRN)  ! P1P2==>B1B3 f2^2/(f1^2-f2^2)*T12-f3^2/(f1^2-f3^2)*T13
                        else  if (freq_comb=='L2L3') then
                            Sat_Clk=Sat_Clk +13.8010*DCBBSX(1,PRN)- 11.2553*DCBBSX(2,PRN)  ! P1P2==>B2B3 f2^2/(f1^2-f2^2)*T12+f2^2/(f2^2-f3^2)*T12-f3^2/(f2^2-f3^2)*T13
                        end if
                    elseif (PRN>GNum+RNum+CNum) then ! Galileo
                        if (freq_comb=='L1L3') then
                            Sat_Clk=Sat_Clk + 1.2606*DCBBSX(1,PRN)- 1.4220*DCBBSX(2,PRN)  ! P1P2==>B1B3 f2^2/(f1^2-f2^2)*T12-f3^2/(f1^2-f3^2)*T13
                        else  if (freq_comb=='L2L3') then
                            Sat_Clk=Sat_Clk -17.6593*DCBBSX(1,PRN)+19.9199*DCBBSX(2,PRN)  ! P1P2==>B2B3 f2^2/(f1^2-f2^2)*T12+f2^2/(f2^2-f3^2)*T12-f3^2/(f2^2-f3^2)*T13
                        end if                    
                    end if
            elseif ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0))then
                if ((PRN<=GNum) .or. (PRN>GNum+RNum+CNum+NumE)) then ! GPS/QZSS
                    Sat_Clk=Sat_Clk + 1.5457*DCBBSX(1,PRN)          ! P1P2==>P1  ! f1**2/(f1+f2)/(f1-f2)
                elseif (PRN>GNum+RNum+CNum) then ! Galileo
                    Sat_Clk=Sat_Clk + 1.2606*DCBBSX(1,PRN)          ! P1P2==>P1
                end if
            elseif ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0))then
                if ((PRN<=GNum) .or. (PRN>GNum+RNum+CNum+NumE)) then ! GPS/QZSS
                    Sat_Clk=Sat_Clk + 2.5457*DCBBSX(1,PRN)           ! P1P2==>P2  ! f1**2/(f1+f2)/(f1-f2)*DCBBSX(1,PRN)
                elseif (PRN>GNum+RNum+CNum) then ! Galileo
                    Sat_Clk=Sat_Clk + 2.2606*DCBBSX(1,PRN)          ! P1P2==>P2
                end if
            elseif ((index(ObsCombine,"P3")/=0) .or. (index(ObsCombine,"G3")/=0))then
                if ((PRN<=GNum) .or. (PRN>GNum+RNum+CNum+NumE)) then ! GPS/QZSS
                    Sat_Clk=Sat_Clk + 1.5457*DCBBSX(1,PRN) + DCBBSX(2,PRN)   ! P1P2==>P3
                elseif (PRN>GNum+RNum+CNum) then ! Galileo
                    Sat_Clk=Sat_Clk + 1.2606*DCBBSX(1,PRN) + DCBBSX(2,PRN)          ! P1P2==>P3
                end if
            end if
        elseif ((PRN>GNum+RNum) .and. (PRN<=GNum+RNum+CNum) .and. (Sat_Clk /=9999.d0)) then  ! 如果是BeiDou
            if (IorQ==0) then  ! I支路数据，需要加上IQ的差异
                if (SP3Time== -14.d0) then   ! 钟需要加上IQ的差异
                    Sat_Clk=Sat_Clk+DCBIQ(5,PRN-GNum-RNum)
                end if
                TGD=DCBIQ(3:4,PRN-GNum-RNum)
            else
                TGD=DCBIQ(1:2,PRN-GNum-RNum)  ! Q支路
            end if
            if (clktype==1) then
                ! 多星定轨卫星钟TGD改正，IGS/一期是以B1B2的Q支路为基准
                if (index(ObsCombine,"PC")/=0) then
                    if (freq_comb=='L1L3') then
                        Sat_Clk=Sat_Clk - 0.4565*TGD(1)-1.4872*TGD(2)  ! B1B2==>B1B3 f2^2/(f1^2-f2^2)*T12-f3^2/(f1^2-f3^2)*T13
                    else  if (freq_comb=='L2L3') then
                        Sat_Clk=Sat_Clk + 2.4872*TGD(1)+8.1018*TGD(2)  ! B1B2==>B2B3 f2^2/(f1^2-f2^2)*T12+f2^2/(f2^2-f3^2)*T12-f3^2/(f2^2-f3^2)*T13
                    end if
                elseif ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0))then   
                    Sat_Clk=Sat_Clk + 1.4872*(TGD(1) - TGD(2))         ! B1B2==>B1
                elseif ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) then
                    Sat_Clk=Sat_Clk + 2.4872*(TGD(1) -TGD(2))                      ! B1B2==>B2
                elseif((index(ObsCombine,"P3")/=0) .or. (index(ObsCombine,"G3")/=0)) then
                    Sat_Clk=Sat_Clk + 2.4872*TGD(1) - 1.4872*TGD(2)             ! B1B2==>B3
                end if
            elseif (clktype==2) then
                ! 多星定轨卫星钟，分米级是以B3的Q支路为基准
                if (index(ObsCombine,"PC")/=0) then
                    if (freq_comb=='L1L2') then
                        Sat_Clk=Sat_Clk - 2.4872*TGD(1) + 1.4872*TGD(2)    ! B3==>B1B2  f1^2/(f1^2-f2^2)*T13+f2^2/(f1^2-f2^2)*[T23]
                    elseif (freq_comb=='L1L3') then
                        Sat_Clk=Sat_Clk - 2.94368*TGD(1)               ! B3==>B1B3 f1^2/(f1^2-f3^2)*T13
                    else  if (freq_comb=='L2L3') then
                        Sat_Clk=Sat_Clk + 9.589*TGD(2)        ! B3==>B2B3 f2^2/(f2^2-f3^2)*T23
                    end if
                elseif ((index(ObsCombine,"P1")/=0) .or. (index(ObsCombine,"G1")/=0))then   
                    Sat_Clk=Sat_Clk - TGD(1)                       ! B3==>B1
                elseif ((index(ObsCombine,"P2")/=0) .or. (index(ObsCombine,"G2")/=0)) then
                    Sat_Clk=Sat_Clk - TGD(2)                      ! B3==>B2
                end if
            end if
        end if
        

!        ! X31a卫星钟
!    if (freq_comb=='L1L2') then
!        Sat_Clk=Sat_Clk - 1.4872*TGD(1) - TGD(2)   ! B3==>B1B2  f1^2/(f1^2-f2^2)*T13+f2^2/(f1^2-f2^2)*[T13-T12]
!    elseif (freq_comb=='L1L3') then
!        Sat_Clk=Sat_Clk - 2.94368*TGD(2)                    ! B3==>B1B3 f1^2/(f1^2-f3^2)*T13
!    else  if (freq_comb=='L2L3') then
!        Sat_Clk=Sat_Clk + 9.589*(TGD(2) - TGD(1))  ! B3==>B2B3 f2^2/(f1^2-f2^2)*T12+f2^2/(f2^2-f3^2)*T12-f3^2/(f2^2-f3^2)*T13
!    end if

    return

end subroutine
    