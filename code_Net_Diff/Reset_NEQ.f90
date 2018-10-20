! ============= Reset_NEQ===================== 
! PURPOSE:
!            This subroutine is to reset NEQ if reference satellite change.
!

! INPUTS:
!         DD                            double difference structure
!         NEQ                         normal equation structure
!         Epo_NEQ                 current epoch normal equation structure
!         RefSys                     reference satellite system, valid in tightly combined RTK
! OUTPUT:
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! =================End of header=================

subroutine Reset_NEQ(RefSat, DD, NEQ, Epo_NEQ, RefSys)
use MOD_DD
use MOD_NEQ
use MOD_Epo_NEQ
use MOD_constant
use MOD_Var
use MOD_FileID
use MOD_GLO_Fre
implicit none
    type(type_DD) :: DD
    type(type_NEQ) :: NEQ
    type(type_Epo_NEQ) :: Epo_NEQ
    integer :: RefSat(5)
    integer(1) :: RefSys
    ! Local variables
    integer :: i, j, sys, maxsys, PRN, freq, NDIFB, Ref_NDIFB
    real(8) :: ref_f1, ref_f2, ref_f3
    
    if (If_TC) then
        maxsys=1
    else
        maxsys=5
    end if

    j=0
    do i=ParaNum+1,NEQ%N+IonoNum
        if (i-ParaNum>2*SatNum) then
            j=2*SatNum
        elseif (i-ParaNum>SatNum) then
            j=SatNum
        end if
        if (i-j-ParaNum<=GNum) then
            sys=1
        elseif (i-j-ParaNum<=GNum+RNum) then
            sys=2
        elseif (i-j-ParaNum<=GNum+RNum+CNum) then
            sys=3
        elseif (i-j-ParaNum<=GNum+RNum+CNum+NumE) then
            sys=4
        elseif (i-j-ParaNum<=GNum+RNum+CNum+NumE+JNum) then
            sys=1
        elseif (i-j-ParaNum<=GNum+RNum+CNum+NumE+JNum+INum) then
            sys=5
        end if
        if ( (RefSat(sys) /= DD%RefSat(sys))  .and. (DD%RefSat(sys)/=0) .and. (RefSat(sys)/=0) ) then
            if (i<=SatNum+ParaNum) then
                if ((i/=RefSat(sys)+ParaNum) .and. (NEQ%dx(i)/=0.d0)) then   ! For L1
                    NEQ%dx(i)=NEQ%dx(i)-NEQ%dx(RefSat(sys)+ParaNum)  ! other satellite
                    NEQ%InvN(:,i)=NEQ%InvN(:,i)-NEQ%InvN(:,RefSat(sys)+ParaNum)
                    NEQ%InvN(i,:)=NEQ%InvN(i,:)-NEQ%InvN(RefSat(sys)+ParaNum,:)
                    if (ar_mode==3) then ! If fixed and hold mode
                        if ((NEQ%fixed_amb(RefSat(sys))/=0.99d0) .and. (NEQ%fixed_amb(i-ParaNum)/=0.99d0)) then
                            NEQ%fixed_amb(i-ParaNum)=NEQ%fixed_amb(i-ParaNum)-NEQ%fixed_amb(RefSat(sys))
                        else
                            NEQ%fixed_amb(i-ParaNum)=0.99d0
                        end if
                    end if
                end if

                if ((i/=RefSat(sys)+ParaNum) .and. If_Est_Iono .and. IonoNum>0) then
                    if  (Epo_NEQ%dx(i)==0.d0) then
                        cycle
                    end if
                    Epo_NEQ%dx(i)=Epo_NEQ%dx(i)-Epo_NEQ%dx(RefSat(sys)+ParaNum)  ! other satellite
                    Epo_NEQ%InvN(:,i)=Epo_NEQ%InvN(:,i)-Epo_NEQ%InvN(:,RefSat(sys)+ParaNum)
                    Epo_NEQ%InvN(i,:)=Epo_NEQ%InvN(i,:)-Epo_NEQ%InvN(RefSat(sys)+ParaNum,:)
                    if (ar_mode==3) then ! If fixed and hold mode
                        if ((Epo_NEQ%fixed_amb(RefSat(sys))/=0.99d0) .and. (Epo_NEQ%fixed_amb(i-ParaNum)/=0.99d0)) then
                            Epo_NEQ%fixed_amb(i-ParaNum)=Epo_NEQ%fixed_amb(i-ParaNum)-Epo_NEQ%fixed_amb(RefSat(sys))
                        else
                            Epo_NEQ%fixed_amb(i-ParaNum)=0.99d0
                        end if
                    end if
                end if
            elseif (i<=2*SatNum+ParaNum) then
                if ((i/=RefSat(sys)+SatNum+ParaNum) .and. (NEQ%dx(i)/=0.d0)) then   ! For L2
                    NEQ%dx(i)=NEQ%dx(i)-NEQ%dx(RefSat(sys)+SatNum+ParaNum)  ! other satellite
                    NEQ%InvN(:,i)=NEQ%InvN(:,i)-NEQ%InvN(:,RefSat(sys)+SatNum+ParaNum)
                    NEQ%InvN(i,:)=NEQ%InvN(i,:)-NEQ%InvN(RefSat(sys)+SatNum+ParaNum,:)
                    if (ar_mode==3) then ! If fixed and hold mode
                        if ((NEQ%fixed_amb(RefSat(sys)+SatNum)/=0.99d0) .and. (NEQ%fixed_amb(i-ParaNum)/=0.99d0)) then
                            NEQ%fixed_amb(i-ParaNum)=NEQ%fixed_amb(i-ParaNum)-NEQ%fixed_amb(RefSat(sys)+SatNum)
                        else
                            NEQ%fixed_amb(i-ParaNum)=0.99d0
                        end if
                    end if
                end if
                
                if ((i/=RefSat(sys)+SatNum+ParaNum) .and. If_Est_Iono .and. IonoNum>0 ) then
                    if  (Epo_NEQ%dx(i)==0.d0) then
                        cycle
                    end if
                    Epo_NEQ%dx(i)=Epo_NEQ%dx(i)-Epo_NEQ%dx(RefSat(sys)+SatNum+ParaNum)  ! other satellite
                    Epo_NEQ%InvN(:,i)=Epo_NEQ%InvN(:,i)-Epo_NEQ%InvN(:,RefSat(sys)+SatNum+ParaNum)
                    Epo_NEQ%InvN(i,:)=Epo_NEQ%InvN(i,:)-Epo_NEQ%InvN(RefSat(sys)+SatNum+ParaNum,:)
                    if (ar_mode==3) then ! If fixed and hold mode
                        if ((Epo_NEQ%fixed_amb(RefSat(sys)+SatNum)/=0.99d0) .and. (Epo_NEQ%fixed_amb(i-ParaNum)/=0.99d0)) then
                            Epo_NEQ%fixed_amb(i-ParaNum)=Epo_NEQ%fixed_amb(i-ParaNum)-Epo_NEQ%fixed_amb(RefSat(sys)+SatNum)
                        else
                            Epo_NEQ%fixed_amb(i-ParaNum)=0.99d0
                        end if
                    end if
                end if
            else
                if ((i/=RefSat(sys)+2*SatNum+ParaNum) .and. (Epo_NEQ%dx(i)/=0.d0)) then   ! For ionosphere parameter
                    Epo_NEQ%dx(i)=Epo_NEQ%dx(i)-Epo_NEQ%dx(RefSat(sys)+SatNum*2+ParaNum)  ! other satellite
                    Epo_NEQ%InvN(:,i)=Epo_NEQ%InvN(:,i)-Epo_NEQ%InvN(:,RefSat(sys)+SatNum*2+ParaNum)
                    Epo_NEQ%InvN(i,:)=Epo_NEQ%InvN(i,:)-Epo_NEQ%InvN(RefSat(sys)+SatNum*2+ParaNum,:)
                end if
            end if
        end if
    end do
    do sys=1,maxsys
        if ( (RefSat(sys) /= DD%RefSat(sys))  .and. (DD%RefSat(sys)/=0) .and. (RefSat(sys)/=0) ) then
!            write(LogID,'(I6,1X,I3,A4,I3,A15)') sys, DD%RefSat(sys),' -->',RefSat(sys),'ref sat change' 
            if (If_TC) then
                write(LogID,'(7X,I3,A4,I3,A16)') DD%RefSat(sys),' -->',RefSat(sys),'ref sat change'
            else
                if (sys==1) then
                    write(LogID,'(7X,A1,I2,A4,A1,I2,A16)') 'G', DD%RefSat(sys),' -->','G',RefSat(sys),'ref sat change'
                elseif (sys==2) then
                    write(LogID,'(7X,A1,I2,A4,A1,I2,A16)') 'R', DD%RefSat(sys)-GNum,' -->','R',RefSat(sys)-GNum,'ref sat change'
                elseif (sys==3) then
                    write(LogID,'(7X,A1,I2,A4,A1,I2,A16)') 'C', DD%RefSat(sys)-GNum-RNum,' -->','C',RefSat(sys)-GNum-RNum,'ref sat change'
                elseif (sys==4) then
                    write(LogID,'(7X,A1,I2,A4,A1,I2,A16)') 'E', DD%RefSat(sys)-GNum-RNum-CNum,' -->','E',RefSat(sys)-GNum-RNum-CNum,'ref sat change'
                end if
            end if
            NEQ%dx(DD%RefSat(sys)+ParaNum)=0.d0-NEQ%dx(RefSat(sys)+ParaNum)  ! old reference satellite
            NEQ%dx(DD%RefSat(sys)+SatNum+ParaNum)=0.d0-NEQ%dx(RefSat(sys)+SatNum+ParaNum)
            NEQ%InvN(:,DD%RefSat(sys)+ParaNum)= 0.d0 -NEQ%InvN(:,RefSat(sys)+ParaNum)
            NEQ%InvN(:,DD%RefSat(sys)+SatNum+ParaNum)= 0.d0 -NEQ%InvN(:,RefSat(sys)+SatNum+ParaNum)
            NEQ%InvN(DD%RefSat(sys)+ParaNum,:)= 0.d0 -NEQ%InvN(RefSat(sys)+ParaNum,:)
            NEQ%InvN(DD%RefSat(sys)+SatNum+ParaNum,:)= 0.d0 -NEQ%InvN(RefSat(sys)+SatNum+ParaNum,:)
            NEQ%dx(RefSat(sys)+ParaNum)=0.d0   ! new reference satellite
            NEQ%dx(RefSat(sys)+SatNum+ParaNum)=0.d0
            NEQ%InvN(:,RefSat(sys)+ParaNum)=0.d0
            NEQ%InvN(RefSat(sys)+ParaNum,:)=0.d0
            NEQ%InvN(:,RefSat(sys)+SatNum+ParaNum)=0.d0
            NEQ%InvN(RefSat(sys)+SatNum+ParaNum,:)=0.d0
                    
            if (If_Est_Iono .and. IonoNum>0) then
                Epo_NEQ%dx(DD%RefSat(sys)+ParaNum)=0.d0-Epo_NEQ%dx(RefSat(sys)+ParaNum)  ! old reference satellite
                Epo_NEQ%dx(DD%RefSat(sys)+SatNum+ParaNum)=0.d0-Epo_NEQ%dx(RefSat(sys)+SatNum+ParaNum)
                Epo_NEQ%dx(DD%RefSat(sys)+2*SatNum+ParaNum)=0.d0-Epo_NEQ%dx(RefSat(sys)+2*SatNum+ParaNum)
                Epo_NEQ%InvN(:,DD%RefSat(sys)+ParaNum)= 0.d0 -Epo_NEQ%InvN(:,RefSat(sys)+ParaNum)
                Epo_NEQ%InvN(:,DD%RefSat(sys)+SatNum+ParaNum)= 0.d0 -Epo_NEQ%InvN(:,RefSat(sys)+SatNum+ParaNum)
                Epo_NEQ%InvN(:,DD%RefSat(sys)+2*SatNum+ParaNum)= 0.d0 -Epo_NEQ%InvN(:,RefSat(sys)+2*SatNum+ParaNum)
                Epo_NEQ%InvN(DD%RefSat(sys)+ParaNum,:)= 0.d0 -Epo_NEQ%InvN(RefSat(sys)+ParaNum,:)
                Epo_NEQ%InvN(DD%RefSat(sys)+SatNum+ParaNum,:)= 0.d0 -Epo_NEQ%InvN(RefSat(sys)+SatNum+ParaNum,:)
                Epo_NEQ%InvN(DD%RefSat(sys)+2*SatNum+ParaNum,:)= 0.d0 -Epo_NEQ%InvN(RefSat(sys)+2*SatNum+ParaNum,:)
                Epo_NEQ%dx(RefSat(sys)+ParaNum)=0.d0   ! new reference satellite
                Epo_NEQ%dx(RefSat(sys)+SatNum+ParaNum)=0.d0
                Epo_NEQ%dx(RefSat(sys)+2*SatNum+ParaNum)=0.d0
                Epo_NEQ%InvN(:,RefSat(sys)+ParaNum)=0.d0
                Epo_NEQ%InvN(RefSat(sys)+ParaNum,:)=0.d0
                Epo_NEQ%InvN(:,RefSat(sys)+SatNum+ParaNum)=0.d0
                Epo_NEQ%InvN(RefSat(sys)+SatNum+ParaNum,:)=0.d0
                Epo_NEQ%InvN(:,RefSat(sys)+2*SatNum+ParaNum)=0.d0
                Epo_NEQ%InvN(RefSat(sys)+2*SatNum+ParaNum,:)=0.d0
            end if

            if (ar_mode==3) then ! If fixed and hold mode
                if (NEQ%fixed_amb(RefSat(sys))/=0.99d0) then
                    NEQ%fixed_amb(DD%RefSat(sys))=0.d0-NEQ%fixed_amb(RefSat(sys))    ! new fixed ambiguity of old reference satellite
                else
                    NEQ%fixed_amb(DD%RefSat(sys))=0.99d0
                end if
                if (NEQ%fixed_amb(RefSat(sys)+SatNum)/=0.99d0) then
                    NEQ%fixed_amb(DD%RefSat(sys)+SatNum)=0.d0-NEQ%fixed_amb(RefSat(sys)+SatNum)
                else
                    NEQ%fixed_amb(DD%RefSat(sys)+SatNum)=0.99d0
                end if
                NEQ%fixed_amb(RefSat(sys))=0.99d0    ! new fixed ambiguity of new reference satellite
                NEQ%fixed_amb(RefSat(sys)+SatNum)=0.99d0

                if (If_Est_Iono .and. IonoNum>0) then
                    if (Epo_NEQ%fixed_amb(RefSat(sys))/=0.99d0) then
                        Epo_NEQ%fixed_amb(DD%RefSat(sys))=0.d0-Epo_NEQ%fixed_amb(RefSat(sys))    ! new fixed ambiguity of old reference satellite
                    else
                        Epo_NEQ%fixed_amb(DD%RefSat(sys))=0.99d0
                    end if
                    if (Epo_NEQ%fixed_amb(RefSat(sys)+SatNum)/=0.99d0) then
                        Epo_NEQ%fixed_amb(DD%RefSat(sys)+SatNum)=0.d0-Epo_NEQ%fixed_amb(RefSat(sys)+SatNum)
                    else
                        Epo_NEQ%fixed_amb(DD%RefSat(sys)+SatNum)=0.99d0
                    end if
                    Epo_NEQ%fixed_amb(RefSat(sys))=0.99d0    ! new fixed ambiguity of new reference satellite
                    Epo_NEQ%fixed_amb(RefSat(sys)+SatNum)=0.99d0
                end if
            end if
        end if
    end do

    ! ******************** For DISB in tightly combined RTK ***********************
    if (If_TC .and. DD%RefSys/=RefSys) then
        If_Fix_DISB=.true.
        if (RefSys==1) then  ! GPS/QZSS
            ref_f1=f_L1
            ref_f2=f_L2
            ref_f3=f_L5
            Ref_NDIFB=4+INT_SystemUsed(1)*4  ! There is bug if no GPS but use QZSS, should seperate GPS and QZSS
        elseif (RefSys==2) then  ! GLONASS
            freq=Fre_Chann(PRN-GNum)  ! Actually, GLONASS will not choose as reference satellite in tightly combined RTK
            ref_f1=(1602.0d0+freq*0.5625d0)*1.0D6   ! f1=(1602.0d0+K*0.5625d0)*1.0d6
            ref_f2=(1246.0d0+freq*0.4375d0)*1.0D6
        elseif (RefSys==3) then  ! BeiDou
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
            Ref_NDIFB=4+INT_SystemUsed(1)*4+INT_SystemUsed(3)*4 ! End index of DISB
        elseif (RefSys==4) then ! GALILEO
            if (freq_comb=='L1L2') then   ! E1 E5a
                ref_f1=f_E1
                ref_f2=f_E5
                f3=f_E5b
            elseif (freq_comb=='L1L3') then   ! E1 E5b
                ref_f1=f_E1
                ref_f2=f_E5b
            elseif (freq_comb=='L2L3') then   ! E5a E5b
                ref_f1=f_E5a
                ref_f2=f_E5b
            end if
            Ref_NDIFB=4+INT_SystemUsed(1)*4+INT_SystemUsed(3)*4+INT_SystemUsed(4)*4 ! End index of DISB
        elseif (RefSys==5) then   ! IRNSS
            ref_f1=f_L1
            ref_f2=f_S
            Ref_NDIFB=4+INT_SystemUsed(1)*4+INT_SystemUsed(3)*4+INT_SystemUsed(4)*4+INT_SystemUsed(6)*4 ! End index of DISB
        end if
        do sys=1,5
            if (sys==RefSys) cycle
            if (sys==1) then  ! GPS/QZSS
                if ((.not.(SystemUsed(sys))) .and. (.not.(SystemUsed(5))) ) cycle  ! If no GPS and QZSS
                NDIFB=4+INT_SystemUsed(1)*4  ! There is bug if no GPS but use QZSS, should seperate GPS and QZSS
                f1=f_L1
                f2=f_L2
                f3=f_L5
            elseif (sys==2) then  ! GLONASS
                if (.not.(SystemUsed(sys))) cycle   ! If no GLONASS
            elseif (sys==3) then  ! BeiDou
                if (.not.(SystemUsed(sys))) cycle   ! If no BeiDou
                NDIFB=4+INT_SystemUsed(1)*4+INT_SystemUsed(3)*4 ! End index of DISB
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
            elseif (sys==4) then ! GALILEO
                if (.not.(SystemUsed(sys))) cycle   ! If no Galileo
                NDIFB=4+INT_SystemUsed(1)*4+INT_SystemUsed(3)*4+INT_SystemUsed(4)*4 ! End index of DISB
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
                if (.not.(SystemUsed(6))) cycle   ! If no IRNSS
                NDIFB=4+INT_SystemUsed(1)*4+INT_SystemUsed(3)*4+INT_SystemUsed(4)*4+INT_SystemUsed(6)*4 ! End index of DISB
                f1=f_L1
                f2=f_S
            end if
            if (sys/=2) then  ! If not GLONASS
                j=1
            elseif (sys==2) then
                j=GloFreqNum  ! If GLONASS
            end if
            do i=1, j
                if (sys/=2) then
                    NDIFB=NDIFB
                elseif (sys==2) then
                    NDIFB=ParaNum-GloFreqNum*4+i*4  ! End index of GLONASS DISB
                end if
                NEQ%dx(NDIFB-3:NDIFB)=NEQ%dx(NDIFB-3:NDIFB) - NEQ%dx(Ref_NDIFB-3:Ref_NDIFB) ! old reference system
                NEQ%InvN(:,NDIFB-3:NDIFB)= NEQ%InvN(:,NDIFB-3:NDIFB) - NEQ%InvN(:,Ref_NDIFB-3:Ref_NDIFB)
                NEQ%InvN(NDIFB-3:NDIFB,:)= NEQ%InvN(NDIFB-3:NDIFB,:) - NEQ%InvN(Ref_NDIFB-3:Ref_NDIFB,:)
                if (NEQ%InvN(NDIFB-1,NDIFB-1)>0.d0) then   ! Add random walk due to frequency difference, 
                    NEQ%InvN(NDIFB-1,NDIFB-1)=NEQ%InvN(NDIFB-1,NDIFB-1)+(50.d0*(c/ref_f1-c/f1))**2  ! 
                    NEQ%InvN(NDIFB,NDIFB)=NEQ%InvN(NDIFB,NDIFB)+(50.d0*(c/ref_f2-c/f2))**2  ! Max ref. sat.  amb. change is set as 50 cycle
                end if
                if (If_Est_Iono .and. IonoNum>0) then
                    Epo_NEQ%dx(NDIFB-3:NDIFB)=Epo_NEQ%dx(NDIFB-3:NDIFB) - Epo_NEQ%dx(Ref_NDIFB-3:Ref_NDIFB)  ! old reference system
                    Epo_NEQ%InvN(:,NDIFB-3:NDIFB)= Epo_NEQ%InvN(:,NDIFB-3:NDIFB) - Epo_NEQ%InvN(:,Ref_NDIFB-3:Ref_NDIFB)
                    Epo_NEQ%InvN(NDIFB-3:NDIFB,:)= Epo_NEQ%InvN(NDIFB-3:NDIFB,:) - Epo_NEQ%InvN(Ref_NDIFB-3:Ref_NDIFB,:)
                    if (Epo_NEQ%InvN(NDIFB-1,NDIFB-1)>0.d0) then   ! Add random walk due to frequency difference
                        Epo_NEQ%InvN(NDIFB-1,NDIFB-1)=Epo_NEQ%InvN(NDIFB-1,NDIFB-1)+(50.d0*(c/ref_f1-c/f1))**2
                        Epo_NEQ%InvN(NDIFB,NDIFB)=Epo_NEQ%InvN(NDIFB,NDIFB)+(50.d0*(c/ref_f2-c/f2))**2
                    end if
                end if
            end do
        end do
        NEQ%dx(Ref_NDIFB-3:Ref_NDIFB)=0.d0   ! New reference system
        NEQ%InvN(:,Ref_NDIFB-3:Ref_NDIFB)=0.d0
        NEQ%InvN(Ref_NDIFB-3:Ref_NDIFB,:)=0.d0
        if (If_Est_Iono .and. IonoNum>0) then
            Epo_NEQ%dx(Ref_NDIFB-3:Ref_NDIFB)=0.d0
            Epo_NEQ%InvN(:,Ref_NDIFB-3:Ref_NDIFB)= 0.d0
            Epo_NEQ%InvN(Ref_NDIFB-3:Ref_NDIFB,:)= 0.d0
        end if
        DD%RefSys=RefSys
    end if

    
    ! ******************** For loosely combined GLONASS DISB change ************************
    if (.not.(If_TC) .and. GloFreqNum>0) then
        sys=2
        if ( Fre_Chann(RefSat(sys)-GNum) /= Fre_Chann(DD%RefSat(sys)-GNum)) then ! If  reference satellite change
            Ref_NDIFB=ParaNum-GloFreqNum*4+Fre_Chann(RefSat(sys)-GNum)*4+8*4
            do i=1,GloFreqNum
                NDIFB=ParaNum-GloFreqNum*4+i*4  ! End index of GLONASS DISB
                if (NEQ%dx(NDIFB)==0.d0) cycle  ! If the new reference satellite frequency or no satellite in this frequency
                if ( DD%RefSat(sys)/=0 .and. RefSat(sys)/=0 ) then 
                    NEQ%dx(NDIFB-3:NDIFB)=NEQ%dx(NDIFB-3:NDIFB)-NEQ%dx(Ref_NDIFB-3:Ref_NDIFB)  ! other satellite
                    NEQ%InvN(:,NDIFB-3:NDIFB)=NEQ%InvN(:,NDIFB-3:NDIFB)-NEQ%InvN(:,Ref_NDIFB-3:Ref_NDIFB)
                    NEQ%InvN(NDIFB-3:NDIFB,:)=NEQ%InvN(NDIFB-3:NDIFB,:)-NEQ%InvN(Ref_NDIFB-3:Ref_NDIFB,:)
                    if (If_Est_Iono .and. IonoNum>0) then
                        Epo_NEQ%dx(NDIFB-3:NDIFB)=Epo_NEQ%dx(NDIFB-3:NDIFB)-Epo_NEQ%dx(Ref_NDIFB-3:Ref_NDIFB)  ! other satellite
                        Epo_NEQ%InvN(:,NDIFB-3:NDIFB)=Epo_NEQ%InvN(:,NDIFB-3:NDIFB)-Epo_NEQ%InvN(:,Ref_NDIFB-3:Ref_NDIFB)
                        Epo_NEQ%InvN(NDIFB-3:NDIFB,:)=Epo_NEQ%InvN(NDIFB-3:NDIFB,:)-Epo_NEQ%InvN(Ref_NDIFB-3:Ref_NDIFB,:)
                    end if
                end if
            end do
            if ( DD%RefSat(sys)/=0 .and. RefSat(sys)/=0 ) then
                NEQ%dx(Ref_NDIFB-3:Ref_NDIFB)=0.d0   ! new reference satellite
                NEQ%InvN(:,Ref_NDIFB-3:Ref_NDIFB)=0.d0
                NEQ%InvN(Ref_NDIFB-3:Ref_NDIFB,:)=0.d0
                if (If_Est_Iono .and. IonoNum>0) then
                    Epo_NEQ%dx(Ref_NDIFB-3:Ref_NDIFB)=0.d0   ! new reference satellite
                    Epo_NEQ%InvN(:,Ref_NDIFB-3:Ref_NDIFB)=0.d0
                    Epo_NEQ%InvN(Ref_NDIFB-3:Ref_NDIFB,:)=0.d0
                end if
            end if
        end if
    end if

    if (ADmethod=='LS') then
        call InvSqrt(NEQ%InvN, NEQ%N, NEQ%Nbb)
        NEQ%U=MATMUL(NEQ%Nbb, NEQ%dx)   ! New Nbb and U is needed
    end if

!        if (If_Est_Iono .and. IonoNum>0 .and. ADmethod=='LS') then  ! This is only for Least Square method, but we change to transformed Kalman Filter
!            call InvSqrt(Epo_NEQ%InvN, Epo_NEQ%N, Epo_NEQ%Nbb)
!            Epo_NEQ%U=MATMUL(Epo_NEQ%Nbb, Epo_NEQ%dx)   ! New Nbb and U is needed
!        end if

    return
end subroutine