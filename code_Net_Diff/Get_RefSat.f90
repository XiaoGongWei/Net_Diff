! ======================== Get_Refsat ========================
! PURPOSE:
!         To check and get the reference satellite in double difference.
!
! If it is in the first epoch, then the reference satellite can be set as the one
!       with hightest elevation. 
!If it is not the first epoch, which means that a  reference satellite already exists,
!       then we shoule check the reference satellite. If the reference satellite not 
!       appear in this epoch, or its zenith elevation is below 30 degree, or the L1 or
!       L2 carrier phase data not exists, then set a new satellite as the reference.
!
! INPUTS:
!          SD             station difference structure
!         RefSat0      reference satellite in DD structure
!                                 In the first epoch, it is 0, then we can set the 
!                                 highest elevation satellite as the reference satellite.
! OUTPUT:
!          RefSat        reference satellite in this epoch
!                                 
! WRITTEN BY: Yize Zhang
! =====================================================

subroutine Get_RefSat(SD,RefSat0, RefSat)
use MOD_SD
use MOD_VAR
use MOD_CycleSlip
implicit none
    ! Arguments:
    type(type_SD)  ::   SD
    integer             ::   RefSat0(5)   ! For each system
    integer             ::   RefSat(5)
    ! Local variables
    integer(1)        ::   i, sys
    real(8)              ::   maxele

    
    RefSat=0
    do sys=1,5  ! For each system, select a referen e satellite
        if (sys==1) then
            if ((.not.(SystemUsed(sys))) .and. (.not.(SystemUsed(5))) ) cycle  ! If no GPS and QZSS
        elseif (sys==5) then
            if (.not.(SystemUsed(6))) cycle   ! If no IRNSS
        else
            if (.not.(SystemUsed(sys))) cycle
        end if
        if (RefSat0(sys)==0) then   ! if the first epoch
            200    maxEle=30.d0
            do i=1, SD%PRNS
                ! PC
                if (SD%Sys(i)/=Sys) cycle
                if ((SD%PRN(i)>GNum+RNum+CNum+NumE+2) .and. (sys==1)) cycle  ! Not choose QZSS3, due to half cycle slip??
                if ( (Combination(2)==.false.) .and. (SD%Ele(i)>maxEle) .and. (SD%P1(i)/=0.d0)  )  then
                    maxele=SD%Ele(i)
                    RefSat(sys)=SD%PRN(i)
                ! L1L2, in DD
                elseif ( (Combination(2)==.true.) .and.(SD%Ele(i)>maxEle) ) then
                    if ( ((a1/=0.d0) .and. (b2/=0.d0)) .or. ((a2/=0.d0) .and. (b1/=0.d0)) ) then   ! Dual frequency
                        if ( (SD%WL(i)/=0.d0) .and. (SD%W4(i)/=0.d0) ) then
                            maxele=SD%Ele(i)
                            RefSat(sys)=SD%PRN(i)
                        end if
                    elseif ( (SD%WL(i)/=0.d0 .or. SD%W4(i)/=0.d0) ) then  ! L1 or L2 only. For long baseline, it use WL+L1L2, so we use SD%WL here
                        maxele=SD%Ele(i)
                        RefSat(sys)=SD%PRN(i)
!                    elseif ( If_Est_Iono .and. SD%EWL(i)/=0.d0 ) then  ! for triple-frequency and long baseline
!                        maxele=SD%Ele(i)
!                        RefSat(sys)=SD%PRN(i)
                    end if
                end if
            end do
        else      ! If the reference satellite already exsits
!            do i=1, SD%PRNS
!                  ! Check if elevation of reference satellite > 40deg and data available and no cycle slip
!                if ((SD%PRN(i)==RefSat0(sys)) .and. (SD%Ele(i)>40.d0) ) then
!                    if ( (Combination(2)==.false.).and. (SD%P1(i)/=0.d0) ) then   ! PC
!                        RefSat(sys)=RefSat0(sys)
!                        exit
!                    else if ( (Combination(2)==.true.) ) then  ! L1L2, in DD
!                        if ( ((a1/=0.d0) .and. (b2/=0.d0)) .or. ((a2/=0.d0) .and. (b1/=0.d0)) ) then   ! Dual frequency
!                            if ( (SD%L1(i)/=0.d0) .and. (SD%L2(i)/=0.d0) .and. (CycleSlip(1)%Slip(RefSat0(sys))==0) .and. (CycleSlip(2)%Slip(RefSat0(sys))==0) ) then
!                                RefSat(sys)=RefSat0(sys)
!                            end if
!                        elseif ( (SD%WL(i)/=0.d0 .or. SD%W4(i)/=0.d0) .and. (CycleSlip(1)%Slip(RefSat0(sys))==0) .and. (CycleSlip(2)%Slip(RefSat0(sys))==0) ) then  ! L1 or L2 only
!                            RefSat(sys)=RefSat0(sys)
!!                        elseif ( If_Est_Iono .and. (SD%WL(i)/=0.d0 .or. SD%W4(i)/=0.d0) .and. (CycleSlip(1)%Slip(RefSat0(sys))==0) .and. (CycleSlip(2)%Slip(RefSat0(sys))==0) ) then  ! L1 or L2 only
!!                            RefSat(sys)=RefSat0(sys)
!                        end if
!                        exit
!                    end if
!                end if
!            end do
!            if (RefSat(sys)>0) cycle
!            ! If not found the reference satellite or it's elevation is below 40 degrees or cycle slip occurs in reference satellite
!            ! it means that a new reference satellite is found
            maxEle=30.d0
            do i=1, SD%PRNS
                ! PC
                if (SD%Sys(i)/=Sys) cycle
                if ((SD%PRN(i)>GNum+RNum+CNum+NumE+2) .and. (sys==1)) cycle  ! Not choose QZSS3, due to half cycle slip??
                if ( (Combination(2)==.false.) .and. (SD%Ele(i)>maxEle) .and. (SD%P1(i)/=0.d0)  )  then
                    maxele=SD%Ele(i)
                    RefSat(sys)=SD%PRN(i)
                ! L1L2, in DD
                elseif ( (Combination(2)==.true.) .and.(SD%Ele(i)>maxEle) ) then
                    if ( ((a1/=0.d0) .and. (b2/=0.d0)) .or. ((a2/=0.d0) .and. (b1/=0.d0)) ) then   ! Dual frequency
                        if ( (SD%L1(i)/=0.d0) .and. (SD%L2(i)/=0.d0) .and. (CycleSlip(1)%Slip(SD%PRN(i))==0) .and. (CycleSlip(2)%Slip(SD%PRN(i))==0) ) then
                            if (ar_mode==0) then  ! If float solution
                                maxele=SD%Ele(i)
                                RefSat(sys)=SD%PRN(i)
                            elseif (par_PRN(SD%PRN(i))/=1 .and. par_PRN_Epo(SD%PRN(i))/=1) then   ! If continuous solution, don't use the unfixed satellite
                                maxele=SD%Ele(i)
                                RefSat(sys)=SD%PRN(i)
                            end if
                        end if
                    elseif ( (SD%WL(i)/=0.d0 .or. SD%W4(i)/=0.d0) .and. (CycleSlip(1)%Slip(SD%PRN(i))==0) .and. (CycleSlip(2)%Slip(SD%PRN(i))==0) ) then  ! L1 or L2 only
                        if (ar_mode==0) then  ! If float solution
                            maxele=SD%Ele(i)
                            RefSat(sys)=SD%PRN(i)
                        elseif (par_PRN(SD%PRN(i))/=1 .and. par_PRN_Epo(SD%PRN(i))/=1) then   ! If continuous solution, don't use the unfixed satellite
                            maxele=SD%Ele(i)
                            RefSat(sys)=SD%PRN(i)
                        end if
!                    elseif ( If_Est_Iono .and. (SD%WL(i)/=0.d0 .or. SD%W4(i)/=0.d0) .and. (CycleSlip(1)%Slip(SD%PRN(i))==0) .and. (CycleSlip(2)%Slip(SD%PRN(i))==0) ) then  ! L1 or L2 only
!                        maxele=SD%Ele(i)
!                        RefSat(sys)=SD%PRN(i)
                    end if
                end if
            end do
            ! If cycle slip for all satellites 
           if (RefSat(sys)==0) goto 200
        end if
    end do

    return
end subroutine