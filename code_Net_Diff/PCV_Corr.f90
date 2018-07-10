! ====================  PCV_Corr  ===================
!
! PURPOSE:
!         Satellite or station receiver PCV(Phase Center Variation) correction.
!
!  REFERENCES:
!        http://www.epncb.oma.be/ftp/station/general/antex14.txt
!        http://www.navipedia.net/index.php/Antenna_Phase_Centre
!        http://www.navipedia.net/index.php/Satellite_Antenna_Phase_Centre
!        http://www.navipedia.net/index.php/Receiver_Antenna_Phase_Centre
! 
! INPUTS:
!          PRN          PRN
!          Ele            Satellite zenith elevation
!          Azi            Azimuth
!
! OUTPUTS:
!          PCV           Phase Center Variation in observed direction
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ====================  End of header  ======================

subroutine PCV_Corr(PRN, Ele, Azi, PCV)       
use MOD_Ant
use MOD_constant
implicit none
    ! Intent in
    integer :: PRN
    real(8) :: Ele, Azi
    ! Intent out
    real(8) :: PCV(2)
    ! Local variables
    integer :: i
    real(8) :: PCV11, PCV12, PCV21, PCV22, PCV1, PCV2, dazi

    PCV11=0.d0; PCV12=0.d0; PCV21=0.d0; PCV22=0.d0
    if (Ant(PRN)%nazi==0) then  ! if non-azimuth-dependent
        call Get_PCV_Ele(PRN, Ele, 1 ,1, PCV1) ! Ferq 1
        call Get_PCV_Ele(PRN, Ele, 1 , 2, PCV2) ! Freq 2
        PCV=(f1*f1*PCV1-f2*f2*PCV2)/(f1+f2)/(f1-f2)
    else
        do i=1,Ant(PRN)%nazi
            dazi=Azi-Ant(PRN)%Azimuth(i)
            if ( (dazi>=0.d0) .and. (dazi<Ant(PRN)%dazi) ) then
                call Get_PCV_Ele(PRN, Ele, i+1, 1, PCV11) ! Ferq 1
                call Get_PCV_Ele(PRN, Ele, i+1, 2, PCV21) ! Freq 2
            else if ( (dazi<=0.d0) .and. (dazi>-Ant(PRN)%dazi) )  then
                call Get_PCV_Ele(PRN, Ele, i+1, 1, PCV12) ! Ferq 1
                call Get_PCV_Ele(PRN, Ele, i+1, 2, PCV22) ! Freq 2
                exit
            end if
        end do
        PCV1=PCV12+dazi/Ant(PRN)%dazi*(PCV12-PCV11)
        PCV2=PCV22+dazi/Ant(PRN)%dazi*(PCV22-PCV21)
        PCV(1)=PCV1
        PCV(2)=PCV2
    end if
end subroutine


! ================= Get_PCV_Ele ===================
! Function:
!             Get PCV according to elevation in a specific azimuth
! Inputs:
!            PRN
!             Ele              elevation
!             or_azi          the order of specific azimuth
!             freq             frequency channel
! Output:
!            PCV            PCV in the specific elevation
! 
! Written by: Yize Zhang
! ===================================================

subroutine Get_PCV_Ele(PRN, Ele, or_azi, freq, PCV)
use MOD_Ant
    implicit none
    real(8) :: Ele, PCV
    integer :: PRN, or_azi, freq
    ! Local Variables
    integer :: i
    real(8) :: PCV1, PCV2, dzen

    PCV1=0.d0
    PCV2=0.d0
    do i=1,Ant(PRN)%nzen
        dzen=Ele-Ant(PRN)%Zenith(i)
        if ( (dzen>=0.d0) .and. (dzen<Ant(PRN)%dzen) ) then
            PCV1=Ant(PRN)%PCV(i,or_azi, freq)
        else if ( (dzen<=0.d0) .and. (dzen>-Ant(PRN)%dzen) ) then
            PCV2=Ant(PRN)%PCV(i,or_azi, freq)
            exit
        end if
    end do
    PCV=PCV2+dzen/Ant(PRN)%dzen*(PCV2-PCV1)
end subroutine