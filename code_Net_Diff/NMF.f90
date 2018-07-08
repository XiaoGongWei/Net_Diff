!===================   Map_NMF   ====================
!
! PURPOSE:
!      Calculate the Slant Total Delay of tropsphere using Niell Mapping Function
!
! REFERENCES:
!   http://www.navipedia.net/index.php/Mapping_of_Niell
!   Niell, A., 1996. Global mapping functions for the atmosphere delay at radio wavelengths. 
!       Journal of Geophysical Research. 101, pp. 3227-3246.
!
! INPUTS:
!      Ele       Satllite elevation(degree)
!      B          Latitude(degree)
!      H          Height(m)
!      DOY     Day of year
!
! OUTPUTS:
!      MFh     Wet mapping fator of tropsphere parameter
!      MFw    Hydrostatic mapping fator of tropsphere parameter
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ==================  End of header  ==================

subroutine Map_NMF(Ele, B, H, DOY, MFh, MFw)
implicit none
    ! Intent in
    real(8) :: Ele
    real(8) :: B, H, absB
    integer :: DOY
    ! Intent out
    real(8) :: STD
    
    integer :: doymin
    real(8), parameter :: pi=3.141592653589d0
    real(8) :: k, ah, bh, ch, aw,bw, cw, m
    real(8) :: aht, bht, cht, MFh, MFw
    
    if (B>0.0) then
        doymin=28 ! mininum day in a year
    else
        doymin=211
    end if

    k=dcos(2.0d0*pi*(DOY-doymin)/365.25d0)
    absB=dabs(B);
    if (absB<15.0d0) then
        ah=1.2769934D-3
        bh=2.9153695D-3
        ch=62.610505D-3
        aw=5.8021897D-4
        bw=1.4275268D-3
        cw=4.3472961D-2
    elseif (absB<30d0) then
        m=(absB-15.0d0)/15.0d0
        ah=1.2769934D-3+m*(1.2683230d0-1.2769934d0)*1D-3-(0.d0+m*(1.2709626d0-0.d0))*1D-5*k
        bh=2.9153695D-3+m*(2.9152299d0-2.9153695d0)*1D-3-(0.d0+m*(2.1414979d0-0.d0))*1D-5*k
        ch=62.610505D-3+m*(62.837393d0-62.610505d0)*1D-3-(0.d0+m*(9.012840d0-0.d0))*1D-5*k
        aw=5.8021897D-4+m*(5.6794847d0-5.8021897d0)*1D-4
        bw=1.4275268D-3+m*(1.5138625d0-1.4275268d0)*1D-3
        cw=4.3472961D-2+m*(4.6729510d0-4.3472961d0)*1D-2
    elseif (absB<45.0d0) then
        m=(absB-30.0d0)/15.0d0
        ah=1.2683230D-3+m*(1.2465397d0-1.2683230d0)*1D-3-(1.2709626d0+m*(2.65236662d0-1.2709626d0))*1D-5*k
        bh=2.9152299D-3+m*(2.9288445d0-2.9152299d0)*1D-3-(2.1414979d0+m*(3.0160779d0-2.1414979d0))*1D-5*k
        ch=62.837393D-3+m*(63.721774d0-62.837393d0)*1D-3-(9.012840d0+m*(4.349703d0-9.012840d0))*1D-5*k
        aw=5.6794847D-4+m*(5.8118019d0-5.6794847d0)*1D-4
        bw=1.5138625D-3+m*(1.4572752d0-1.5138625d0)*1D-3
        cw=4.6729510D-2+m*(4.3908931d0-4.6729510d0)*1D-2
    elseif (absB<60.0d0) then
        m=(absB-45.0d0)/15.0d0
        ah=1.2465397D-3+m*(1.2196049d0-1.2465397d0)*1D-3-(2.65236662d0+m*(3.4000452d0-2.65236662d0))*1D-5*k
        bh=2.9288445D-3+m*(2.9022565d0-2.9288445d0)*1D-3-(3.0160779d0+m*(7.2562722d0-3.0160779d0))*1D-5*k
        ch=63.721774D-3+m*(63.824265d0-63.721774d0)*1D-3-(4.349703d0+m*(84.795348d0-4.349703d0))*1D-5*k
        aw=5.8118019D-4+m*(5.9727542d0-5.8118019d0)*1D-4
        bw=1.4572752D-3+m*(1.5007428d0-1.4572752d0)*1D-3
        cw=4.3908931D-2+m*(4.4626982-4.3908931)*1D-2
    elseif (absB<75.0d0) then
        m=(absB-60.0d0)/15.0d0
        ah=1.2196049D-3+m*(1.2045996d0-1.2196049d0)*1D-3-(3.4000452d0+m*(4.1202191d0-3.4000452d0))*1D-5*k
        bh=2.9022565D-3+m*(2.9024912d0-2.9022565d0)*1D-3-(7.2562722d0+m*(11.723375d0-7.2562722d0))*1D-5*k
        ch=63.824265D-3+m*(64.258455d0-63.824265d0)*1D-3-(84.795348d0+m*(170.37206d0-84.795348d0))*1D-5*k
        aw=5.9727542D-4+m*(6.1641693d0-5.9727542d0)*1D-4
        bw=1.5007428D-3+m*(1.7599082d0-1.5007428d0)*1D-3
        cw=4.4626982D-2+m*(5.4736038d0-4.4626982d0)*1D-2
    else
        ah=1.2045996D-3-4.1202191D-5*k
        bh=2.9024912D-3-11.723375D-5*k
        ch=64.258455D-3-170.37206D-5*k
        aw=6.1641693D-4
        bw=1.7599082D-3
        cw=5.4736038D-2
    end if

    aht=2.53D-5;bht=5.49D-3;cht=1.14D-3
    ! Dry and wet mapping function
    MFh=(1.0d0+(ah/(1.0d0+bh/(1.0d0+ch))))/(dsind(Ele)+(ah/(dsind(Ele)+bh/(dsind(Ele)+ch))))+ &
        (1.0d0/dsind(Ele)-(1.0d0+(aht/(1.0d0+bht/(1.0d0+cht))))/(dsind(Ele)+(aht/(dsind(Ele)+bht/(dsind(Ele)+cht)))))*H/1000.0d0
    MFw=(1.0d0+(aw/(1.0d0+bw/(1.0d0+cw))))/(dsind(Ele)+(aw/(dsind(Ele)+bw/(dsind(Ele)+cw))))

    return
end subroutine
    
   