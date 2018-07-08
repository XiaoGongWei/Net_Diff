! ================ IonoCompensate ==============
! Function: 
!              This subroutine is to Compensate the residuals
!  in double difference(DD) OMC, which mainly comes
!  from ionosohere difference in two stations
!
! INPUTS:
!       Week           ObsWeek
!       Sow             ObsSow
!      RefSat          ObsData type
!       DD              DD struct
!      Epo_NEQ     Epo_NEQ struct
!
! WRITTEN BY: Yize Zhang, zhyize@163.com
! ================ End of Header ===============

subroutine IonoCompensate(Week, Sow,RefSat, DD,Epo_NEQ)
use MOD_DD
use MOD_Epo_NEQ
use MOD_FileID
use MOD_IonoDDRes
use MOD_VAR
implicit none
    type(type_DD) :: DD
    type(type_Epo_NEQ) :: Epo_NEQ
    integer :: Week
    real(8) :: Sow
    integer :: RefSat(5)
    integer(2) :: i, sys, PRN
    real(8) :: dL1, dL2

    ! 1. Change DD residuals if reference satelite change
    do PRN=1, MaxPRN
        if (PRN<=GNum) then
            sys=1
        elseif (PRN<=GNum+RNum) then
            sys=2
        elseif (PRN<=GNum+RNum+CNum) then
            sys=3
        elseif (PRN<=GNum+RNum+CNum+NumE) then
            sys=4
        elseif (PRN<=GNum+RNum+CNum+NumE+JNum) then
            sys=1
        elseif (PRN<=GNum+RNum+CNum+NumE+JNum+INum) then
            sys=5
        end if
        if (RefSat(sys) /= Last_RefSat(sys) .and. Last_RefSat(sys)/=0 .and. RefSat(sys)/=0) then
            if (PRN/=RefSat(sys)) then
                do i=1,5 !Ubound(IonoDDRes(PRN)%L1)
                    if (IonoDDRes(PRN)%L1(i)/=0.d0 .and. IonoDDRes(RefSat(sys))%L1(i)/=0.d0 ) then
                        IonoDDRes(PRN)%L1(i)=IonoDDRes(PRN)%L1(i)-IonoDDRes(RefSat(sys))%L1(i)
                        IonoDDRes(PRN)%amb1(i)=IonoDDRes(PRN)%amb1(i)-IonoDDRes(RefSat(sys))%amb1(i)
                    else
                        IonoDDRes(PRN)%L1(i)=0.d0
                        IonoDDRes(PRN)%amb1(i)=0.d0
                    end if
                    if (IonoDDRes(PRN)%L2(i)/=0.d0 .and. IonoDDRes(RefSat(sys))%L2(i)/=0.d0 ) then
                        IonoDDRes(PRN)%L2(i)=IonoDDRes(PRN)%L2(i)-IonoDDRes(RefSat(sys))%L2(i)
                        IonoDDRes(PRN)%amb2(i)=IonoDDRes(PRN)%amb2(i)-IonoDDRes(RefSat(sys))%amb2(i)
                    else
                        IonoDDRes(PRN)%L2(i)=0.d0
                        IonoDDRes(PRN)%amb2(i)=0.d0
                    end if
                end do
            end if
        end if
    end do
    do sys=1,5
        if (RefSat(sys) /= Last_RefSat(sys) .and. Last_RefSat(sys)/=0 .and. RefSat(sys)/=0) then
            IonoDDRes(RefSat(sys))%L1=0.d0
            IonoDDRes(RefSat(sys))%L2=0.d0
            IonoDDRes(RefSat(sys))%amb1=0.d0
            IonoDDRes(RefSat(sys))%amb2=0.d0
            Last_RefSat(sys)=RefSat(sys)
        end if
    end do    

    ! 2. Add the rediduals to carrier phase obsevatioin by linear model
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'dL1'
    do i=1,DD%PRNS
        PRN=DD%PRN(i)
        if (DD%L1(i)/=0.d0) then
            call Linear(5, IonoDDRes(PRN)%Sow, IonoDDRes(PRN)%L1, Sow, IonoDDRes(PRN)%dL1)
            if (any(IonoDDRes(PRN)%amb1-IonoDDRes(PRN)%amb1(1)/=0.d0)) then
                IonoDDRes(PRN)%dL1=0.d0   ! If fixed ambiguity not stable
            end if
            DD%L1(i)=DD%L1(i)+IonoDDRes(PRN)%dL1
!            DD%WL(i)=DD%WL(i)+IonoDDRes(PRN)%dL1
            write(unit=LogID,fmt='(I4,F10.3)',advance='no') PRN, IonoDDRes(PRN)%dL1
        end if
    end do
    write(LogID,'(A)') ''
    write(unit=LogID,fmt='(5X,A5)',advance='no') 'dL2'
    do i=1,DD%PRNS
        PRN=DD%PRN(i)
        if (DD%L2(i)/=0.d0) then
            call Linear(5,  IonoDDRes(PRN)%Sow, IonoDDRes(PRN)%L2, Sow, IonoDDRes(PRN)%dL2)
            if (any(IonoDDRes(PRN)%amb2-IonoDDRes(PRN)%amb2(1)/=0.d0)) then
                IonoDDRes(PRN)%dL2=0.d0   ! If fixed ambiguity not stable
            end if
            DD%L2(i)=DD%L2(i)+IonoDDRes(PRN)%dL2
!            DD%W4(i)=DD%W4(i)+IonoDDRes(PRN)%dL2
            write(unit=LogID,fmt='(I4,F10.3)',advance='no') PRN, IonoDDRes(PRN)%dL2
        end if
    end do
    write(LogID,'(A)') ''
    do i=1,DD%PRNS
        PRN=DD%PRN(i)
        if ( (a1*f1+a2*f2/=0.d0) .and. (b1*f1+b2*f2/=0.d0) ) then  ! If two frequency combination
            if (DD%L1(i)/=0.d0 .and. DD%L2(i)/=0.d0) then
                DD%WL(i)=DD%WL(i)+(a1*f1*IonoDDRes(PRN)%dL1+a2*f2*IonoDDRes(PRN)%dL2)/(a1*f1+a2*f2)
            end if
        elseif ( (DD%L1(i)/=0.d0) .and. (a1/=0.d0) ) then  ! If L1 frequency combination
            DD%WL(i)=DD%WL(i)+IonoDDRes(PRN)%dL1
        elseif ( (DD%L2(i)/=0.d0) .and. (a2/=0.d0) ) then  ! If L2 frequency combination
            DD%WL(i)=DD%WL(i)+IonoDDRes(PRN)%dL2     
        end if
        if ((a1*f1+a2*f2/=0.d0) .and. (b1*f1+b2*f2/=0.d0) ) then  ! If two frequency combination
            if (DD%L1(i)/=0.d0 .and. DD%L1(i)/=0.d0) then
                DD%W4(i)=DD%W4(i)+(b1*f1*IonoDDRes(PRN)%dL1+b2*f2*IonoDDRes(PRN)%dL2)/(b1*f1+b2*f2)
            end if
        elseif ( (DD%L1(i)/=0.d0) .and. (b1/=0.d0) ) then  ! If L1 frequency combination
            DD%W4(i)=DD%W4(i)+IonoDDRes(PRN)%dL1
        elseif ( (DD%L2(i)/=0.d0) .and. (b2/=0.d0) ) then  ! If L2 frequency combination
            DD%W4(i)=DD%W4(i)+IonoDDRes(PRN)%dL2
        end if
    end do

end subroutine


! ================== Linear ===============
! PURPOSE:
!    Linear fitting
!
! INPUTS:
!   N        Nunber of points
!   X        X series
!   Y        Y series
!   Xp      Interpolated X
!
! OUTPUTS:
!   Yp       Interpolated Y
!
! WRITTEN BY: Yize Zhang, zhyize@163.com
! ================ End of Header ============

subroutine Linear(N, X,Y,Xp,Yp)
implicit none
    integer :: N
    real(8) :: X(N), Y(N), Xp
    real(8)  :: Yp
    integer :: i, K, Num
    real(8) :: A(N,2), L(N), Nbb(2,2), InvN(2,2), U(2), dx(2), V(N)
    real(8) :: sigma, maxV, SumY
    integer(1) :: maxL
    logical :: AD_Flag

    K=0
    Yp=0.d0

    SumY=0.d0
    do i=1,N
        if (Y(i)==0.d0) cycle
        if (Xp-X(i)>300.d0) cycle ! Maximum time difference is within 15min
        K=K+1
        SumY=SumY+Y(i)
    end do
    if (K<5) return
    Yp=SumY/K   ! Constant model
    return

    do i=1,N
        if (Y(i)==0.d0) cycle
        if (Xp-X(i)>900.d0) cycle ! Maximum time difference is within 15min
        K=K+1
        A(K,:)=(/1.d0, X(i)-Xp/)
        L(K)=Y(i)
    end do
    if (K<5) return

    Num=K
    AD_Flag=.true.
    do while(AD_Flag)
        AD_Flag=.false.
        Nbb=MATMUL(TRANSPOSE(A(1:K,:)),A(1:K,:))
        U=MATMUL(TRANSPOSE(A(1:K,:)),L(1:K))
        call InvSqrt(Nbb, 2, InvN)
        dx=MATMUL(InvN, U)
        V(1:K)=MATMUL(A(1:K,:),dx)-L(1:K)
        sigma=dsqrt(DOT_PRODUCT(V(1:K), V(1:K))/(Num-2) )
        maxV=maxval(dabs(V(1:K)))
        maxL=maxloc(dabs(V(1:K)),dim=1)
        if ((dabs(maxV)>3.d0*sigma)  ) then
            A(maxL,:)=0.d0
            L(maxL)=0.d0
            Num=Num-1
            if (Num<5) return
!            if (K-Num>4) return   ! If more than 4 epoches were deleted
            AD_Flag=.true.
        end if
    end do
    Yp=dx(1) !+dx(2)*Xp+dx(3)*Xp**2

end subroutine