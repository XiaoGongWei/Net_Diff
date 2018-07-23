! ============= Reset_NEQ===================== 
! PURPOSE:
!            This subroutine is to reset NEQ if reference satellite change.
!

! INPUTS:
!         DD                            double difference structure
!         NEQ                         normal equation structure
!         Epo_NEQ                 current epoch normal equation structure
! OUTPUT:
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! =================End of header=================

subroutine Reset_NEQ(RefSat, DD, NEQ, Epo_NEQ)
use MOD_DD
use MOD_NEQ
use MOD_Epo_NEQ
use MOD_constant
use MOD_Var
use MOD_FileID
implicit none
    type(type_DD) :: DD
    type(type_NEQ) :: NEQ
    type(type_Epo_NEQ) :: Epo_NEQ
    integer :: RefSat(5)
    ! Local variables
    integer :: i, j, sys

    j=0
    do i=ParaNum+1,NEQ%N+IonoNum
        if (i-ParaNum>2*MaxPRN) then
            j=2*MaxPRN
        elseif (i-ParaNum>MaxPRN) then
            j=MaxPRN
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
            if (i<=MaxPRN+ParaNum) then
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
            elseif (i<=2*MaxPRN+ParaNum) then
                if ((i/=RefSat(sys)+MaxPRN+ParaNum) .and. (NEQ%dx(i)/=0.d0)) then   ! For L2
                    NEQ%dx(i)=NEQ%dx(i)-NEQ%dx(RefSat(sys)+MaxPRN+ParaNum)  ! other satellite
                    NEQ%InvN(:,i)=NEQ%InvN(:,i)-NEQ%InvN(:,RefSat(sys)+MaxPRN+ParaNum)
                    NEQ%InvN(i,:)=NEQ%InvN(i,:)-NEQ%InvN(RefSat(sys)+MaxPRN+ParaNum,:)
                    if (ar_mode==3) then ! If fixed and hold mode
                        if ((NEQ%fixed_amb(RefSat(sys)+MaxPRN)/=0.99d0) .and. (NEQ%fixed_amb(i-ParaNum)/=0.99d0)) then
                            NEQ%fixed_amb(i-ParaNum)=NEQ%fixed_amb(i-ParaNum)-NEQ%fixed_amb(RefSat(sys)+MaxPRN)
                        else
                            NEQ%fixed_amb(i-ParaNum)=0.99d0
                        end if
                    end if
                end if
                
                if ((i/=RefSat(sys)+MaxPRN+ParaNum) .and. If_Est_Iono .and. IonoNum>0 ) then
                    if  (Epo_NEQ%dx(i)==0.d0) then
                        cycle
                    end if
                    Epo_NEQ%dx(i)=Epo_NEQ%dx(i)-Epo_NEQ%dx(RefSat(sys)+MaxPRN+ParaNum)  ! other satellite
                    Epo_NEQ%InvN(:,i)=Epo_NEQ%InvN(:,i)-Epo_NEQ%InvN(:,RefSat(sys)+MaxPRN+ParaNum)
                    Epo_NEQ%InvN(i,:)=Epo_NEQ%InvN(i,:)-Epo_NEQ%InvN(RefSat(sys)+MaxPRN+ParaNum,:)
                    if (ar_mode==3) then ! If fixed and hold mode
                        if ((Epo_NEQ%fixed_amb(RefSat(sys)+MaxPRN)/=0.99d0) .and. (Epo_NEQ%fixed_amb(i-ParaNum)/=0.99d0)) then
                            Epo_NEQ%fixed_amb(i-ParaNum)=Epo_NEQ%fixed_amb(i-ParaNum)-Epo_NEQ%fixed_amb(RefSat(sys)+MaxPRN)
                        else
                            Epo_NEQ%fixed_amb(i-ParaNum)=0.99d0
                        end if
                    end if
                end if
            else
                if ((i/=RefSat(sys)+2*MaxPRN+ParaNum) .and. (Epo_NEQ%dx(i)/=0.d0)) then   ! For ionosphere parameter
                    Epo_NEQ%dx(i)=Epo_NEQ%dx(i)-Epo_NEQ%dx(RefSat(sys)+MaxPRN*2+ParaNum)  ! other satellite
                    Epo_NEQ%InvN(:,i)=Epo_NEQ%InvN(:,i)-Epo_NEQ%InvN(:,RefSat(sys)+MaxPRN*2+ParaNum)
                    Epo_NEQ%InvN(i,:)=Epo_NEQ%InvN(i,:)-Epo_NEQ%InvN(RefSat(sys)+MaxPRN*2+ParaNum,:)
                end if
            end if
        end if
    end do
    do sys=1,5
        if ( (RefSat(sys) /= DD%RefSat(sys))  .and. (DD%RefSat(sys)/=0) .and. (RefSat(sys)/=0) ) then
            write(LogID,'(I6,2X,I2,A3,I2,A15)') sys, DD%RefSat(sys),'-->',RefSat(sys),'ref sat change'  
            NEQ%dx(DD%RefSat(sys)+ParaNum)=0.d0-NEQ%dx(RefSat(sys)+ParaNum)  ! old reference satellite
            NEQ%dx(DD%RefSat(sys)+MaxPRN+ParaNum)=0.d0-NEQ%dx(RefSat(sys)+MaxPRN+ParaNum)
            NEQ%InvN(:,DD%RefSat(sys)+ParaNum)= 0.d0 -NEQ%InvN(:,RefSat(sys)+ParaNum)
            NEQ%InvN(:,DD%RefSat(sys)+MaxPRN+ParaNum)= 0.d0 -NEQ%InvN(:,RefSat(sys)+MaxPRN+ParaNum)
            NEQ%InvN(DD%RefSat(sys)+ParaNum,:)= 0.d0 -NEQ%InvN(RefSat(sys)+ParaNum,:)
            NEQ%InvN(DD%RefSat(sys)+MaxPRN+ParaNum,:)= 0.d0 -NEQ%InvN(RefSat(sys)+MaxPRN+ParaNum,:)
            NEQ%dx(RefSat(sys)+ParaNum)=0.d0   ! new reference satellite
            NEQ%dx(RefSat(sys)+MaxPRN+ParaNum)=0.d0
            NEQ%InvN(:,RefSat(sys)+ParaNum)=0.d0
            NEQ%InvN(RefSat(sys)+ParaNum,:)=0.d0
            NEQ%InvN(:,RefSat(sys)+MaxPRN+ParaNum)=0.d0
            NEQ%InvN(RefSat(sys)+MaxPRN+ParaNum,:)=0.d0
                    
            if (If_Est_Iono .and. IonoNum>0) then
                Epo_NEQ%dx(DD%RefSat(sys)+ParaNum)=0.d0-Epo_NEQ%dx(RefSat(sys)+ParaNum)  ! old reference satellite
                Epo_NEQ%dx(DD%RefSat(sys)+MaxPRN+ParaNum)=0.d0-Epo_NEQ%dx(RefSat(sys)+MaxPRN+ParaNum)
                Epo_NEQ%dx(DD%RefSat(sys)+2*MaxPRN+ParaNum)=0.d0-Epo_NEQ%dx(RefSat(sys)+2*MaxPRN+ParaNum)
                Epo_NEQ%InvN(:,DD%RefSat(sys)+ParaNum)= 0.d0 -Epo_NEQ%InvN(:,RefSat(sys)+ParaNum)
                Epo_NEQ%InvN(:,DD%RefSat(sys)+MaxPRN+ParaNum)= 0.d0 -Epo_NEQ%InvN(:,RefSat(sys)+MaxPRN+ParaNum)
                Epo_NEQ%InvN(:,DD%RefSat(sys)+2*MaxPRN+ParaNum)= 0.d0 -Epo_NEQ%InvN(:,RefSat(sys)+2*MaxPRN+ParaNum)
                Epo_NEQ%InvN(DD%RefSat(sys)+ParaNum,:)= 0.d0 -Epo_NEQ%InvN(RefSat(sys)+ParaNum,:)
                Epo_NEQ%InvN(DD%RefSat(sys)+MaxPRN+ParaNum,:)= 0.d0 -Epo_NEQ%InvN(RefSat(sys)+MaxPRN+ParaNum,:)
                Epo_NEQ%InvN(DD%RefSat(sys)+2*MaxPRN+ParaNum,:)= 0.d0 -Epo_NEQ%InvN(RefSat(sys)+2*MaxPRN+ParaNum,:)
                Epo_NEQ%dx(RefSat(sys)+ParaNum)=0.d0   ! new reference satellite
                Epo_NEQ%dx(RefSat(sys)+MaxPRN+ParaNum)=0.d0
                Epo_NEQ%dx(RefSat(sys)+2*MaxPRN+ParaNum)=0.d0
                Epo_NEQ%InvN(:,RefSat(sys)+ParaNum)=0.d0
                Epo_NEQ%InvN(RefSat(sys)+ParaNum,:)=0.d0
                Epo_NEQ%InvN(:,RefSat(sys)+MaxPRN+ParaNum)=0.d0
                Epo_NEQ%InvN(RefSat(sys)+MaxPRN+ParaNum,:)=0.d0
                Epo_NEQ%InvN(:,RefSat(sys)+2*MaxPRN+ParaNum)=0.d0
                Epo_NEQ%InvN(RefSat(sys)+2*MaxPRN+ParaNum,:)=0.d0
            end if

            if (ar_mode==3) then ! If fixed and hold mode
                if (NEQ%fixed_amb(RefSat(sys))/=0.99d0) then
                    NEQ%fixed_amb(DD%RefSat(sys))=0.d0-NEQ%fixed_amb(RefSat(sys))    ! new fixed ambiguity of old reference satellite
                else
                    NEQ%fixed_amb(DD%RefSat(sys))=0.99d0
                end if
                if (NEQ%fixed_amb(RefSat(sys)+MaxPRN)/=0.99d0) then
                    NEQ%fixed_amb(DD%RefSat(sys)+MaxPRN)=0.d0-NEQ%fixed_amb(RefSat(sys)+MaxPRN)
                else
                    NEQ%fixed_amb(DD%RefSat(sys)+MaxPRN)=0.99d0
                end if
                NEQ%fixed_amb(RefSat(sys))=0.99d0    ! new fixed ambiguity of new reference satellite
                NEQ%fixed_amb(RefSat(sys)+MaxPRN)=0.99d0

                if (If_Est_Iono .and. IonoNum>0) then
                    if (Epo_NEQ%fixed_amb(RefSat(sys))/=0.99d0) then
                        Epo_NEQ%fixed_amb(DD%RefSat(sys))=0.d0-Epo_NEQ%fixed_amb(RefSat(sys))    ! new fixed ambiguity of old reference satellite
                    else
                        Epo_NEQ%fixed_amb(DD%RefSat(sys))=0.99d0
                    end if
                    if (Epo_NEQ%fixed_amb(RefSat(sys)+MaxPRN)/=0.99d0) then
                        Epo_NEQ%fixed_amb(DD%RefSat(sys)+MaxPRN)=0.d0-Epo_NEQ%fixed_amb(RefSat(sys)+MaxPRN)
                    else
                        Epo_NEQ%fixed_amb(DD%RefSat(sys)+MaxPRN)=0.99d0
                    end if
                    Epo_NEQ%fixed_amb(RefSat(sys))=0.99d0    ! new fixed ambiguity of new reference satellite
                    Epo_NEQ%fixed_amb(RefSat(sys)+MaxPRN)=0.99d0
                end if
            end if

            call InvSqrt(NEQ%InvN, NEQ%N, NEQ%Nbb)
            NEQ%U=MATMUL(NEQ%Nbb, NEQ%dx)   ! New Nbb and U is needed

            if (If_Est_Iono .and. IonoNum>0) then
                call InvSqrt(Epo_NEQ%InvN, Epo_NEQ%N, Epo_NEQ%Nbb)
                Epo_NEQ%U=MATMUL(Epo_NEQ%Nbb, Epo_NEQ%dx)   ! New Nbb and U is needed
            end if
        end if
    end do

    return
end subroutine