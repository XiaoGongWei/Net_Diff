! ============= Sovle_NEQ===================== 
! PURPOSE:
!            Solve the normal equation and get the coordinate.
!
! INPUTS:
!         NEQ                         normal equation structure, used in Step 1
!
! OUTPUT:
!         Coor                        Coordinate
!
! WRITTEN BY: Yize Zhang
! =============================================

subroutine Solve_NEQ_Net(NEQ, Coor)
use MOD_FileID
use MOD_NEQ
use MOD_Epo_NEQ
use MOD_constant
implicit none
    type(type_NEQ) :: NEQ
    type(type_Epo_NEQ) :: Epo_NEQ
    integer :: i, PRN, N
    logical :: AD_flag
    real(8) :: Coor(3)
    real(8) :: maxV
    integer :: maxL

    Ad_Flag=.true.
    NEQ%maxV=0.d0
    NEQ%maxL=0
    N=NEQ%PRNS
    do while(AD_flag)
        Ad_Flag=.false.
        call Invsqrt(NEQ%Nbb, NEQ%N, NEQ%InvN)
        NEQ%dx=matmul(NEQ%InvN, NEQ%U)   ! In distance(meter)

        ! =================== Outliers Detect =====================
        NEQ%Vp1(1:N)=matmul(NEQ%Ap1(1:N, 1:5), NEQ%dx(1:5)) - NEQ%Lp1(1:N)
        NEQ%maxV(1:1)=maxval(dabs(NEQ%Vp1(1:N)))
        NEQ%maxL(1:1)=maxloc(dabs(NEQ%Vp1(1:N)))
        NEQ%Vp2(1:N)=matmul(NEQ%Ap2(1:N, 1:5), NEQ%dx(1:5)) - NEQ%Lp2(1:N)
        NEQ%maxV(2:2)=maxval(dabs(NEQ%Vp2(1:N)))
        NEQ%maxL(2:2)=maxloc(dabs(NEQ%Vp2(1:N)))
        NEQ%Vwl(1:N)=matmul(NEQ%Awl(1:N, :), NEQ%dx) - NEQ%Lwl(1:N)
        NEQ%maxV(3:3)=maxval(dabs(NEQ%Vwl(1:N)))
        NEQ%maxL(3:3)=maxloc(dabs(NEQ%Vwl(1:N)))
        NEQ%Vw4(1:N)=matmul(NEQ%Aw4(1:N, :), NEQ%dx) - NEQ%Lw4(1:N)
        NEQ%maxV(4:4)=maxval(dabs(NEQ%Vw4(1:N)))
        NEQ%maxL(4:4)=maxloc(dabs(NEQ%Vw4(1:N)))

        maxV=maxval(dabs(NEQ%maxV))
        maxL=maxloc(dabs(NEQ%maxV),dim=1)
        
        if ( (dabs(maxV)>2.d0)  ) then
            if ((maxL)==1) then   ! maxV in P1
                call Minus_NEQ( NEQ%Nbb(1:5,1:5), NEQ%U(1:5), NEQ%Ap1(1:N,:), NEQ%Lp1(1:N), &
                       NEQ%P(1:N, 1:N), 5,  N, NEQ%maxL(1), NEQ%SumN)
            elseif  ((maxL)==2) then   ! maxV in P2
                call Minus_NEQ( NEQ%Nbb(1:5,1:5), NEQ%U(1:5), NEQ%Ap2(1:N,:), NEQ%Lp2(1:N), &
                       NEQ%P(1:N, 1:N), 5,  N, NEQ%maxL(2), NEQ%SumN)
            elseif  ((maxL)==3) then
                call Minus_NEQ( NEQ%Nbb, NEQ%U, NEQ%Awl(1:N,:), NEQ%Lwl(1:N), &
                       NEQ%P(1:N, 1:N), NEQ%N,  N, NEQ%maxL(3), NEQ%SumN)
            elseif  ((maxL)==4) then
                call Minus_NEQ( NEQ%Nbb, NEQ%U, NEQ%Aw4(1:N,:), NEQ%Lw4(1:N), &
                       NEQ%P(1:N, 1:N), NEQ%N,  N, NEQ%maxL(4), NEQ%SumN)
            end if
!            if (NEQ%SumN<3+N*2) exit
            Ad_Flag=.true.
        end if 
        ! =================== End of Outliers Detect =====================
    end do

    Coor=NEQ%dx(1:3)
    return
    
end subroutine