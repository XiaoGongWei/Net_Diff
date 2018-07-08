! ===================Elimi_Para================
! 
! PURPOSE:
!           Eliminate parameter in the NEQ(Normal EQuation) 
!   
! INPUTS:
!            Nbb, U               The normal equation  Nbb*x=U, N is a low triangular matrix
!            N                        The size of the normal equation
!            t                          Eliminate the parameter in the order of t
!
! OUTPUTS:
!             Nbb, U                 The new normal equation.
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! =================End of header=================

subroutine Elimi_Para(Nbb,U,N,t)
implicit none
    ! Intent in and out
    real(8) :: Nbb(N,N), U(N)
    integer :: N, t
    ! Local variables
    real(8) :: temp_Nbb(N),temp_coe,temp_U,coe
    integer :: j, k
    
    if (Nbb(t,t)==0.d0) return
    temp_Nbb(1:t)=Nbb(t,1:t)
    temp_Nbb(t:N)=Nbb(t:N,t)
    temp_coe=temp_Nbb(t)
    temp_U=U(t)
    do j=1,N
        coe=temp_Nbb(j)/temp_coe
        if (coe==0.d0) cycle
        do k=j,N
            if (temp_Nbb(k)==0.d0) cycle
            Nbb(k,j)=Nbb(k,j)- temp_Nbb(k)*coe 
            if (dabs(Nbb(k,j))<1.e-11) Nbb(k,j)=0.d0
        end do
        U(j)=U(j)-temp_U*coe
    end do
    return
end subroutine