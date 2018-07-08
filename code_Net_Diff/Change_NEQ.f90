! ===================  Change_NEQ  ================
! 
! PURPOSE:
!         Change the Normal Equation by adding or minus a observation.
!    
! INPUTS:
!           Nbb, U               The normal equation  Nbb*x=U, N is a low triangular matrix
!           N                        The size of the normal equation
!           A                        The vector of the observation
!           L                         The constant part of the observation
!           operation           The operation of the obseration, "add" or "cut"
!
! OUTPUT:
!           Nbb, U                 The new normal equation.
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ================= End of header  =================


subroutine Change_NEQ(Nbb, U, N, A, L,operation)
implicit none
    ! Intent in and out
    real(8) :: Nbb(N,N), U(N), A(N), L
    integer :: N
    character(3) :: operation
    ! Local Variables
    integer :: i, j

    if (index(operation,"add") /= 0) then
        do i=1,N
            if (A(i)==0.d0) cycle
            do j=i,N
                if (A(j)==0.d0) cycle
                Nbb(j,i)=Nbb(j,i)+A(i)*A(j)
            end do
            U(i)=U(i)+A(i)*L
        end do
    else if (index(operation,"cut") /= 0) then
        do i=1,N
            if (A(i)==0.d0) cycle
            do j=i,N
                if (A(j)==0.d0) cycle
                Nbb(j,i)=Nbb(j,i)-A(i)*A(j)
            end do
            U(i)=U(i)-A(i)*L
        end do
    else
        write(*,*) "operation method not found: ",operation
        pause
        stop
    end if

end