! ==================== InvSqrt ==================
! 
! PURPOSE:
!             Inverse of a symmetric positive matrix, usually used in a normal function.
!
! METHOD: 
!             Square root Factorization(Cholesky Factorization)
!
! NOTICE:
!             If all of the elements in row i and line i are zeros, matrix A is rank defect.
!         Here we set A(i,i)=1.d0 to make A undefect, so in the inverse of A (InvA), 
!         all the elements in line i and low i are zeros except InvA(i,i)=1.d0.
! 
! INPUTS:
!     A             a symmetric positive matrix, type: real(8).
!     n             the size of A is n*n, type: integer.
! 
! OUTPUTS:
!     InvA        Inverse matrix of A, type: real(8).
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
!==================  End of Header  =====================
 
subroutine InvSqrt(AA,n,InvA)
implicit none
    ! Intent in
    integer :: n
    real(8):: A(n,n)
    real(8) :: AA(n,n)
    ! Intent out
    real(8) :: InvA(n,n)
    ! Local variable
    integer :: i, j, s
    real(8) :: GCho(n,n), C(n,n), gg, g, gy(n), gx(n)
    real(8) :: b(n)
    
    Gcho=0.0d0
    A=AA
    InvA=0.d0
    ! Cholesky Factorization
    do i=1,n
!        if  ( (all(dabs(A(i:n,i))<1.d-6)) .and. (all(dabs(A(i,1:i))<1.d-6)) ) then
        if  ( (all(dabs(A(i:n,i))<1.d-9)) .and. (all(dabs(A(i,1:i))<1.d-9)) ) then
            A(i,i)=1.d0    ! Check that if matrix A is rank defect.
        end if
        do j=i,n
            if (i==j) then
                g=0.d0
                do s=1, (i-1)
                    g=g+GCho(i,s)*GCho(i,s)
                end do
                gg=dsqrt(A(j,i) - g)
                GCho(j,i)=gg
            else
                g=0.d0
                do s=1,(i-1)
                    g=g+GCho(i,s)*GCho(j,s)
                end do
                GCho(j,i)=(A(j,i) - g)/gg
                !GCho(i,j)=GCho(j,i)
            end if
        end do
    end do
    
    ! GCho*y=I, get y
    do i=1,n
        gy=0.d0
        b=0.d0
        b(i)=1.d0
        do j=1,(i-1)
            gy=gy+GCho(i,j)*C(:,j)
        end do
        C(:,i)=(b-gy)/GCho(i,i)
    end do

    ! GCho'*InvA=y, get InvA
    do i=n,1,-1
        gx=0.d0
        do j=(i+1),n
            gx=gx+GCho(j,i)*InvA(:,j)
        end do
        InvA(:,i)=(C(:,i)-gx)/GCho(i,i)
    end do

    do i=1,n
        if  ( (all(dabs(AA(i:n,i))<1.d-11)) .and. (all(dabs(AA(i,1:i))<1.d-11)) ) then
            InvA(i,i)=0.d0    ! If A is rank defect, restroe InvA
        end if
    end do

    return
end subroutine