! 
! ===================Minus_NEQ================
! Purpose:
!                 Change the Normal Equation by minus an error observation.
!    
! Inputs: 
!               Nbb, U               The normal equation  Nbb*x=U
!               A                        The matrix of the error equations
!               L                         The constant part of the error equations
!               N                        The size of the normal equation
!               M                       The number of the error equations
!               MaxL                  The location of error equation to be deleted
!               Sum                    The  total number of error equations
!
!  Outputs:
!              Nbb, U                 The new normal equation.
!              Sum                      The new total number of error equations
!
! Written by: Yize Zhang
! =================End of header=================

subroutine Minus_NEQ( Nbb, U, A, L, P, N, M, maxL, Sum)
implicit none
    real(8) :: Nbb(N,N), U(N)
    real(8) :: A(M, N), L(M), P(M, M)
    integer :: N, M
    integer :: maxL, Sum

    Nbb=Nbb - matmul(  matmul( transpose(A), P ), A  )
    U=U - matmul(  matmul( transpose(A), P ), L  )
    A( maxL, :)=0.d0
    L( maxL )=0.d0
    Nbb=Nbb + matmul(  matmul( transpose(A), P ), A  )
    U=U + matmul(  matmul( transpose(A), P ), L  )
    Sum=Sum-1
end subroutine

