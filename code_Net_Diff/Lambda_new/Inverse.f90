
subroutine Inverse(n, A, IA)

implicit none
integer :: n
real(8) :: A(n,n), IA(n,n)
real(8) :: B(n,n)
integer :: i, j

forall(i=1:N, j=1:N, i==j) IA=1.d0
forall(i=1:N, j=1:N, i/=j) IA=0.d0
B=A

call Upper(B,IA,N)
call Lower(B,IA,N)

forall(i=1:N) IA(i,:)=IA(i,:)/B(i,i)

return

return

end subroutine

subroutine Upper(M, S, N)
implicit none
integer :: N
real(8) :: M(n,n), S(n,n)
integer :: i, j
real(8) :: E

do i=1, N-1
    do j=i+1,N
        E=M(j,i)/M(i,i)
        M(j,i:N)=M(j,i:N)-M(i,i:N)*E
        S(j,:)=S(j,:)-S(i,:)*E
    end do
end do

return
end subroutine


subroutine Lower(M, S, N)
implicit none
integer :: N
real(8) :: M(n,n), S(n,n)
integer :: i, j
real(8) :: E

do i=N, 2, -1
    do j=i-1, 1, -1
        E=M(j,i)/M(i,i)
        M(j,1:N)=M(j,1:N)-M(i,1:N)*E
        S(j,:)=S(j,:)-S(i,:)*E
    end do
end do

return
end subroutine
