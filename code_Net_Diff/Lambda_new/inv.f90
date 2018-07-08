
subroutine inv(A,invA,N)
!  Purpose   :  计算逆矩阵
!-----------------------------------------------------
!  Input  parameters  :
!       1. A: 需要求逆的矩阵
!       2. N: 矩阵的维度
!  Output parameters  :
!       1. invA: A的逆矩阵
implicit real*8(a-z)
integer::n
integer::i
real*8::A(n,n),invA(n,n),E(n,n)
E=0.d0
!设置E为单位矩阵
do i=1,n
   E(i,i)=1.d0
end do
call mateq(A,E,invA,N,N)
end subroutine inv


subroutine mateq(A,B,X,N,M)
!  Purpose   :  高斯列主元消去法计算矩阵方程
!                 AX=B
implicit real*8(a-z)
integer::N,M,i
real*8::A(N,N),B(N,M),X(N,M)
real*8::btemp(N),xtemp(N)
do i=1,M
    btemp=B(:,i)
    call elgauss(A,btemp,xtemp,N)
    X(:,i)=xtemp
end do
end subroutine mateq


subroutine elgauss(A,b,x,N)
!  Purpose   :  高斯列主元消去法
!                 Ax=b
implicit real*8(a-z)
integer::i,k,N
integer::id_max  !主元素标号
real*8::A(N,N),b(N),x(N)
real*8::Aup(N,N),bup(N)
!Ab为增广矩阵  [Ab]
real*8::Ab(N,N+1)
real*8::vtemp1(N+1),vtemp2(N+1)
Ab(1:N,1:N)=A
Ab(:,N+1)=b
do k=1,N-1
    elmax=dabs(Ab(k,k))
    id_max=k
    do i=k+1,n
      if (dabs(Ab(i,k))>elmax) then
         elmax=dabs(Ab(i,k))
         id_max=i
      end if          
    end do
    vtemp1=Ab(k,:)
    vtemp2=Ab(id_max,:)
    Ab(k,:)=vtemp2
    Ab(id_max,:)=vtemp1   
!#########################################################
   do i=k+1,N
     temp=Ab(i,k)/Ab(k,k)
     Ab(i,:)=Ab(i,:)-temp*Ab(k,:)
   end do
end do
Aup(:,:)=Ab(1:N,1:N)
bup(:)=Ab(:,N+1)
call uptri(Aup,bup,x,n)
end subroutine elgauss


subroutine uptri(A,b,x,N)
!  Purpose   :  上三角方程组的回带方法
!                 Ax=b
implicit real*8(a-z)
integer::i,j,N
real*8::A(N,N),b(N),x(N)
x(N)=b(N)/A(N,N)
!回带部分
do i=n-1,1,-1
    x(i)=b(i)
   do j=i+1,N
    x(i)=x(i)-a(i,j)*x(j)
   end do
    x(i)=x(i)/A(i,i)
end do
end subroutine uptri