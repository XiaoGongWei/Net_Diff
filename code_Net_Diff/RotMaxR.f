Copyright (c) Tongji University 2006.
c    All rights reserved.
c    函数名：RotMaxR
c
c    功能： 求旋转矩阵
c
c    Created by Jiexian Wang April 2006
c    Modified by zStone JUNE 30 2006
c
c     input:
c         real*8 Phi, Lam       旋转角度
c
c      output:
c         real*8 R(3,3)         旋转矩阵
c         integer*4 Ier         错误标志
c********************************************************************

      subroutine RotMaxR(Phi, Lam, R,Ier)
      IMPLICIT none

c.... 输入输出变量
      real*8 Phi, Lam
      real*8 R(3,3)
      integer*4 Ier

c.... 内部变量            
      real*8 R2(3, 3), R3(3, 3)
      
      Ier =0

c.... y轴旋转phi      
      call  RotMat1(2, Phi, R2,Ier)
      if(Ier .ne. 0) goto 999
      
c.... z轴旋转-lam      
      call RotMat1(3, -Lam, R3,Ier) 
      if(Ier .ne. 0) goto 999 
      
c.... 旋转矩阵      
      call  MultiplyRR(R, R3, R2,Ier)
      if(Ier .ne. 0) goto 999
      
 999  return
      end
      
Copyright (c) Tongji University 2006.
c    All rights reserved.
c
c    函数名： MultiplyRR
c        
c    功能：   两个3*3 矩阵相乘
c
c    Created by Jiexian Wang April 2006
c
c    Modified by zStone JULY 7 2006
c
c     input:
c        real*8 R1(3,3),R2(3, 3)   总个数 移动个数
c        integer*4 ier
c
c     output:
c        real*8 R(3, 3)
c
c****************************************************************** 
      subroutine MultiplyRR(R, R1, R2,Ier)
      IMPLICIT none
      Integer*4 I,J,K,Ier
      real*8 R(3,3),R1(3, 3), R2(3, 3)

      Ier=0

      do 10 i = 1 , 3
	    do 10 j = 1 , 3
		  R(i, j) = 0.0d0
	      do 10 k = 1 , 3
		     R(i, j) = R(i, j) + R1(i, k) * R2(k, j)
 10   continue

      return
      end

Copyright (c) Tongji University 2006.
c    All rights reserved.
c    函数名：RotMat
c
c    功能： 求旋转矩阵
c
c    Created by Jiexian Wang April 2006
c    Modified by zStone JUNE 30 2006
c
c     input:
c         integer*4 N           旋转控制
c         real*8 Alf            旋转角度 弧度
c
c      output:
c         real*8 R(3,3)         旋转矩阵
c         integer*4 Ier         错误标志
c********************************************************************

      subroutine RotMat1(N, Alf, R,Ier)
      IMPLICIT none
      integer*4 N,Ier
      real*8 R(3,3),Alf

      Ier=0

c.... 以x轴为旋转轴      
      if ( N .eq. 1) then
	    R(1, 1) = 1.0d0
	    R(1, 2) = 0.0d0
	    R(1, 3) = 0.0d0
	    R(2, 1) = 0.0d0
	    R(2, 2) = DCOS(Alf)
	    R(2, 3) = DSIN(Alf)
	    R(3, 1) = 0.0d0
	    R(3, 2) = -DSIN(Alf)
	    R(3, 3) = DCOS(Alf)

c.... 以y轴为旋转轴
      elseif ( N .eq. 2) then
	    R(1, 1) = DCOS(Alf)
	    R(1, 2) = 0.0d0
	    R(1, 3) = -DSIN(Alf)
	    R(2, 1) = 0.0d0
	    R(2, 2) = 1.0d0
	    R(2, 3) = 0.0d0
	    R(3, 1) = Sin(Alf)
	    R(3, 2) = 0.0d0
	    R(3, 3) = DCOS(Alf)

c.... 以z轴为旋转轴
      elseif ( N .eq.3) then
	    R(1, 1) = DCOS(Alf)
	    R(1, 2) = DSIN(Alf)
	    R(1, 3) = 0.0d0
	    R(2, 1) = -DSIN(Alf)
	    R(2, 2) = DCOS(Alf)
	    R(2, 3) = 0.0d0
	    R(3, 1) = 0.0d0
	    R(3, 2) = 0.0d0
	    R(3, 3) = 1.0d0
      else
		ier=1 
      endif
      return
      end      
      