Copyright (c) Tongji University 2006.
c    All rights reserved.
c    函数名：XYZUEN
c
c    功能： 计算高度角、方位角
c
c    Created by Jiexian Wang April 2006
c    Modified by zStone JUNE 30 2006
c
c     input:
c        real*8 Xsta, Ysta, Zsta         测站坐标
c        real*8 Xsat, Ysat, Zsat         卫星坐标
c
c     output:
c        real*8 Elevation, Azimuth       站星高度角 方位角
c        integer*4 ier                   错误标志
c********************************************************************

      subroutine XYZUEN(Xsta, Ysta, Zsta, Xsat, Ysat, Zsat
     . , Elevation, Azimuth,ier)
      IMPLICIT none
      !include 'wadpdef.fti'	
      	   
c     输入输出参数
      real*8 Xsta, Ysta, Zsta
      real*8 Xsat, Ysat, Zsat
      real*8 Elevation, Azimuth
      integer*4 ier
      real(8) :: pi

c     内部变量     
      real*8 Phi, Lam,XYZ(3), UEN(3), R(3, 3), H
      Integer*4 I,J

      pi=3.14159265358d0
      ier=0

c.... 站星方向
      XYZ(1) = Xsat - Xsta
      XYZ(2) = Ysat - Ysta
      XYZ(3) = Zsat - Zsta

c.... 测站大地坐标
      !call XYZ2BLH(Phi, Lam, H, Xsta, Ysta, Zsta,Ier)
      call XYZ2BLH(Xsta, Ysta, Zsta,Phi, Lam, H)  !这个函数是以度数为单位，需要转成弧度
      Phi=Phi*pi/180.d0
      Lam=Lam*pi/180.d0
      if (Ier .ne. 0) then 
         print *, 'Error XYZUEN.f call XYZ2BLH error ', phi,' ', lam
      endif
      
c.... 地心系 转为 站心系
      call  RotMaxR(Phi, Lam, R,Ier)
      if (Ier .ne. 0) then 
         print *, 'Error XYZUEN.f call RotMaxR error', R  
      endif

c.... 计算站星站心地平坐标
      do 10 i = 1 , 3
         UEN(i) = 0.0d0
         do 10 j = 1 , 3
            UEN(i) = UEN(i) + R(j, i) * XYZ(j)
 10   continue
 
c.... 高度角
      Elevation = datan(UEN(1)/dsqrt(UEN(2)*UEN(2) + UEN(3)*UEN(3)))

c.... 方位角
      Azimuth = datan2(UEN(2), UEN(3))
      if (Azimuth .lt.0.0d0) Azimuth=Azimuth+pi*2.0d0
      
 999  continue
      return
	end