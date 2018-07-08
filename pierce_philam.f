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
c        real*8 XSTA, YSTA, ZSTA     测站 地固系直角坐标
c        real*8 XSAT, YSAT, ZSAT     卫星 地固系直角坐标
c
c     output:
c        real*8 Phi, Lam             穿刺点球面坐标
c        real*8 RcIono               倾斜因子
c        real*8 Elevation            高度角
c        real*8 Azimuth              方位角
c        Integer*4 Ier               错误标志        
c
c****************************************************************** 
      subroutine pierce_philam(Xsta, Ysta , Zsta
     . , Xsat, Ysat, Zsat ,Phi, Lam, RcIono, Elevation, Azimuth,Ier) 

      IMPLICIT none
      
      !include 'wadpdef.fti'	
c.... 输入、输出变量
      real*8 XSTA, YSTA, ZSTA
      real*8 XSAT, YSAT, ZSAT
      real*8 Phi, Lam
      real*8 RcIono
      real*8 Elevation
      real*8 Azimuth
      Integer*4 Ier
      real(8) Earth_Radius 
      real(8) Iono_Hei
           
c.... 内部使用变量
      real*8 DIST
      real*8 EX
      real*8 EY
      real*8 EZ
      real*8 B
      real*8 C
      real*8 X,Y,Z
      EARTH_RADIUS=6378.137d0
      Iono_Hei = 375.d0
c.... 初始化错误标志
      ier = 0
c.... 计算站星距离
      Dist = (Xsat - Xsta)**2 + (Ysat - Ysta)**2 + (Zsat - Zsta)**2
      Dist = dsqrt(Dist)

c.... 站星方向矢量
      ex = (Xsat - Xsta) / Dist
      ey = (Ysat - Ysta) / Dist
      ez = (Zsat - Zsta) / Dist

c.... 求解IPP位置
c      print *, 'Iono_Hei', Iono_Hei, 'Earth_Radius', Earth_Radius
      B = 2.0d0 * (Xsta * ex + Ysta * ey + Zsta * ez)
      c = Xsta * Xsta + Ysta * Ysta + Zsta * Zsta 
     . - (Earth_Radius  + Iono_Hei)*(Earth_Radius + Iono_Hei)
     . *1000000.0d0
      Dist = -B + dsqrt(B * B - 4.0d0 * c)
      Dist = Dist / 2.0d0

c.... IPP三维直角坐标
      X = Xsta + Dist * ex
      Y = Ysta + Dist * ey
      Z = Zsta + Dist * ez

c.... 化为球面坐标
      Phi = datan(Z / dsqrt(X * X + Y * Y))
      Lam = datan2(Y, X)

c.... 倾斜因子
      RcIono = ex * X + ey * Y + ez * Z
      RcIono = RcIono / dsqrt(X * X + Y * Y + Z * Z)

c.... 计算高度角、方位角    
 	  call  XYZUEN(Xsta, Ysta, Zsta, Xsat, Ysat, Zsat
     . , Elevation, Azimuth,Ier) 
 

      return
      end