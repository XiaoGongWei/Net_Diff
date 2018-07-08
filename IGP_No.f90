!********************************************************************
!  函数名： IGP_No
!
!  用途：   格网点编号与经纬度转换
!           Mode=0 由经纬度求IGP号
!           Mode=1 由IGP号求经纬度
!
!  输入、输出参数：
!           integer*4 phi, lam 格网点经纬度 角度
!           integer*4 ipgno    格网点编号
!           integer*4 mode     转换控制
!           integer*4 ier      错误标志
! 2012-1-21 4:07 将格网点经纬度改为实数，适应格网加密
!********************************************************************

      subroutine IGP_No(Phi, Lam, IGPNo,  Mode,Ier) 
      IMPLICIT none
!      include 'wadpdef.fti'
!.... I/O变量
!      integer*4 phi, lam
      real*8     phi, lam
      integer*4 IGPNo
      integer*4 IGPNosequent
      integer*4 mode
      integer*4 ier 

! .... 格网划分 行、列数     
      integer*4 N_Row 
      parameter( N_Row = 20)     
      integer*4 N_Col            
      parameter( N_Col = 16)

! .... 格网划分 起始经纬度 角度
      real*8  B0
      parameter( B0 = 7.5) 
      real*8 L0
      parameter( L0 = 70.0)

! .... 格网划分 间隔 角度      
      real*8 Delta_B      
      parameter( Delta_B = 2.5)
      real*8 Delta_L      
      parameter( Delta_L = 5.0)
! .... 格网划分 截止经纬度 角度
      integer*4 B1
      parameter( B1 = 55)
      integer*4 L1
      parameter( L1 = 145)

!.... 内部使用变量
      Integer*4 i_row,i_col, ntol

!.... 初始化错误标志      
      Ier = 0

!.... 总格网点个数      
      ntol = n_row * n_col

!.... 经纬度 --> 编号
      if (mode .eq. 0) then
          if (phi .lt. b0 .or. phi .gt. b1 .or. lam .lt. l0 .or. lam .gt. l1) then
             igpno = 0
             ier = 1
          else
             i_row = int((phi - b0) / delta_b) + 1
             i_col = int((lam - l0) / delta_l) + 1
             IGPNosequent = (i_col - 1) * n_row + i_row               
             if(mod(IGPNosequent,2) .eq. 0)then
!...igpno是偶数的(2,4,6,...,320)==>1-160
               igpno=IGPNosequent/2
             else
!...igpno是奇数的(1,3,5,...,319)==>161-320
               igpno=(IGPNosequent+1)/2+160
             endif
          end if
!.... 编号 --> 经纬度
      elseif (mode .eq. 1) then 
          if (igpno .le. 0 .or. igpno .gt. ntol) then
             phi = 0
             lam = 0
             ier=1
          else
             if(igpno .le. ntol/2)then
               IGPNosequent=igpno*2
             else
               IGPNosequent=(igpno-ntol/2)*2-1
             endif
             i_col = int((IGPNosequent -1)/ n_row ) + 1
             i_row = IGPNosequent - (i_col - 1) * n_row
             phi = b0 + (i_row - 1) * delta_b
             lam = l0 + (i_col - 1) * delta_l
          end if
!.... mode错误
      else
          ier=1
      endif

      return
      end
      