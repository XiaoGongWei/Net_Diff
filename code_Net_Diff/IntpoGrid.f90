
      subroutine IntpoGrid(Phi , Lam , VDelay ,IGP_Delay0 , IGP_Delay_SiGma0,Ier)
      IMPLICIT none
!      include 'wadpdef.fti'
!     zstone 2006/07/07
!     input:
!            real*8 phi					IPP坐标
!            real*8 lam
!            real*8 IGP_Delay0			电离层格网计算结果
!            real*8 IGP_Delay_SiGma0	电离层格网中误差
!
!      output:
!            real*8 VDelay				输入IPP点电离层垂直延迟
!            integer*4 Ier				错误信息
   
! 输入输出变量
      real*8 Phi
      real*8 Lam
      real*8 IGP_Delay0(320)
      real*8 IGP_Delay_SiGma0(320)
      real*8 VDelay
      Integer*4 Ier
      
! 内部使用变量   
      Integer*4 i
      Integer*4 i_Col
      Integer*4 i_Row
      Integer*4 NN(4)
      Integer*4 N
      Integer*4 Ntol
      Integer*4 I_Perice
      real*8 Xmin
      real*8 Xmax
      real*8 Ymin
      real*8 Ymax
      real*8 Ww(0:4)
      real*8 Xpp
      real*8 Ypp
      real*8 Get_UDRE_P 

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
       
        Ier = 0
                 
!      parameter( N_Row = 10
      
       Ntol = N_Row * N_Col * 2  

      
      VDelay = 0.0d0      
! 穿刺点左下角格网点的行、列数        
      i_Row = Int((Phi - B0) / Delta_B) + 1
      i_Col = Int((Lam - L0) / Delta_L) + 1

! 穿刺点与格子左下角的差   Xpp Ypp 范围 0～1       
      Xpp = (Phi - Int(Phi / Delta_B) * Delta_B) / Delta_B
      Ypp = (Lam - Int(Lam / Delta_L) * Delta_L) / Delta_L

! zstone  四个格网点编号  4(左上)    1(右上)
!                         3(左下)    2(右下)
      call IGP_No((i_Row + 0) * Delta_B + B0, (i_Col + 0) * Delta_L + L0, NN(1), 0,Ier)  
      if(Ier .ne. 0) print *, 'IntpoGrid: IGP_No error'
      call IGP_No((i_Row - 1) * Delta_B + B0, (i_Col + 0)  * Delta_L + L0, NN(2), 0,Ier)  
      if(Ier .ne. 0) print *, 'IntpoGrid: IGP_No error'
      call IGP_No((i_Row - 1) * Delta_B + B0, (i_Col - 1)      * Delta_L + L0, NN(3), 0,Ier)  
      if(Ier .ne. 0) print *, 'IntpoGrid: IGP_No error'
      call IGP_No((i_Row + 0) * Delta_B + B0, (i_Col - 1)   * Delta_L + L0, NN(4), 0,Ier) 
      if(Ier .ne. 0) print *, 'IntpoGrid: IGP_No error'

      N = 0
      do i = 1 , 4
        If (NN(i) .gt. 0 ) N = N + 1
      enddo
!      print *, 'IntpoGrid .f  NN',NN
!      print *,'IntpoGrid .f  N' ,N
      If (N .lt. 3) GoTo 999
            
      Ww(1) = Xpp * Ypp
      Ww(2) = (1.0d0 - Xpp) * Ypp
      Ww(3) = (1.0d0 - Xpp) * (1.0d0 - Ypp)
      Ww(4) = Xpp * (1.0d0 - Ypp)


      VDelay = 0.0d0
      Ww(0) = 0.0d0
      N = 0
      do i = 1 , 4
        If ((IGP_Delay_SiGma0(NN(i)) .gt. 0.0d0) .and. (IGP_Delay_SiGma0(NN(i)) .lt. 5.0d0)) Then
          N = N + 1
          VDelay = VDelay + IGP_Delay0(NN(i)) * Ww(i)
          Ww(0) = Ww(0) + Ww(i)
        End If
      enddo

      If (N .ge. 3) Then
        VDelay = VDelay / Ww(0)
      Else
      	VDelay = 0.d0
!       print *,'IntpoGrid: N < 3 error'
        GoTo 999
      End If
      

!      Ww(0) = 0.0d0      
!      N = 0
!      do i = 1 , 4
!        If (IGP_Delay_SiGma0(NN(i)+ Ntol/2) .gt. 0.0d0) Then
!          N = N + 1
!          VDelay_phs =VDelay_phs +
!     .     IGP_Delay0(NN(i)+Ntol/2) * Ww(i)
!          Ww(0) = Ww(0) + Ww(i)
!        End If
!      enddo
! 
!      If (N .ge. 3) Then
!        VDelay_phs = VDelay_phs / Ww(0)
!      Else
! c        print *,'IntpoGrid: N < 3 error'
!        VDelay_phs = 0.0d0
!      End If
      
      return

 999  ier=1
      return
      end
