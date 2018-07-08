c  利用14参数模型计算L1上电离层延迟iondely
cccc by gxq 2014/1/5 15:57
c*****************************************

       subroutine cal_iondelay(ISow,RStaPos,RSatPos,iono_14para
     . ,Rionvdelay,I1outerr)
      implicit none
     
      integer(kind=1)  :: I1outerr
      integer(kind=4)  :: Week0,Ierr,i,ISow
      real(kind=8)     :: RStaPos(3),RSatPos(3)      
      real(kind=8)     :: RPhi,RLam,RcIono,RElevation,RAzimuth
      real(kind=8)     :: RA1,RB,iono_14para(14),Rionvdelay
      real(kind=8)     :: Ralpha(4),Rbeta(4),Rgama(4)     
                    
      I1outerr=1   
c************************************************  
c*********电离层14参*****************************     
!      iono_14para(1)=8.750940e+00
!      iono_14para(2)=3.243312e+00
!      iono_14para(3)=1.416835e-08   
!      iono_14para(4)=1.942698e-07
!      iono_14para(5)=-1.837431e-06
!      iono_14para(6)=3.579759e-06
!      iono_14para(7)=1.274286e+05   
!      iono_14para(8)=-5.908882e+05
!      iono_14para(9)=5.955631e+06
!      iono_14para(10)=-1.138271e+07
!      iono_14para(11)=-3.561169e+03   
!      iono_14para(12)=1.652139e+05
!      iono_14para(13)=-1.257419e+06
!      iono_14para(14)=2.655165e+06                       
c**************************************************       
      RA1 = iono_14para(1)
      RB  = iono_14para(2)
      
      do i=1,4
        Ralpha(i)= iono_14para(i+2)
        Rbeta(i) = iono_14para(i+6)
        Rgama(i) = iono_14para(i+10)
      enddo
      
      call Pierce_PHiLam(RStaPos(1),RStaPos(2),RStaPos(3),
     .RSatPos(1),RSatPos(2),RSatPos(3), RPhi, RLam, RcIono, 
     .RElevation, RAzimuth,Ierr)
      if(Ierr .ne. 0) then
         print *, 'integerity_cor.f Pierce_PHiLam error'
         return
      endif
 
      call iono_14(RA1, RB, Ralpha, Rbeta, Rgama, 
     .                   RPhi, RLam, dble(ISow), Rionvdelay)
            if(Rionvdelay .lt. 0)Rionvdelay=1.5   !2.25改为1.5m 9TECu对应B1频点的电离层延迟
        ! 将天顶延迟改为斜路径延迟   基于北斗一期B1频点，单位是：m
         
       Rionvdelay = Rionvdelay/RcIono   
      I1outerr=0      
      return
      end
      
c####################################################################c
      subroutine iono_14(A1, B, alpha, beta, gama,
     .                   lat, lon, sow, iono_delay)
c---------------------------------------------------------------------+
c说明：信息处理系统电离层十四参数模型                                 |
c                                                                     |
c输入：十四参数模型 夜间平场        A1          纳秒                  |
c      十四参数模型 夜间斜率        B           纳秒/pi               |
c      十四参数模型 白天振幅        alpha(4)    秒/pi^0~3             |
c      十四参数模型 余弦周期        beta(4)     秒/pi^0~3             |
c      十四参数模型 余弦初相        gama(4)     秒/pi^0~3             |
c      穿刺点纬度                   lat         弧度                  |
c      穿刺点经度                   lon         弧度                  |
c      北斗时周秒                   sow         秒                    |
c                                                                     |
c输出：穿刺点B1频点天顶延迟         iono_delay  米                    |
c                                                                     |
c日期：2007-9-21                                                      |
c---------------------------------------------------------------------+
      implicit none
     
      real*8 A1, B, alpha(4), beta(4), gama(4)
      real*8 lat, lon, sow
      real*8 iono_delay
      real*8 A2, A3, A4, t, PI, VLIGHT
      integer*4 i
c.... Frequency L1 L2
      real*8,parameter ::  Freq1  = 1561.098d6
      real*8,parameter ::  Freq2  = 1207.14d6
      real*8,parameter ::  Freq3  = 1268.52d6
     
      data PI/3.1415926535897932384626433832795d0/
      data VLIGHT/299792458.0d0/
      
c.... 单位转换      
	  lat = lat / PI

c.... 穿刺点地方时
      t = sow + (lon / PI) * 43200.d0
      t = mod(t, 86400.d0)

c.... 白天电离层延迟余弦曲线振幅
      A2= alpha(1)

c.... 白天电离层延迟余弦曲线初相位
c      A3= gama(1)   根据新接口更改为 gama(1)+ 50400.d0  cyl&wxl 2010-6-17 1:34 
      A3= gama(1) + 50400.d0

c.... 白天电离层延迟余弦曲线周期
      A4 = beta(1)
      
      do i=1,3
        A2 = A2 + alpha(i+1)*(lat**i)
        A3 = A3 + gama(i+1)*(lat**i)
        A4 = A4 + beta(i+1)*(lat**i)
      enddo      
      
c.... 白天电离层延迟余弦曲线振幅取值判断
      if(A2 .lt. 0.d0) A2 = 0.d0

c.... 白天电离层延迟余弦曲线初相位取值判断
      if(A3 .gt. 55800.d0) then
        A3 = 55800.d0
      else 
        if(A3 .lt. 43200.d0) then
           A3 = 43200.d0
        endif
      endif

c.... 白天电离层延迟余弦曲线周期取值判断
c      if(A4 .gt. 158400.d0) then
c        A4 = 158400.d0
c      else 
c        if(A4 .lt. 86400.d0) then
c           A4 = 86400.d0
c        endif                
c      endif
cc.... 全天电离层延迟余弦曲线周期取值判断   
c 根据新的接口，改为 72000 到172800 s (20h - 48h) cyl&wxl 2010-6-17 1:34
   
      if(A4 .gt. 172800.d0) then
        A4 = 172800.d0
      else   
        if(A4 .lt. 72000.d0) then
           A4 = 72000.d0
      endif        
        
      endif
c.... 把当前时刻归算到以初相位为中心的一整天内
      do while(t .lt. (A3 - 43200.d0))
        t = t + 86400.d0
      end do
      
c  去除判断条件 与用户保持一致 cyl&wxl 2010-6-17 1:34
c      do while(t .gt. (A3 + 43200.d0))   cyl&wxl
c        t = t - 86400.d0
c      end do

c.... 计算天顶延迟
      iono_delay = (A1 - B * lat) * VLIGHT * 1e-9
c..... 如果夜间电离层值计算结果为负值，根据经验则置为1.5m
      if(0.1d0 > iono_delay) then
      iono_delay =1.5
      endif
      if( (A4 - 4.d0 * dabs(t - A3)) .gt. 0.d0) then
        iono_delay=iono_delay + A2*dcos(2.d0*PI*(t-A3)/A4)*VLIGHT
      endif
      
!!   根据2011-4-27 7:47 接口，电离层参数修改为基于B1频点（原来为基于B3频点）
!! 根据14参计算的电离层延迟为B1频点的，在格网点电离层延迟基于B3频点
!! cyl 2011-4-27 7:47
!       iono_delay =iono_delay*Freq1*Freq1/(Freq3*Freq3)
       
            
      end        