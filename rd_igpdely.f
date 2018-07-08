      module Pcorfile_RD
      implicit none
            
c.... 读取的格网点电离层的最大历元数 86400/180*1天 =480
      integer*4 MAX_GRID_RDEPOCH
      parameter ( MAX_GRID_RDEPOCH = 500)
      
c.... 电离层格网点个数
      integer*4 IONOGRIDPOINT_NUM
      parameter (IONOGRIDPOINT_NUM = 320)      
 
      integer*4,save :: neph_igp     

c ... 格网改正数
      type igpdata
      real*8 Igpdely
      real*8 Igpdely_sigma
      end type igpdata

      type read_igp
      integer*1 IorQ
      integer*2 WN
      integer*4 Sow
      integer*4 igpcount
      type(igpdata) :: igp(IONOGRIDPOINT_NUM)
      end type read_igp      
      
c..   统计读取的格网点电离层延迟   
      type(read_igp),save  :: igp_all(MAX_GRID_RDEPOCH) 

      end module Pcorfile_RD


! ********************************************************************
!  rd_igpdely_online
!  write by : cyl 2011-6-4 13:31
! ********************************************************************

      subroutine rd_igpdely_online(yr,mn,dy) 
      use Pcorfile_RD
      use MOD_FileDir
      use MOD_FileID
      implicit none

      integer*4 yr,mn,dy
      integer*4 ier, ierr
      integer*4 Idigpfile,Idtime,inum,i,j,k,ilast
      character*611 :: rdline
      character*60 :: igpfn
      integer*2    :: wk
      integer*4    :: isow,stid,dt,sowrd,wkrd
      integer*4    :: ii(16),wklast,isowlast,ineph
      real*8       :: ipp(16),ippsig(16)
      integer*4    :: week0,isow0,iiwk,iisw
      logical      :: firstime
      real*8,allocatable       :: igprd(:),igpsigmrd(:)  
      
      allocate(igprd(IONOGRIDPOINT_NUM),STAT=ier)
      allocate(igpsigmrd(IONOGRIDPOINT_NUM),STAT=ier)
           
          
      neph_igp = 0
      inum = 0
      firstime = .true.
      
      do i =1,MAX_GRID_RDEPOCH
      	igp_all(i).igpcount = 0
      enddo
      	
      	
      
      igprd = 0.d0
      igpsigmrd = -1.d0
      
      Idigpfile = FileID_Mark
      FileID_Mark = FileID_Mark + 1

      write(igpfn,'(a,a,i4.4,i2.2,i2.2,a4)') 
     . trim(IonDir),'GIVEI@',
     .    yr,mn,dy, '.dat'
      
      open(unit = Idigpfile,file = trim(igpfn),
     .     status = 'old',iostat = ierr)
      if(ierr .ne. 0 ) then
      	print *,'read GIVE@yyyymmdd online Err!'
      pause
      	stop
      endif  
      
       

 5    continue
   
      
      do i = 1,20
      	
 6       read(Idigpfile,'(a611)',err=99, end= 100) rdline

         read(rdline,'(i4,i8,3x,16(i6,f8.3,8X,f8.3,7x))',err=99) 
     .       wkrd,sowrd,(ii(j),ipp(j),ippsig(j),j=1,16)
     
         
             
         do k = 1,16
         	  igprd(ii(k)) = ipp(k)
         	  igpsigmrd(ii(k)) = ippsig(k)
         enddo

      enddo

      
8     continue

      if( firstime ) then
      	week0 = wkrd
      	isow0 = sowrd   
      	firstime = .false.
      endif
      

      neph_igp = neph_igp + 1
      
      if(neph_igp .eq. 1) then
      	
      	igp_all(neph_igp).IorQ = 2
      	igp_all(neph_igp).WN   = wkrd
      	igp_all(neph_igp).Sow  = sowrd
      	
      	do i = 1,160
      	 igp_all(neph_igp).igp(i).Igpdely = igprd(i)
      	 igp_all(neph_igp).igp(i).Igpdely_sigma = igpsigmrd(i)
      	enddo  
      	    	
      else
      	
      	ilast = neph_igp-1
      	wklast = igp_all(ilast).WN       
      	isowlast = igp_all(ilast).Sow      	
      	
      	
      	dt = (wkrd -wklast)*604800 + (sowrd-isowlast)
      	
      	if( dt .lt. 0 ) then
      		goto 5
      	endif
      		
      	if (dt .ne. 180) then
      		
      		print *,'last time',wklast,isowlast
      		print *, 'time now',wkrd,sowrd
      		inum = int(dt / 180)
c 如果时间相差不足180秒，按180秒保存数组  
c 如果时间相差超过180秒，将中间数组补齐，使数组等间隔存储，间隔18秒    		
      	 if (inum .le. 1 ) then
      	      igp_all(neph_igp).IorQ = 2
      	      igp_all(neph_igp).WN   = wklast
      	      igp_all(neph_igp).Sow  = isowlast + 180
      	      
      	      do i = 1,320
      	       igp_all(neph_igp).igp(i).Igpdely = igprd(i)
      	       igp_all(neph_igp).igp(i).Igpdely_sigma = igpsigmrd(i)
      	      enddo
      	        
          else
          	
          	do i = 1,inum -1
          	  ineph = neph_igp+i-1
          	  
      	      igp_all(ineph).IorQ = 2
      	      igp_all(ineph).WN   = wklast
      	      igp_all(ineph).Sow  = isowlast + 180*i
      	      
      	      do k = 1,320
      	      	
      	       igp_all(ineph).igp(k).Igpdely = 
     .         igp_all(ilast).igp(k).Igpdely
      	       igp_all(ineph).igp(k).Igpdely_sigma = 
     .         igp_all(ilast).igp(k).Igpdely_sigma
     
      	      enddo
      	      
          	enddo 
          	
          	neph_igp = neph_igp + inum -1
          	
      	    igp_all(neph_igp).IorQ = 2
      	    igp_all(neph_igp).WN   = wklast
      	    igp_all(neph_igp).Sow  = isowlast + 180*(inum)
      	    
      	    do i = 1,160
      	     igp_all(neph_igp).igp(i).Igpdely = igprd(i)
      	     igp_all(neph_igp).igp(i).Igpdely_sigma = igpsigmrd(i)
      	    enddo         	
          	
          endif  ! inum =1
          
        else   ! dt ==180
        	  
      	      igp_all(neph_igp).IorQ = 2
      	      igp_all(neph_igp).WN   = wkrd
      	      igp_all(neph_igp).Sow  = sowrd
      	      
      	      do i = 1,320
      	       igp_all(neph_igp).igp(i).Igpdely = igprd(i)
      	       igp_all(neph_igp).igp(i).Igpdely_sigma = igpsigmrd(i)
      	      enddo
          
      	endif	! dt ~/==180
      	
      	
      endif ! neph_igp ~=1      
      
      goto 5
      
         
99    print *,  rdline
      goto 6    

 
      
100   continue  
  
      close(Idigpfile) 

      
c      do i = 1, neph_igp
c      	 do j = 1,160
c      	 	igprd(j) = igp_all(i).igp(j).Igpdely 
c      	 enddo
c        write(155,'(I4,2X,I6)') igp_all(i).WN,igp_all(i).Sow
c        write(155,'(10(f8.4,1x))') igprd(1:160)     
c      enddo

      deallocate(igprd)
      deallocate(igpsigmrd) 
           
      return
      end subroutine



      
!********************************************************************
!  get_igpdely
!  write by : cyl 2011-6-4 13:31
!********************************************************************
           
      subroutine get_igpdely(tinweek,tinsow,
     .           IGP_Delay,IGP_Delay_SiGma,ierr)
      use Pcorfile_RD
      IMPLICIT none  
!      include './input/wadpdef.fti'
      
      real*8,parameter :: hsec=180.d0 
      
      integer*4,intent(in) :: tinweek
      real*8,intent(in) :: tinsow
      integer*4,intent(out) :: ierr  
      REAL*8,intent(out) :: IGP_Delay(IONOGRIDPOINT_NUM)
      real*8,intent(out) :: IGP_Delay_SiGma(IONOGRIDPOINT_NUM) 
      
      integer*4 :: ndata,idt,ii,ieph
      integer*4 :: idtfirst,idtend
      integer*4 :: firstwk,firstsow,endwk,endsow      
      
      ierr = -1
      
      firstwk = igp_all(1).WN
      firstsow = igp_all(1).sow
      
      endwk = igp_all(neph_igp).WN
      endsow = igp_all(neph_igp).sow
      
     
      idtfirst =int((tinweek - firstwk)*604800 + 
     .  int(tinsow- firstsow) )    
      
      if(idtfirst .lt. 0 ) then
      	 print *, 'input time before first time in get_igpdely!'
      	 print *, 'input time',tinweek,tinsow
      	 print *, 'first time',firstwk,firstsow
      	 ierr = -1
      	 goto 124
      endif
      
      idtend = int(( tinweek-endwk)*604800 +
     .  int(tinsow-endsow ) )
      
      if(idtend .gt. 200) then
      	 print *, 'input time after last time in get_igpdely!'
      	 print *, 'input time',tinweek,tinsow
      	 print *, 'end time',endwk,endsow
      	 ierr = -1
      	 goto 124
      endif
      
      ieph =  int(idtfirst/180) + 1
      
      if( ieph .eq. 0) then  
      	 print *, 'input time befor first time in get_pcor!'
      	 ierr = -1  
      	 goto 124
      	
      elseif (ieph .gt. 0 .and. ieph .le. neph_igp) then
      	
       
        do ii = 1,IONOGRIDPOINT_NUM
           IGP_Delay(ii) = igp_all(ieph).igp(ii).Igpdely  
           IGP_Delay_SiGma(ii) = igp_all(ieph).igp(ii).Igpdely_sigma
        enddo   	

        
        idt =int((tinweek - igp_all(ieph).WN)*604800 + 
     .  int(tinsow - igp_all(ieph).Sow) + 0.001)
        
        if(idt .lt. 0 .or. idt .ge. 180) then
        	print *,'time err in time =',tinweek,tinsow
          print *, igp_all(ieph).WN,igp_all(ieph).Sow
          write(*,'(10(f8.4,2x))') IGP_Delay(1:160)
          write(*,'(10(f8.4,2x))') IGP_Delay_SiGma(1:160)
          stop
        endif
        

      else
      	 print *, 'input time after last time in get_pcor!'
      	 ierr = -1  
      	 goto 124
      endif
      	
      ierr = 0	
     
      
 124  continue   
      
      
      return
      end subroutine