!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!                                                                               
! This is a package of subroutines to read CCSM T85 and T42 L26 data in NETCDF format   
! and to prepare Initial and boundary conditions for RegCM3. Both Global and Window of
! CCSM data are acceptable.                        
! Written By Moetasim Ashfaq Dec-2005 @ PURDUE.EDU
!                                                                               
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!                                                                               
!     SUBROUTINE CAM85                                                          
!     Read unpacked CCSM NETCDF T85 L26 (six hourly) data and save into data
!	arrays. Each varaible is read in seperate monthly data files.  
!     
!     SUBROUTINE CAM42                                          
!     Read unpacked CCSM NETCDF T42 L26 (six hourly) data and save into data 
!     arrays. Each varaible is read in seperate yearly data files.    
!     
!     SUBROUTINE HEADER_CAM85 & HEADER_CAM42				                        
!     Define pressure levels, Ak and Bk coeffcients and global grid dimensions
!     In CCSM, the vertical coordinate is a hybrid sigma-pressure system
!     The pressure is defined as:
!     P=Ak*Po+Bk*PS(i,j)
!     All 3D fields required for ICBC are at mid-points, so
!     Ak refers to hyam, the hybrid A coefficient at midpoints, and
!     Bk refers to hybm, the hybrid B coefficient at midpoints
!     Po=1000mb  
!                                                  				
!     SUBROUTINE GET_CAM85	                                                
!     Main subroutine to read data arrays from CAM85 and prepare ICBCs at RCM Grid 
!     
!     SUBROUTINE GET_CAM42	                                                
!     Main subroutine to read data arrays from CAM42 and prepare ICBCs at RCM Grid 
!     
!     SUBROUTINE INITDATE3			                                        
!     Initialize CCSM 365 days calendar (No leap years)
!	                        
!     SUBROUTINE CAMCLNDR                                                           
!     Subroutine for SST preparation with No leap year
!     
!     SUBROUTINE HANDLE_ERR                                                        
!     Handle error for NETCDF calls                                             
!                                                                               
!	SUBROUTINES CALLED FROM ICBC.f                                                 
!     1) INTLIN 	
!     2) INTLOG               
!     3) HUMIF1FV    
!     3) BILINX2CR 	
!     4) BILINX2DT 	
!     5) TOP2BTM 					
!     6) UVROT4 						
!     7) INTGTB 
!     8) P1P2                                                         
!     8) INTPSN 	
!     9) INTV3 
!     10) MKSST 
!     10) HUMID2FV 								
!     10) INTV1 	
!     11) WRITEF 							        
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!     Common Blocks								
!     CAMvars for arrays containing data on GCM grid and sigma levels 	
!     CAMvars1 for arrays containing data on GCM grid and pressure levels
!     CAMvars2 for arrays containing data on RCM grid and pressure levels  
!     CAMvars4 for arrays containing data on RCM grid and sigma levels 	
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!     DATA PREPARATION
!     Dataset required to use this code can be preapred using NCO utilitiles such as
!     NCKS, NCRCAT etc.
!     Prepare:
!     Monthly data files for CAM85 
!     Yearly data files for CAM42  
!     For example:
!     To extract global data of CAM42 for Specific Humidity
!     ncks -v Q input.nc ccsm.shum.nyear.nc		,and	
!     to extract a subset (window) of CAM85 data for Specific Humidity
!     ncks -d lat,min,max -d lon,min,max -v Q input.nc ccsm.shumJAN.nyear.nc
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!     NAMING CONVENTION (Global Data Files) CAM85 
!     (MONTH) = JAN/FEB/MAR/APR/MAY/JUN/JUL/AUG/SEP/OCT/NOV/DEC
!     
!     ccsm.air(MONTH).nyear.nc	 	for 'T'     (Temperature)			        
!     ccsm.hgt(MONTH).nyear.nc 		for 'Z3'    (Geopotential Height)		        
!     ccsm.shum(MONTH).nyear.nc  	for 'Q'     (Specific Humidity)		
!     ccsm.uwnd(MONTH).nyear.nc		for 'U'     (Zonal Wind)			
!     ccsm.vwnd(MONTH).nyear.nc		for 'V'     (Meridonial Wind)			
!     ccsm.pres(MONTH).nyear.nc		for 'PS'    (Surface Pressure)
!		      
!     PATH /DATA/CAM85/NYEAR/							
!     ccsm_ht.nc			for 'PHIS'  (Surface Geopotential-static field)
!		
!     PATH /DATA/CAM85/
!	
!     NAMING CONVENTION (Global Data Files) CAM42						        	
!     ccsm.air.nyear            	for 'T'    (Temperature)			        
!     ccsm.hgt.nyear.nc 		for 'Z3'   (Geopotential Height)		        
!     ccsm.shum.nyear.nc 		for 'Q'	   (Specific Humidity)		
!     ccsm.uwnd.nyear.nc		for 'U'	   (Zonal Wind)			
!     ccsm.vwnd.nyear.nc		for 'V'	   (Meridonial Wind)			
!     ccsm.pres.nyear.nc		for 'PS'	(Surface Pressure)
!		      
!     PATH /DATA/CAM42/NYEAR/			
!     ccsm_ht.nc			for 'PHIS'	(Surface Geopotential-static field)
!		
!     PATH /DATA/CAM42/
!	
!       		                                                                           
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
      subroutine CAM85(idate,idate0,xlon,xlat,glon,glat,jx,iy
     &     ,i0,i1,j0,j1)							
      implicit	none										
      include 	'netcdf.inc'					
      integer 	idate,idate0									
      
      integer	ilon,jlat,klev
      parameter (ilon=256)
      parameter (jlat=128)
      parameter (klev=26)
      
      character*3 nmonth(12)
      data nmonth/'JAN','FEB','MAR','APR','MAY','JUN','JUL'
     &     ,      'AUG','SEP','OCT','NOV','DEC'/

      real,	allocatable:: work(:,:,:)
      real,	allocatable:: work1(:)
      real,	allocatable:: work2(:)
      real	Tvar,Hvar,Qvar,Uvar,Vvar,PSvar,ZSvar
      common /CAM85vars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &     ,    	Qvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &     ,    	Vvar(ilon,jlat,klev), PSvar(ilon,jlat)
     &     ,    	ZSvar(ilon,jlat)

      character*25 inname
      character*39 pathaddname
      character*2  varname(6)
      data varname/'T','Z3','Q','U','V','PS'/
      integer 	nyear,month,nday,nhour
      integer 	inet6,ivar6,start,count
      common /CAM85open/ inet6(6),ivar6(6),start(10),count(10)
      real		checklat,checklon,checktim
      integer	starta,counta
      common /CAM85open2/ starta(10),counta(10)
      common /CAM85check/ checklat(6),checklon(6),checktim(6)
      integer 	status,icode
      integer 	kkrec,nlev,ilev,inet,ivar,it
      integer 	lonid,latid,timid,latlen,lonlen,timlen
      integer	i,j,k,ii,jj
      integer	i0,i1,j0,j1,rec
      integer	jx,iy,jmin,jmax
      real	nlon0,nlon1,nlat0,nlat1
      real    	glat(jlat), glon(ilon)
      real      xlon(jx,iy),xlat(jx,iy)
      logical 	there
      
      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-NYEAR*10000-MONTH*100
      NHOUR=IDATE-NYEAR*1000000-MONTH*10000-NDAY*100
      do kkrec=1,6
         if(kkrec.le.5) nlev=klev
         if(kkrec.eq.6) nlev=0
         if(idate.eq.idate0.or.
     &        ((mod(idate,100000).eq.10100.and.
     &        mod(idate,1000000).ne.110100)).or.
     &        mod(idate,10000).eq.100) then
            if(kkrec.eq.1) then
               write(inname,100) nyear,  'air',nmonth(month),nyear
            else if(kkrec.eq.2) then
               write(inname,100) nyear,  'hgt',nmonth(month),nyear
            else if(kkrec.eq.3) then
               write(inname,200) nyear, 'shum',nmonth(month),nyear
            else if(kkrec.eq.4) then
               write(inname,200) nyear, 'uwnd',nmonth(month),nyear
            else if(kkrec.eq.5) then
               write(inname,200) nyear, 'vwnd',nmonth(month),nyear
            else if(kkrec.eq.6) then
               write(inname,200) nyear, 'pres',nmonth(month),nyear
            endif
 100        format(I4,'/','ccsm.',A3,A3,'.',I4,'.nc')
 200        format(I4,'/','ccsm.',A4,A3,'.',I4,'.nc')
            
            pathaddname = '../DATA/CAM85/'//inname
            inquire(file=pathaddname,exist=there)
            if(.not.there) then
               print*,pathaddname,' is not available'
               stop
            endif
            
            status=nf_open(pathaddname,nf_nowrite,inet6(kkrec))
            if (status.ne.nf_noerr) call handle_err(status)
            status=nf_inq_varid(inet6(kkrec),varname(kkrec)
     &           ,ivar6(kkrec))
            if (status.ne.nf_noerr) call handle_err(status)
            write(*,*) inet6(kkrec),pathaddname,ivar6(kkrec),
     &           varname(kkrec)
         endif
         
         status=nf_inq_dimid(inet6(kkrec),'lon',lonid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimid(inet6(kkrec),'lat',latid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimid(inet6(kkrec),'time',timid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimlen(inet6(kkrec),lonid,lonlen)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimlen(inet6(kkrec),latid,latlen)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimlen(inet6(kkrec),timid,timlen)
         if (status.ne.nf_noerr) call handle_err(status)
         
         checklon(kkrec)=lonlen
         checklat(kkrec)=latlen
         checktim(kkrec)=timlen
         if (kkrec.gt.1) then
            if(checklat(kkrec).ne.checklat(kkrec-1).or. 
     &           checklon(kkrec).ne.checklon(kkrec-1).or.
     &           checktim(kkrec).ne.checktim(kkrec-1)) then
               print*, 'DOMAIN DIMENSIONS DO NOT MATCH'
               print*, 'LAT for',varname(kkrec+1),'=',checklat(kkrec)
               print*, 'LAT for',varname(kkrec),'=',checklat(kkrec)
               print*, 'LON for',varname(kkrec+1),'=',checklon(kkrec)
               print*, 'LON for',varname(kkrec),'=',checklon(kkrec)
               print*, 'TIME for',varname(kkrec+1),'=',checktim(kkrec)
               print*, 'TIME for',varname(kkrec),'=',checktim(kkrec)
               STOP    'Check CCSM Data Files'
            endif
         endif
         if(kkrec.eq.1) then
            allocate (work(lonlen,latlen,klev))
            allocate (work1(lonlen))
            allocate (work2(latlen))
         endif
         if (xlon(1,1).lt.0.0) then
            nlon0=xlon(1,1)+360.
         else
            nlon0=xlon(1,1)
         endif
         if (xlon(jx,1).lt.0.0) then
            nlon1=xlon(jx,1)+360.
         else
            nlon1=xlon(jx,1)
         endif
         nlat0=xlat(1,1)
         nlat1=xlat(1,iy)
         starta(1)=1
         counta(1)=lonlen
         status=nf_inq_varid(inet6(kkrec),'lon',lonid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_varid(inet6(kkrec),'lat',latid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_get_vara_real(inet6(kkrec),lonid,starta,counta
     &        ,work1)
         counta(1)=latlen
         status=nf_get_vara_real(inet6(kkrec),latid,starta,counta
     &        ,work2)
         
         if(nlon0.lt.work1(1).or.nlon1.gt.work1(lonlen).or.
     &        nlat0.lt.work2(1).or.nlat1.gt.work2(latlen))then
            print*, 'DOMAIN DIMENSIONS DO NOT MATCH'
            print*, 'CCSM Window LON min=',work1(1),'max=',work1(lonlen)
            print*, 'RCM  Domain LON min=',nlon0,'max=',nlon1
            print*, 'CCSM Window LAT min=',work2(1),'max=',work2(latlen)
            print*, 'RCM  Domain LAT min=',nlat0,'max=',nlat1	
            STOP    'Correct Domain Parameters in domain.param'
         endif
         
         if (idate.eq.idate0) then
            i0=work1(1)/1.40625+1
            i1=work1(lonlen)/1.40625+1
            do 300 i=1, jlat
               jmin=nint(glat(i)-work2(1))
               jmax=nint(glat(i)-work2(latlen))
               if (jmin.eq.0) j0=i
               if (jmax.eq.0) j1=i
 300        continue	   
         endif
         
         it=(NDAY-1)*4+NHOUR/6+1
         if(MONTH.eq.2) it=it+ 31*4
         if(MONTH.eq.3) it=it+ 59*4
         if(MONTH.eq.4) it=it+ 90*4
         if(MONTH.eq.5) it=it+120*4
         if(MONTH.eq.6) it=it+151*4
         if(MONTH.eq.7) it=it+181*4
         if(MONTH.eq.8) it=it+212*4
         if(MONTH.eq.9) it=it+243*4
         if(MONTH.eq.10)it=it+273*4
         if(MONTH.eq.11)it=it+304*4
         if(MONTH.eq.12)it=it+334*4
         
         count(1)=lonlen
         count(2)=latlen
         count(3)=klev
         count(4)=1
         start(1)=1
         start(2)=1
         start(3)=1
         start(4)=it
         inet=inet6(kkrec)
         ivar=ivar6(kkrec)
         if(nlev.gt.0) then
            count(3)=nlev
            status=nf_get_vara_real(inet,ivar,start,count,work)
            do ilev = 1, nlev
               if(kkrec.eq.1) then
                  do jj=j0,j1
                     j=jj-j0+1
                     if(i0.gt.i1) then
                        do ii=i0,ilon
                           i=ii-i0+1
                           Tvar(ii,jj,ilev) = work(i,j,ilev)  
                        enddo
                        do ii=1,i1
                           i=ii+(ilon-i0)+1
                           Tvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     else
                        do ii=i0,i1
                           i=ii-i0+1
                           Tvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     endif 
                  enddo
               else if(kkrec.eq.2) then
                  do jj=j0,j1
                     j=jj-j0+1
                     if(i0.gt.i1) then
                        do ii=i0,ilon
                           i=ii-i0+1
                           Hvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                        do ii=1,i1
                           i=ii+(ilon-i0)+1
                           Hvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     else
                        do ii=i0,i1
                           i=ii-i0+1
                           Hvar(ii,jj,ilev) = work(i,j,ilev)				
                        enddo
                     endif
                  enddo
               else if(kkrec.eq.3) then
                  do jj=j0,j1
                     j=jj-j0+1
                     if(i0.gt.i1) then
                        do ii=i0,ilon
                           i=ii-i0+1
                           Qvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                        do ii=1,i1
                           i=ii+(ilon-i0)+1
                           Qvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     else
                        do ii=i0,i1
                           i=ii-i0+1
                           Qvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     endif
                  enddo		
               else if(kkrec.eq.4) then
                  do jj=j0,j1
                     j=jj-j0+1
                     if(i0.gt.i1) then
                        do ii=i0,ilon
                           i=ii-i0+1
                           Uvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                        do ii=1,i1
                           i=ii+(ilon-i0)+1
                           Uvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     else
                        do ii=i0,i1
                           i=ii-i0+1
                           Uvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     endif
                  enddo	
               else if(kkrec.eq.5) then
                  do jj=j0,j1
                     j=jj-j0+1
                     if(i0.gt.i1) then
                        do ii=i0,ilon
                           i=ii-i0+1
                           Vvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                        do ii=1,i1
                           i=ii+(ilon-i0)+1
                           Vvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     else
                        do ii=i0,i1
                           i=ii-i0+1
                           Vvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     endif
                  enddo
               endif
            enddo
            
         else if(nlev.eq.0)then
            count(3)=1
            start(3)=it
            status=nf_get_vara_real(inet,ivar,start,count,work)
            if(kkrec.eq.6) then
               do jj=j0,j1
                  j=jj-j0+1
                  if(i0.gt.i1) then
                     do ii=i0,ilon
                        i=ii-i0+1
                        PSvar(ii,jj) = work(i,j,1)*0.01
                     enddo
                     do ii=1,i1
                        i=ii+(ilon-i0)+1
                        PSvar(ii,jj) = work(i,j,1)*0.01
	             enddo
                  else
                     do ii=i0,i1
                        i=ii-i0+1
                        PSvar(ii,jj) = work(i,j,1)*0.01
                     enddo
                  endif
               enddo
            endif
         endif
      enddo
      return
      end
      
      
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
      
      subroutine CAM42(idate,idate0,xlon,xlat,glon,glat,jx,iy
     &     ,i0,i1,j0,j1)
      implicit	none
      include 	'netcdf.inc'
      integer 	idate,idate0
      
      integer	ilon,jlat,klev
      parameter (ilon=128)
      parameter (jlat=64)
      parameter (klev=26)
      
      real,	allocatable:: work(:,:,:)
      real,	allocatable:: work1(:)
      real,	allocatable:: work2(:)
      real	Tvar,Hvar,Qvar,Uvar,Vvar,PSvar,ZSvar
      common /CAM42vars/ 	Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &     ,    	Qvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &     ,    	Vvar(ilon,jlat,klev), PSvar(ilon,jlat)
     &     ,    	ZSvar(ilon,jlat)
      
      character*25 inname
      character*39 pathaddname
      character*2  varname(6)
      data varname/'T','Z3','Q','U','V','PS'/
      integer 	nyear,month,nday,nhour
      integer 	inet6,ivar6,start,count
      common /CAM42open/ inet6(6),ivar6(6),start(10),count(10)
      integer	starta,counta
      common /CAM42open2/ starta(10),counta(10)
      real	checklat,checklon,checktim 
      common /CAM42check/checklat(6),checklon(6),checktim(6)
      integer 	status,icode
      integer 	kkrec,nlev,ilev,inet,ivar,it
      integer 	lonid,latid,latlen,lonlen,timid,timlen
      integer	i,j,k,ii,jj
      integer	i0,i1,j0,j1,rec
      integer	jx,iy,jmin,jmax
      real	nlon0,nlon1,nlat0,nlat1
      real    	glat(jlat), glon(ilon)
      real      xlon(jx,iy),xlat(jx,iy)
      logical 	there
	
      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-NYEAR*10000-MONTH*100
      NHOUR=IDATE-NYEAR*1000000-MONTH*10000-NDAY*100
      do kkrec=1,6
         if(kkrec.le.5) nlev=klev
         if(kkrec.eq.6) nlev=0
         if(idate.eq.idate0.or.
     &        ((mod(idate,100000).eq.10100).and.
     &        mod(idate,1000000).ne.110100)) then
            if(kkrec.eq.1) then
               write(inname,100) nyear,  'air.',nyear 
            else if(kkrec.eq.2) then
               write(inname,100) nyear,  'hgt.',nyear
            else if(kkrec.eq.3) then
               write(inname,200) nyear, 'shum.',nyear
            else if(kkrec.eq.4) then
               write(inname,200) nyear, 'uwnd.',nyear
            else if(kkrec.eq.5) then
               write(inname,200) nyear, 'vwnd.',nyear
            else if(kkrec.eq.6) then
               write(inname,200) nyear, 'pres.',nyear
            endif
 100        format(I4,'/','ccsm.',A4,I4,'.nc')
 200        format(I4,'/','ccsm.',A5,I4,'.nc')
            
            pathaddname = '../DATA/CAM42/'//inname
            inquire(file=pathaddname,exist=there)
            if(.not.there) then
               print*,pathaddname,' is not available'
               stop
            endif
            
            status=nf_open(pathaddname,nf_nowrite,inet6(kkrec))
            if (status.ne.nf_noerr) call handle_err(status)
            status=nf_inq_varid(inet6(kkrec),varname(kkrec)
     &           ,ivar6(kkrec))
            if (status.ne.nf_noerr) call handle_err(status)
            write(*,*) inet6(kkrec),pathaddname,ivar6(kkrec),
     &           varname(kkrec)
         endif
         
         status=nf_inq_dimid(inet6(kkrec),'lon',lonid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimid(inet6(kkrec),'lat',latid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimid(inet6(kkrec),'time',timid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimlen(inet6(kkrec),lonid,lonlen)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimlen(inet6(kkrec),latid,latlen)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimlen(inet6(kkrec),timid,timlen)
         if (status.ne.nf_noerr) call handle_err(status)
         
         
         checklon(kkrec)=lonlen
         checklat(kkrec)=latlen
         checktim(kkrec)=timlen
         if (kkrec.gt.1) then
            if(checklat(kkrec).ne.checklat(kkrec-1).or. 
     &           checklon(kkrec).ne.checklon(kkrec-1).or.
     &           checktim(kkrec).ne.checktim(kkrec-1)) then
               print*, 'DOMAIN DIMENSIONS DO NOT MATCH'
               print*, 'LAT for',varname(kkrec+1),'=',checklat(kkrec)
               print*, 'LAT for',varname(kkrec),'=',checklat(kkrec)
               print*, 'LON for',varname(kkrec+1),'=',checklon(kkrec)
               print*, 'LON for',varname(kkrec),'=',checklon(kkrec)
               print*, 'TIME for',varname(kkrec+1),'=',checktim(kkrec)
               print*, 'TIME for',varname(kkrec),'=',checktim(kkrec)
               STOP    'Check CCSM Data Files'
            endif
         endif
         if(kkrec.eq.1) then
            allocate (work(lonlen,latlen,klev))
            allocate (work1(lonlen))
            allocate (work2(latlen))
         endif
         if (xlon(1,1).lt.0.0) then
            nlon0=xlon(1,1)+360.
         else
            nlon0=xlon(1,1)
         endif
         if (xlon(jx,1).lt.0.0) then
            nlon1=xlon(jx,1)+360.
         else
            nlon1=xlon(jx,1)
         endif
         nlat0=xlat(1,1)
         nlat1=xlat(1,iy)
         starta(1)=1
         counta(1)=lonlen
         status=nf_inq_varid(inet6(kkrec),'lon',lonid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_varid(inet6(kkrec),'lat',latid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_get_vara_real(inet6(kkrec),lonid,starta,counta
     &        ,work1)
         counta(1)=latlen
         status=nf_get_vara_real(inet6(kkrec),latid,starta,counta
     &        ,work2)
         
         if(nlon0.lt.work1(1).or.nlon1.gt.work1(lonlen).or.
     &        nlat0.lt.work2(1).or.nlat1.gt.work2(latlen))then
            print*, 'DOMAIN DIMENSIONS DO NOT MATCH'
            print*, 'CCSM Window LON min=',work1(1),'max=',work1(lonlen)
            print*, 'RCM  Domain LON min=',nlon0,'max=',nlon1
            print*, 'CCSM Window LAT min=',work2(1),'max=',work2(latlen)
            print*, 'RCM  Domain LAT min=',nlat0,'max=',nlat1	
            STOP    'Correct Domain Parameters in domain.param'
         endif
         
         if (idate.eq.idate0) then
            i0=work1(1)/2.8125+1
            i1=work1(lonlen)/2.8125+1
            do 300 i=1, jlat
               jmin=nint(glat(i)-work2(1))
               jmax=nint(glat(i)-work2(latlen))
               if (jmin.eq.0) j0=i
               if (jmax.eq.0) j1=i
 300        continue	   
         endif
         
         it=(NDAY-1)*4+NHOUR/6+1
         if(MONTH.eq.2) it=it+ 31*4
         if(MONTH.eq.3) it=it+ 59*4
         if(MONTH.eq.4) it=it+ 90*4
         if(MONTH.eq.5) it=it+120*4
         if(MONTH.eq.6) it=it+151*4
         if(MONTH.eq.7) it=it+181*4
         if(MONTH.eq.8) it=it+212*4
         if(MONTH.eq.9) it=it+243*4
         if(MONTH.eq.10)it=it+273*4
         if(MONTH.eq.11)it=it+304*4
         if(MONTH.eq.12)it=it+334*4
         
         count(1)=lonlen
         count(2)=latlen
         count(3)=klev
         count(4)=1
         start(1)=1
         start(2)=1
         start(3)=1
         start(4)=it
         inet=inet6(kkrec)
         ivar=ivar6(kkrec)
         if(nlev.gt.0) then
            count(3)=nlev
            status=nf_get_vara_real(inet,ivar,start,count,work)
            do ilev = 1, nlev
               if(kkrec.eq.1) then
                  do jj=j0,j1
                     j=jj-j0+1
                     if(i0.gt.i1) then
                        do ii=i0,ilon
                           i=ii-i0+1
                           Tvar(ii,jj,ilev) = work(i,j,ilev)  
                        enddo
                        do ii=1,i1
                           i=ii+(ilon-i0)+1
                           Tvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     else
                        do ii=i0,i1
                           i=ii-i0+1
                           Tvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     endif 
                  enddo
				 
               else if(kkrec.eq.2) then
                  do jj=j0,j1
                     j=jj-j0+1
                     if(i0.gt.i1) then
                        do ii=i0,ilon
                           i=ii-i0+1
                           Hvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                        do ii=1,i1
                           i=ii+(ilon-i0)+1
                           Hvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     else
                        do ii=i0,i1
                           i=ii-i0+1
                           Hvar(ii,jj,ilev) = work(i,j,ilev)				
                        enddo
                     endif
                  enddo
               else if(kkrec.eq.3) then
                  do jj=j0,j1
                     j=jj-j0+1
                     if(i0.gt.i1) then
                        do ii=i0,ilon
                           i=ii-i0+1
                           Qvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                        do ii=1,i1
                           i=ii+(ilon-i0)+1
                           Qvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     else
                        do ii=i0,i1
                           i=ii-i0+1
                           Qvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     endif
                  enddo		
               else if(kkrec.eq.4) then
                  do jj=j0,j1
                     j=jj-j0+1
                     if(i0.gt.i1) then
                        do ii=i0,ilon
                           i=ii-i0+1
                           Uvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                        do ii=1,i1
                           i=ii+(ilon-i0)+1
                           Uvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     else
                        do ii=i0,i1
                           i=ii-i0+1
                           Uvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     endif
                  enddo	
               else if(kkrec.eq.5) then
                  do jj=j0,j1
                     j=jj-j0+1
                     if(i0.gt.i1) then
                        do ii=i0,ilon
                           i=ii-i0+1
                           Vvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                        do ii=1,i1
                           i=ii+(ilon-i0)+1
                           Vvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     else
                        do ii=i0,i1
                           i=ii-i0+1
                           Vvar(ii,jj,ilev) = work(i,j,ilev)
                        enddo
                     endif
                  enddo
               endif
            enddo
         else if(nlev.eq.0)then
            count(3)=1
            start(3)=it
            status=nf_get_vara_real(inet,ivar,start,count,work)
            if(kkrec.eq.6) then
               do jj=j0,j1
                  j=jj-j0+1
                  if(i0.gt.i1) then
                     do ii=i0,ilon
                        i=ii-i0+1
                        PSvar(ii,jj) = work(i,j,1)*0.01
                     enddo
                     do ii=1,i1
                        i=ii+(ilon-i0)+1
                        PSvar(ii,jj) = work(i,j,1)*0.01
	             enddo
                  else
                     do ii=i0,i1
                        i=ii-i0+1
                        PSvar(ii,jj) = work(i,j,1)*0.01
                     enddo
                  endif
               enddo
            endif
         endif
      enddo
      return
      end


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine GET_CAM42(IDATE)
      implicit  none
      integer   idate
      include   'icbc.param'
      include   'netcdf.inc'
      integer   ilon,jlat,klev
      parameter (ilon=128)
      parameter (jlat=64)
      parameter (klev=26)
		
      integer	 npl
      parameter	(npl=18)

      real    	glat,glon,sigma1,sigmar,pplev,ak,bk,psref
      common /globalCAM42/glat(jlat),glon(ilon),sigma1(npl),sigmar(npl)
     &     ,		pplev(npl),ak(klev+1),bk(klev+1),psref
      real      	Z1(ilon,jlat,klev)
      real    	Tvar,Hvar,Qvar,Uvar,Vvar,PSvar,ZSvar
      common /CAM42vars/	Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &     ,		Qvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &     ,		Vvar(ilon,jlat,klev),PSvar(ilon,jlat)
     &     ,		ZSvar(ilon,jlat)
	    	   
      real    	Tp,Qp,Hp,Up,Vp,PP3d
      
      common /CAM42var2/ 	Tp(ilon,jlat,npl),Hp(ilon,jlat,npl)
     &     ,		Qp(ilon,jlat,npl), Up(ilon,jlat,npl)
     &     ,		Vp(ilon,jlat,npl),PP3d(ilon,jlat,klev)
     
	    
      real		b2(ilon,jlat,npl*3)
      equivalence	(b2(1,1,1),Tp(1,1,1))
      real		d2(ilon,jlat,npl*2)
      equivalence	(d2(1,1,1),Up(1,1,1))
      
      real    	T3,H3,Q3,U3,V3,B3PD
      
      common /CAM42var3/ 	T3(jx,iy,npl), H3(jx,iy,npl)
     &     ,		Q3(jx,iy,npl), U3(jx,iy,npl)
     &     ,	  	V3(jx,iy,npl), B3PD(jx,iy)

      real		b3(jx,iy,npl*3)
      equivalence	(b3(1,1,1),T3(1,1,1))
      real		d3(jx,iy,npl*2)
      equivalence	(d3(1,1,1),U3(1,1,1))
      
      real    	T4,H4,Q4,U4,V4,PS4,TS4
      common /CAM42var4/ 	T4(jx,iy,kz), H4(jx,iy,kz)
     &     ,		Q4(jx,iy,kz), U4(jx,iy,kz)
     &     ,	  	V4(jx,iy,kz), PS4(jx,iy)
     &     ,	  	TS4(jx,iy)
      
          
      real    	xlon,xlat,dlon,dlat,coriol,xlandu,snowcv,topogm,toposdgm
      real	  	msfx,sigma2,sigmaf,dsigma
      common /DOMAIN/ 	xlon(jx,iy),xlat(jx,iy),dlon(jx,iy),dlat(jx,iy)
     &     ,		coriol(jx,iy),xlandu(jx,iy),snowcv(jx,iy),topogm(jx,iy)
     &     ,		toposdgm(jx,iy),msfx(jx,iy),sigma2(kz),sigmaf(kz+1)
     &     ,		dsigma(kz)
      character*6 lgtype
      real		ptop,clat,clon,delx,grdfac,plat,plon
      integer 	igrads,ibigend
      common /LGRID2/ 	ptop,clat,clon,plat,plon,delx,grdfac
     &     ,		igrads,ibigend,lgtype
	   
!     Dummy arrays for interpolation
      real	dum1(jx,iy,npl)
      real    	dum2(jx,iy,npl,2)
      real    	c1(npl), c2(kz)
		
      integer 	nyrp,nmop
      real    	wt,w,pp,h2
 
!     Arrays for new calculation of P
      real      pa(jx,iy),za(jx,iy)
      real  	tlayer(jx,iy)
      real  	sst1(jx,iy), sst2(jx,iy)
      

      character*4 varname
      character*38 pathaddname
      data varname/'PHIS'/
      real,  	allocatable:: work(:,:)
      integer   inet1,ivar1,icode,status,counta,starta
      common /CAM42open1/	starta(10),counta(10)
      integer	lonid,latid,lonlen,latlen
      integer   checklat,checklon
      integer 	i,j,k,ii,jj
      integer	i0,i1,j0,j1
      integer	imin,imax
      logical 	there
 
      call CAM42(idate,idate1,xlon,xlat,glon,glat,jx,iy,
     &     i0,i1,j0,j1)   
      if(idate.eq.idate1) then
         pathaddname= '../DATA/CAM42/ccsm_ht.nc'
         inquire(file=pathaddname,exist=there)
         if(.not.there)then
            print*,pathaddname, 'is not available'
            stop
         endif
         status=nf_open(pathaddname,nf_nowrite,inet1)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_varid(inet1,varname,ivar1)
         if (status.ne.nf_noerr) call handle_err(status)
         write(*,*)inet1,pathaddname,ivar1,varname
         status=nf_inq_dimid(inet1,'lon',lonid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimid(inet1,'lat',latid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimlen(inet1,lonid,lonlen)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimlen(inet1,latid,latlen)
         if (status.ne.nf_noerr) call handle_err(status)
         imin=i0
         imax=i1
         if (i0.gt.i1) imax=i1+ilon
         checklon=(imax-imin)+1
         checklat=(j1-j0)+1
         print*, 'i0,i1',i0,i1
         if(checklon.ne.lonlen.or.checklat.ne.latlen) then
            print*, 'DOMAIN DIMENSIONS DO NOT MATCH'
            print*, 'LAT for 3D Variables = ',checklat
            print*, 'LAT for'  ,varname,'= ',latlen
            print*, 'LON for 3D Variables = ',checklon
            print*, 'LON for'  ,varname,'= ',lonlen
            STOP    'Check Dimensions in CCSM Data Files'
         endif
         allocate (work(lonlen,latlen))
         counta(1)=lonlen
         counta(2)=latlen
         counta(3)=1
         starta(1)=1
         starta(2)=1
         starta(3)=1
         status=nf_get_vara_real(inet1,ivar1,starta,counta,work)
         if (status.ne.nf_noerr) call handle_err(status)
         do jj=j0,j1
            j=jj-j0+1
            if(i0.gt.i1) then
               do ii=i0,ilon
                  i=ii-i0+1
                  ZSvar(ii,jj)=work(i,j)/9.80616
               enddo
               do ii=1,i1
                  i=ii+(ilon-i0)
                  ZSvar(ii,jj)=work(i,j)/9.80616
               enddo
            else
               do ii=i0,i1
                  i=ii-i0+1
                  ZSvar(ii,jj)=work(i,j)/9.80616
               enddo
            endif
         enddo
      endif
      write(*,*) 'Read in fields at Date: ',idate

      do K = 1,klev
         do J = 1,jlat
            do I = 1,ilon
               if(PSvar(I,J).GT.-9995.) then
                  PP3d(I,J,K) = PSvar(I,J)*0.5*(bk(k)+bk(k+1))
     &                 +          0.5*(ak(k)+ak(k+1))
               else
                  PP3d(I,J,K) = -9999.0
               endif
            enddo
         enddo
      enddo
      
      CALL HEIGHT(Hp,Hvar,Tvar,PSvar,PP3D,ZSvar,ilon,jlat,klev,pplev
     &     ,npl)
      
      CALL INTLIN(UP,Uvar,PSvar,PP3D,ilon,jlat,klev,pplev,npl)
      CALL INTLIN(VP,Vvar,PSvar,PP3D,ilon,jlat,klev,pplev,npl)
      
      CALL INTLOG(TP,Tvar,PSvar,PP3D,ilon,jlat,klev,pplev,npl)

      CALL HUMID1FV(Tvar,Qvar,PP3D,ilon,jlat,klev)
      CALL INTLIN(QP,Qvar,PSvar,PP3d,ilon,jlat,klev,pplev,npl)
      
      CALL BILINXCX(B3,B2,xlon,xlat,glon,glat,jx,iy,ilon,jlat,npl)
      CALL BILINXDT(D3,D2,dlon,dlat,glon,glat,jx,iy,ilon,jlat,npl)
      CALL UVROT4(U3,V3,dlon,dlat,clon,clat,grdfac,jx,iy,npl
     &           ,plon,plat,lgtype)

      CALL TOP2BTM(T3,jx,iy,npl)
      CALL TOP2BTM(Q3,jx,iy,npl)
      CALL TOP2BTM(H3,jx,iy,npl)
      CALL TOP2BTM(U3,jx,iy,npl)
      CALL TOP2BTM(V3,jx,iy,npl)

      CALL INTGTB(PA,ZA,tlayer,topogm,T3,H3
     &     ,sigmar,jx,iy,npl,dum1,dum2)
 
      CALL INTPSN(PS4,topogm,PA,ZA,tlayer,ptop,jx,iy)
      CALL P1P2(B3PD,PS4,jx,iy)

      CALL INTV3(TS4,T3,PS4,sigmar,ptop,jx,iy,npl,dum1)     
      CALL CAMCLNDR( idate, nyrp, nmop, wt )
      CALL MKSST(TS4,SST1,SST2,topogm,xlandu,jx,iy,nyrp,nmop,wt)
      CALL INTV1(U4,U3,B3PD,sigma2,sigmar,ptop,jx,iy,kz
     &     ,npl,1,1,dum2,c1,c2)
      CALL INTV1(V4,V3,B3PD,sigma2,sigmar,ptop,jx,iy,kz
     &     ,npl,1,1,dum2,c1,c2)
      CALL INTV2(T4,T3,PS4,sigma2,sigmar,ptop,jx,iy,kz
     &     ,npl,1,dum2,c1,c2)
 
      CALL INTV1(Q4,Q3,PS4,sigma2,sigmar,ptop,jx,iy,kz
     &     ,npl,1,0,dum2,c1,c2)
      CALL HUMID2FV(T4,Q4,PS4,ptop,sigma2,jx,iy,kz)

      CALL HYDROST(H4,T4,topogm,PS4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)

   
      CALL WRITEF(U4,V4,T4,Q4,PP,W,PS4,TS4,ptop,jx,iy,kz,idate
     &     ,iomega)
      return
      end

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
      subroutine GET_CAM85(IDATE)
      implicit  none
      integer   idate
      include   'icbc.param'
      include   'netcdf.inc'
      integer   ilon,jlat,klev
      parameter (ilon=256)
      parameter (jlat=128)
      parameter (klev=26)
		
      integer	 npl
      parameter	(npl=18)

      real    	glat,glon,sigma1,sigmar,pplev,ak,bk,psref
      common /globalCAM85/glat(jlat),glon(ilon),sigma1(npl),sigmar(npl)
     &     ,		pplev(npl),ak(klev+1),bk(klev+1),psref
      real      	Z1(ilon,jlat,klev)
      real    	Tvar,Hvar,Qvar,Uvar,Vvar,PSvar,ZSvar
      common /CAM85vars/	Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &     ,		Qvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &     ,		Vvar(ilon,jlat,klev),PSvar(ilon,jlat)
     &     ,		ZSvar(ilon,jlat)
	    	   
      real    	Tp,Qp,Hp,Up,Vp,PP3d
      
      common /CAM85var2/ 	Tp(ilon,jlat,npl),Hp(ilon,jlat,npl)
     &     ,		Qp(ilon,jlat,npl), Up(ilon,jlat,npl)
     &     ,		Vp(ilon,jlat,npl),PP3d(ilon,jlat,klev)
     
	    
      real		b2(ilon,jlat,npl*3)
      equivalence	(b2(1,1,1),Tp(1,1,1))
      real		d2(ilon,jlat,npl*2)
      equivalence	(d2(1,1,1),Up(1,1,1))
      
      real    	T3,H3,Q3,U3,V3,B3PD
      
      common /CAM85var3/ 	T3(jx,iy,npl), H3(jx,iy,npl)
     &     ,		Q3(jx,iy,npl), U3(jx,iy,npl)
     &     ,	  	V3(jx,iy,npl), B3PD(jx,iy)

      real		b3(jx,iy,npl*3)
      equivalence	(b3(1,1,1),T3(1,1,1))
      real		d3(jx,iy,npl*2)
      equivalence	(d3(1,1,1),U3(1,1,1))
      
      real    	T4,H4,Q4,U4,V4,PS4,TS4
      common /CAM85var4/ 	T4(jx,iy,kz), H4(jx,iy,kz)
     &     ,		Q4(jx,iy,kz), U4(jx,iy,kz)
     &     ,	  	V4(jx,iy,kz), PS4(jx,iy)
     &     ,	  	TS4(jx,iy)
      
          
      real    	xlon,xlat,dlon,dlat,coriol,xlandu,snowcv,topogm,toposdgm
      real	msfx,sigma2,sigmaf,dsigma
      common /DOMAIN/ 	xlon(jx,iy),xlat(jx,iy),dlon(jx,iy),dlat(jx,iy)
     &     ,		coriol(jx,iy),xlandu(jx,iy),snowcv(jx,iy),topogm(jx,iy)
     &     ,		toposdgm(jx,iy),msfx(jx,iy),sigma2(kz),sigmaf(kz+1)
     &     ,		dsigma(kz)
      character*6 lgtype
      real	ptop,clat,clon,delx,grdfac,plat,plon
      integer 	igrads,ibigend
      common /LGRID2/ 	ptop,clat,clon,plat,plon,delx,grdfac
     &     ,		igrads,ibigend,lgtype
	   
!     Dummy arrays for interpolation
      real	dum1(jx,iy,npl)
      real    	dum2(jx,iy,npl,2)
      real    	c1(npl), c2(kz)
		
      integer 	nyrp,nmop
      real    	wt,w,pp,h2
 
!     Arrays for new calculation of P
      real      pa(jx,iy),za(jx,iy)
      real  	tlayer(jx,iy)
      real  	sst1(jx,iy), sst2(jx,iy)
      

      character*4 varname
      character*38 pathaddname
      data varname/'PHIS'/
      real,  	allocatable:: work(:,:)
      integer   inet1,ivar1,icode,status,counta,starta
      common /CAM85open1/	starta(10),counta(10)
      integer	lonid,latid,lonlen,latlen
      integer   checklat,checklon
      integer 	i,j,k,ii,jj
      integer	i0,i1,j0,j1
	integer	imin,imax
      logical 	there
 
      call CAM85(idate,idate1,xlon,xlat,glon,glat,jx,iy,
     &     i0,i1,j0,j1)   
      if(idate.eq.idate1) then
         pathaddname= '../DATA/CAM85/ccsm_ht.nc'
         inquire(file=pathaddname,exist=there)
         if(.not.there)then
            print*,pathaddname, 'is not available'
            stop
         endif
         status=nf_open(pathaddname,nf_nowrite,inet1)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_varid(inet1,varname,ivar1)
         if (status.ne.nf_noerr) call handle_err(status)
         write(*,*)inet1,ivar1,pathaddname,ivar1,varname
         status=nf_inq_dimid(inet1,'lon',lonid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimid(inet1,'lat',latid)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimlen(inet1,lonid,lonlen)
         if (status.ne.nf_noerr) call handle_err(status)
         status=nf_inq_dimlen(inet1,latid,latlen)
         if (status.ne.nf_noerr) call handle_err(status)
         imin=i0
         imax=i1
         if (i0.gt.i1) imax=i1+ilon
         checklon=(imax-imin)+1
         checklat=(j1-j0)+1
         if(checklon.ne.lonlen.or.checklat.ne.latlen) then
            print*, 'DOMAIN DIMENSIONS DO NOT MATCH'
            print*, 'LAT for 3D Variables = ',checklat
            print*, 'LAT for'  ,varname,'= ',latlen
            print*, 'LON for 3D Variables = ',checklon
            print*, 'LON for'  ,varname,'= ',lonlen
            STOP    'Check Dimensions in CCSM Data Files'
         endif
         allocate (work(lonlen,latlen))
         counta(1)=lonlen
         counta(2)=latlen
         counta(3)=1
         starta(1)=1
         starta(2)=1
         starta(3)=1
         status=nf_get_vara_real(inet1,ivar1,starta,counta,work)
         if (status.ne.nf_noerr) call handle_err(status)
         do jj=j0,j1
            j=jj-j0+1
            if(i0.gt.i1) then
               do ii=i0,ilon
                  i=ii-i0+1
                  ZSvar(ii,jj)=work(i,j)/9.80616
               enddo
               do ii=1,i1
                  i=ii+(ilon-i0)
                  ZSvar(ii,jj)=work(i,j)/9.80616
               enddo
            else
               do ii=i0,i1
                  i=ii-i0+1
                  ZSvar(ii,jj)=work(i,j)/9.80616
               enddo
            endif
         enddo
      endif
      write(*,*) 'Read in fields at Date: ',idate
      
      do K = 1,klev
         do J = 1,jlat
            do I = 1,ilon
               if(PSvar(I,J).GT.-9995.) then
                  PP3d(I,J,K) = PSvar(I,J)*0.5*(ak(k)+bk(k+1))
     &                 +          0.5*(ak(k)+bk(k+1))
               else
                  PP3d(I,J,K) = -9999.0
               endif
            enddo
         enddo
      enddo
      
      CALL HEIGHT(Hp,Hvar,Tvar,PSvar,PP3D,ZSvar,ilon,jlat,klev,pplev
     &     ,npl)

      CALL INTLIN(UP,Uvar,PSvar,PP3D,ilon,jlat,klev,pplev,npl)
      CALL INTLIN(VP,Vvar,PSvar,PP3D,ilon,jlat,klev,pplev,npl)

      CALL INTLOG(TP,Tvar,PSvar,PP3D,ilon,jlat,klev,pplev,npl)

      CALL HUMID1FV(Tvar,Qvar,PP3D,ilon,jlat,klev)
      CALL INTLIN(QP,Qvar,PSvar,PP3d,ilon,jlat,klev,pplev,npl)

      CALL BILINXCX(B3,B2,xlon,xlat,glon,glat,jx,iy,ilon,jlat,npl)
      CALL BILINXDT(D3,D2,dlon,dlat,glon,glat,jx,iy,ilon,jlat,npl)
      CALL UVROT4(U3,V3,dlon,dlat,clon,clat,grdfac,jx,iy,npl
     &     ,plon,plat,lgtype)

      CALL TOP2BTM(T3,jx,iy,npl)
      CALL TOP2BTM(Q3,jx,iy,npl)
      CALL TOP2BTM(H3,jx,iy,npl)
      CALL TOP2BTM(U3,jx,iy,npl)
      CALL TOP2BTM(V3,jx,iy,npl)

      CALL INTGTB(PA,ZA,tlayer,topogm,T3,H3
     &     ,sigmar,jx,iy,npl,dum1,dum2)
 
      CALL INTPSN(PS4,topogm,PA,ZA,tlayer,ptop,jx,iy)
      CALL P1P2(B3PD,PS4,jx,iy)

      CALL INTV3(TS4,T3,PS4,sigmar,ptop,jx,iy,npl,dum1)     
      CALL CAMCLNDR( idate, nyrp, nmop, wt )
      CALL MKSST(TS4,SST1,SST2,topogm,xlandu,jx,iy,nyrp,nmop,wt)
      CALL INTV1(U4,U3,B3PD,sigma2,sigmar,ptop,jx,iy,kz
     &     ,npl,1,1,dum2,c1,c2)
      CALL INTV1(V4,V3,B3PD,sigma2,sigmar,ptop,jx,iy,kz
     &     ,npl,1,1,dum2,c1,c2)
      CALL INTV2(T4,T3,PS4,sigma2,sigmar,ptop,jx,iy,kz
     &     ,npl,1,dum2,c1,c2)
 
      CALL INTV1(Q4,Q3,PS4,sigma2,sigmar,ptop,jx,iy,kz
     &     ,npl,1,0,dum2,c1,c2)
      CALL HUMID2FV(T4,Q4,PS4,ptop,sigma2,jx,iy,kz)

      CALL HYDROST(H4,T4,topogm,PS4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)

   
      CALL WRITEF(U4,V4,T4,Q4,PP,W,PS4,TS4,ptop,jx,iy,kz,idate
     &     ,iomega)
      return
      end

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      subroutine HEAD_CAM42
      implicit 	none
      
      integer	ilon,jlat,klev
      parameter (ilon=128)
      parameter (jlat=64)
      parameter (klev=26)
		
      integer	npl
      parameter	(npl=18)
      real    	glat, glon,sigma1,sigmar,pplev,ak,bk,psref
      common /globalCAM42/ glat(jlat),glon(ilon),sigma1(npl),sigmar(npl)
     &     ,		 pplev(npl),ak(klev+1),bk(klev+1),psref
      integer 	i,j,k,kr
      data glat/
     &  -87.8638,-85.0965,-82.3129,-79.5256,-76.7369,-73.9475,-71.1578,
     &  -68.3678,-65.5776,-62.7874,-59.9970,-57.2066,-54.4162,-51.6257,
     &  -48.8352,-46.0447,-43.2542,-40.4636,-37.6731,-34.8825,-32.0919,
     &  -29.3014,-26.5108,-23.7202,-20.9296,-18.1390,-15.3484,-12.5578,
     &  -9.76715,-6.97653,-4.18592,-1.39531, 1.39531, 4.18592, 6.97653,
     &   9.76715, 12.5578, 15.3484, 18.1390, 20.9296, 23.7202, 26.5108,
     &   29.3014, 32.0919, 34.8825, 37.6731, 40.4636, 43.2542, 46.0447,
     &   48.8352, 51.6257, 54.4162, 57.2066, 59.9970, 62.7874, 65.5776,
     &   68.3678, 71.1578, 73.9475, 76.7369, 79.5256, 82.3129, 85.0965,
     &   87.8638/
      do i = 1,ilon
         glon(i) = float(i-1)*2.8125
      enddo
      
      PPLEV(1)	=  30.
      PPLEV(2)	=  50.
      PPLEV(3)	=  70.
      PPLEV(4)	= 100.
      PPLEV(5)	= 150.	
      PPLEV(6)	= 200.
      PPLEV(7)	= 250.
      PPLEV(8)	= 300.
      PPLEV(9)	= 350.
      PPLEV(10)	= 420.
      PPLEV(11)	= 500.
      PPLEV(12)	= 600.
      PPLEV(13)	= 700.
      PPLEV(14)	= 780.
      PPLEV(15)	= 850.
      PPLEV(16)	= 920.
      PPLEV(17)	= 960.
      PPLEV(18)	=1000.
             
		
      do k = 1, npl
         sigmar(k) = pplev(k)*0.001
      enddo
    		
      do 116 k = 1, npl
         kr = npl-k+1 
         sigma1(k) = sigmar(kr)
 116  continue
      psref =   1000.
      
!     Ak is acatually Ak*Po

      ak(1) =   2.194067
      ak(2) =   4.895209
      ak(3) =   9.882418
      ak(4) =   18.05201
      ak(5) =   29.83724
      ak(6) =   44.62334
      ak(7) =   61.60587
      ak(8) =   78.51243
      ak(9) =   77.31271
      ak(10)=   75.90131
      ak(11)=   74.24086
      ak(12)=   72.28744
      ak(13)=   69.98932
      ak(14)=   67.28574
      ak(15)=   64.10509
      ak(16)=   60.36322
      ak(17)=   55.96111
      ak(18)=   50.78225
      ak(19)=   44.68960
      ak(20)=   37.52190
      ak(21)=   29.08949
      ak(22)=   20.84739
      ak(23)=   13.34443
      ak(24)=   7.084990
      ak(25)=   2.521360
      ak(26)=   0.0
      ak(27)=   0.0
                
      bk(1) =   0.0
      bk(2) =   0.0
      bk(3) =   0.0
      bk(4) =   0.0
      bk(5) =   0.0
      bk(6) =   0.0
      bk(7) =   0.0
      bk(8) =   0.0
      bk(9) =   0.0150530
      bk(10)=   0.0327622
      bk(11)=   0.0535962
      bk(12)=   0.0781062
      bk(13)=   0.1069411
      bk(14)=   0.1408637
      bk(15)=   0.1807720
      bk(16)=   0.2277220
      bk(17)=   0.2829562
      bk(18)=   0.3479364
      bk(19)=   0.4243822
      bk(20)=   0.5143168
      bk(21)=   0.6201202
      bk(22)=   0.7235355
      bk(23)=   0.8176768
      bk(24)=   0.8962153
      bk(25)=   0.9534761
      bk(26)=   0.9851122
      bk(27)=   1.

              
      return
      end

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
      subroutine HEAD_CAM85
      implicit 	none
      
      integer	ilon,jlat,klev
      parameter (ilon=256)
      parameter (jlat=128)
      parameter (klev=26)
		
      integer	npl
      parameter	(npl=18)
      real    	glat, glon,sigma1,sigmar,pplev,ak,bk,psref
      common /globalCAM85/ glat(jlat),glon(ilon),sigma1(npl),sigmar(npl)
     &     ,		 pplev(npl),ak(klev+1),bk(klev+1),psref
      integer 	i,j,k,kr
      data glat/
     & -88.928, -87.539, -86.141, -84.742, -83.343, -81.942,
     & -80.543, -79.142, -77.742, -76.341, -74.941, -73.539,
     & -72.139, -70.739, -69.338, -67.937, -66.536, -65.136,
     & -63.735, -62.334, -60.934, -59.534, -58.132, -56.731,
     & -55.331, -53.930, -52.529, -51.128, -49.728, -48.327,
     & -46.926, -45.525, -44.125, -42.724, -41.323, -39.922,
     & -38.522, -37.121, -35.720, -34.319, -32.919, -31.518,
     & -30.117, -28.716, -27.315, -25.915, -24.514, -23.113,
     & -21.712, -20.312, -18.911, -17.510, -16.109, -14.709,
     & -13.308, -11.907, -10.506,  -9.105,  -7.705,  -6.304,
     &  -4.903,  -3.502,  -2.102,  -0.701,   0.701,   2.102,
     &   3.502,   4.903,   6.304,   7.705,   9.105,  10.506,
     &  11.907,  13.308,  14.709,  16.109,  17.510,  18.911,
     &  20.312,  21.712,  23.113,  24.514,  25.915,  27.315,
     &  28.716,  30.117,  31.518,  32.919,  34.319,  35.720,
     &  37.121,  38.522,  39.922,  41.323,  42.724,  44.125,
     &  45.525,  46.926,  48.327,  49.728,  51.128,  52.529,
     &  53.930,  55.331,  56.731,  58.132,  59.533,  60.934,
     &  62.334,  63.735,  65.136,  66.536,  67.937,  69.338,
     &  70.739,  72.139,  73.539,  74.941,  76.341,  77.742,
     &  79.142,  80.543,  81.942,  83.343,  84.742,  86.141,
     &  87.539,  88.928/
 
      do i = 1,ilon
         glon(i) = float(i-1)*1.40625
      enddo
      
      PPLEV(1)	=  30.
      PPLEV(2)	=  50.
      PPLEV(3)	=  70.
      PPLEV(4)	= 100.
      PPLEV(5)	= 150.	
      PPLEV(6)	= 200.
      PPLEV(7)	= 250.
      PPLEV(8)	= 300.
      PPLEV(9)	= 350.
      PPLEV(10)	= 420.
      PPLEV(11)	= 500.
      PPLEV(12)	= 600.
      PPLEV(13)	= 700.
      PPLEV(14)	= 780.
      PPLEV(15)	= 850.
      PPLEV(16)	= 920.
      PPLEV(17)	= 960.
      PPLEV(18)	=1000.
             
		
      do k = 1, npl
         sigmar(k) = pplev(k)*0.001
      enddo
    		
      do 116 k = 1, npl
         kr = npl-k+1 
         sigma1(k) = sigmar(kr)
 116  continue
      psref =   1000.
      
!     Ak is acatually Ak*Po

      ak(1) =   2.194067
      ak(2) =   4.895209
      ak(3) =   9.882418
      ak(4) =   18.05201
      ak(5) =   29.83724
      ak(6) =   44.62334
      ak(7) =   61.60587
      ak(8) =   78.51243
      ak(9) =   77.31271
      ak(10)=   75.90131
      ak(11)=   74.24086
      ak(12)=   72.28744
      ak(13)=   69.98932
      ak(14)=   67.28574
      ak(15)=   64.10509
      ak(16)=   60.36322
      ak(17)=   55.96111
      ak(18)=   50.78225
      ak(19)=   44.68960
      ak(20)=   37.52190
      ak(21)=   29.08949
      ak(22)=   20.84739
      ak(23)=   13.34443
      ak(24)=   7.084990
      ak(25)=   2.521360
      ak(26)=   0.0
      ak(27)=   0.0
                
      bk(1) =   0.0
      bk(2) =   0.0
      bk(3) =   0.0
      bk(4) =   0.0
      bk(5) =   0.0
      bk(6) =   0.0
      bk(7) =   0.0
      bk(8) =   0.0
      bk(9) =   0.0150530
      bk(10)=   0.0327622
      bk(11)=   0.0535962
      bk(12)=   0.0781062
      bk(13)=   0.1069411
      bk(14)=   0.1408637
      bk(15)=   0.1807720
      bk(16)=   0.2277220
      bk(17)=   0.2829562
      bk(18)=   0.3479364
      bk(19)=   0.4243822
      bk(20)=   0.5143168
      bk(21)=   0.6201202
      bk(22)=   0.7235355
      bk(23)=   0.8176768
      bk(24)=   0.8962153
      bk(25)=   0.9534761
      bk(26)=   0.9851122
      bk(27)=   1.
              
      return
      end

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      subroutine INITDATE3
      implicit	none
      integer 	mdate
      common /DATENUM/	mdate(289280)
      integer 	mbase,nbase,nrec,nyear,mon,nday,i,m
      nrec=0
      do nyear=1948,2045
         mbase = nyear*1000000
         do mon=1,12
            mbase = mbase+10000
            if(mon.eq.1.or.mon.eq.3.or.mon.eq.5.or.mon.eq.7
     &           .or.mon.eq.8.or.mon.eq.10.or.mon.eq.12) then
               nday=31
            else if(mon.eq.4.or.mon.eq.6.or.mon.eq.9.or.mon.eq.11)then
               nday=30
            else
               nday=28
            endif
            nbase = mbase
            do I=1,nday
               nbase = nbase+100
               do m=1,4
                  nrec=nrec+1
                  if(m.eq.1) then
                     mdate(nrec)=nbase
                  else if(m.eq.2) then
                     mdate(nrec)=nbase+6
                  else if(m.eq.3) then
                     mdate(nrec)=nbase+12
                  else
                     mdate(nrec)=nbase+18
                  endif
               enddo
            enddo
         enddo
      enddo
      write(*,*) 'nrec = ',nrec
      return
      end
      
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!

      subroutine CAMCLNDR(mdate, nyrp, nmop, wt)
      implicit none
 
      integer mdate,nyrp,nmop
      real    wt
      integer lenmon(12), midmon(12), julmid(12), jprev(12)
 
      data    lenmon/31,28,31,30,31,30,31,31,30,31,30,31/
      data    midmon/16,14,16,15,16,15,16,16,15,16,15,16/
      integer idate,iyr,imo,iday,ileap,j,julday,nyr,nmo
      real    fnumer,fdemon
 
!     Initialize nmop, nyrp										!
      nmop = 1
      nyrp = 0
 
      idate = mdate/100
      iyr = idate/10000
      imo = ( idate - iyr*10000 )/100
      iday = mod(idate, 100)
      
      jprev(1) = 0
      do 100 j=2,12
         jprev(j)  = jprev(j-1) + lenmon(j-1)
 100  continue
      do 200 j=1,12
         julmid(j) = jprev(j) + midmon(j)
 200  continue
      julday = iday + jprev(imo)
      
      do 300  nyr=1000, 2100
         do 300  nmo=1,12
 
            if( (nyr.eq.iyr) .and. (julmid(nmo).gt.julday) ) goto 400
            if  (nyr.gt.iyr) goto 400
            
            nmop = nmo
            nyrp = nyr
            
 300     continue
         
 400  continue
      fnumer = float(julday  - julmid(nmop))
      if(fnumer.lt.0.) fnumer = fnumer + 365.
      fdemon = float(julmid(nmo) - julmid(nmop))
      if(fdemon.le.0.) fdemon = fdemon + 365.
      wt = fnumer / fdemon
      
      return
      end

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!

!     Error Handler for NETCDF Calls					
      subroutine HANDLE_ERR(STATUS)
      implicit none
      include 	'netcdf.inc'
      integer	status
      if (status.ne.nf_noerr) then
         print*, nf_strerror(status)
         stop 'Netcdf Error'
      endif
      return
      end 
	




