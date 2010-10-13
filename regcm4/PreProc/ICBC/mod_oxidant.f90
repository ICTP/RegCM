      module  mod_oxidant
      use mod_dynparam

      implicit none

!
! PARAMETER definitions
! 
      integer, parameter :: ilon = 128, jlat = 64, ilev = 26, itime = 12

      integer nyear,month,nday,nhour
      integer k,l
      integer k0


   
      real, dimension(ilon)                 ::  t42lon
      real, dimension(jlat)                 ::  t42lat
      real, dimension(ilev)                 ::  t42hyam,t42hybm
      real, dimension(ilon,jlat,itime)      ::  t42ps
      real, dimension(ilon,ilev,jlat,itime) ::  t42oh,t42ho2,t42o3,t42no3,t42h2o2

      real, dimension(ilon,jlat)            :: POXID_2
      
      real, dimension(ilon,jlat,ilev)       :: oh2,ho22,o32,no32,h2o22




      real    :: PRCM,PMPI,PMPJ
      public  getmozart


      contains
      subroutine getmozart(idate)
      use mod_grid
      use mod_wrtoxd
      use mod_interp, only:bilinx2
      use netcdf
      implicit none

      integer         :: i,j,k,k0
      integer idate
      real, allocatable, dimension(:,:)     :: POXID_3
      real, allocatable, dimension(:,:,:)   :: oh3,ho23,o33,no33,h2o23
      integer, dimension(10) ::  start,count
      integer, dimension(10) ::  dimids
      integer :: dimid,ncid,status,xtype,ndim,natts,len

      allocate(POXID_3(jx,iy))

      allocate(oh3(jx,iy,ilev))
      allocate(ho23(jx,iy,ilev))
      allocate(o33(jx,iy,ilev))
      allocate(no33(jx,iy,ilev))
      allocate(h2o23(jx,iy,ilev))


!      call finddate_icbc(nstart,globidate1)
       write(*,*)'START',iy,jx


         status=nf90_open                                  &
         ('/lustre/fsolmon/regcm/RegCM-4.0.OX/PreProc/ICBC/oxid_3d_64x128_L26_c030722.nc', &
                       nf90_nowrite,ncid)
         if ( status/=nf90_noerr ) write (*,*) nf90_strerror(status)

         do i=1,128
            t42lon(i) = float(i-1)*2.8125
         enddo

         t42lat(1)  = -87.86380       
         t42lat(2)  = -85.09653       
         t42lat(3)  = -82.31291       
         t42lat(4)  = -79.52560
         t42lat(5)  = -76.73690       
         t42lat(6)  = -73.94752       
         t42lat(7)  = -71.15775       
         t42lat(8)  = -68.36776
         t42lat(9)  = -65.57761       
         t42lat(10) = -62.78735       
         t42lat(11) = -59.99702       
         t42lat(12) = -57.20663
         t42lat(13) = -54.41620       
         t42lat(14) = -51.62573       
         t42lat(15) = -48.83524       
         t42lat(16) = -46.04473
         t42lat(17) = -43.25420       
         t42lat(18) = -40.46365       
         t42lat(19) = -37.67309       
         t42lat(20) = -34.88252
         t42lat(21) = -32.09195       
         t42lat(22) = -29.30136       
         t42lat(23) = -26.51077       
         t42lat(24) = -23.72017
         t42lat(25) = -20.92957       
         t42lat(26) = -18.13897       
         t42lat(27) = -15.34836       
         t42lat(28) = -12.55776
         t42lat(29) = -9.767145       
         t42lat(30) = -6.976533       
         t42lat(31) = -4.185921       
         t42lat(32) = -1.395307
         t42lat(33) = 1.395307        
         t42lat(34) = 4.185921        
         t42lat(35) = 6.976533        
         t42lat(36) = 9.767145
         t42lat(37) = 12.55776        
         t42lat(38) = 15.34836        
         t42lat(39) = 18.13897        
         t42lat(40) = 20.92957
         t42lat(41) = 23.72017        
         t42lat(42) = 26.51077        
         t42lat(43) = 29.30136        
         t42lat(44) = 32.09195
         t42lat(45) = 34.88252        
         t42lat(46) = 37.67309        
         t42lat(47) = 40.46365        
         t42lat(48) = 43.25420
         t42lat(49) = 46.04473        
         t42lat(50) = 48.83524        
         t42lat(51) = 51.62573        
         t42lat(52) = 54.41620
         t42lat(53) = 57.20663        
         t42lat(54) = 59.99702        
         t42lat(55) = 62.78735        
         t42lat(56) = 65.57761
         t42lat(57) = 68.36776        
         t42lat(58) = 71.15775        
         t42lat(59) = 73.94752        
         t42lat(60) = 76.73690
         t42lat(61) = 79.52560        
         t42lat(62) = 82.31291        
         t42lat(63) = 85.09653        
         t42lat(64) = 87.86380

!c  t42lev:
!c    3.544638        7.388813       13.96721        23.94463
!c   37.23029        53.11460        70.05914        85.43912
!c  100.5147        118.2503        139.1154        163.6620
!c  192.5399        226.5132        266.4810        313.5013
!c  368.8180        433.8952        510.4553        600.5241
!c  696.7963        787.7020        867.1607        929.6490
!c  970.5548        992.5561

         t42hyam(1)  =  3.5446379
         t42hyam(2)  =  7.3888130
         t42hyam(3)  = 13.967210
         t42hyam(4)  = 23.944629
         t42hyam(5)  = 37.230290
         t42hyam(6)  = 53.114600
         t42hyam(7)  = 70.059143
         t42hyam(8)  = 77.912569
         t42hyam(9)  = 76.607011
         t42hyam(10) = 75.071082
         t42hyam(11) = 73.264152
         t42hyam(12) = 71.138389
         t42hyam(13) = 68.637542
         t42hyam(14) = 65.695412
         t42hyam(15) = 62.234160
         t42hyam(16) = 58.162171
         t42hyam(17) = 53.371690
         t42hyam(18) = 47.735929
         t42hyam(19) = 41.105751
         t42hyam(20) = 33.305701
         t42hyam(21) = 24.968440
         t42hyam(22) = 17.095910
         t42hyam(23) = 10.214710
         t42hyam(24) =  4.8031751
         t42hyam(25) =  1.2606850
         t42hyam(26) =  0.0000000

         t42hybm(1)  = 0.0000000
         t42hybm(2)  = 0.0000000
         t42hybm(3)  = 0.0000000
         t42hybm(4)  = 0.0000000     
         t42hybm(5)  = 0.0000000
         t42hybm(6)  = 0.0000000
         t42hybm(7)  = 0.0000000
         t42hybm(8)  = 7.5265500E-03 
         t42hybm(9)  = 2.3907680E-02
         t42hybm(10) = 4.3179251E-02
         t42hybm(11) = 6.5851226E-02
         t42hybm(12) = 9.2523657E-02 
         t42hybm(13) = 0.1239024
         t42hybm(14) = 0.1608178
         t42hybm(15) = 0.2042469
         t42hybm(16) = 0.2553391     
         t42hybm(17) = 0.3154463
         t42hybm(18) = 0.3861593
         t42hybm(19) = 0.4693495
         t42hybm(20) = 0.5672184     
         t42hybm(21) = 0.6718278
         t42hybm(22) = 0.7706061
         t42hybm(23) = 0.8569460
         t42hybm(24) = 0.9248458     
         t42hybm(25) = 0.9692941
         t42hybm(26) = 0.9925561
!C   t42time
!C  15.5    45.    74.5    105.    135.5   166.   196.5
!C  227.5  258.   288.5    319.    349.5

         do i=4,10
            start(i) = 0
         enddo

         do i=1,3
            start(i) = 1
         enddo

         count(1) = 128
         count(2) = 64
         count(3) = 12

!C  Retrieve data for Variable 'PS'
         status=nf90_get_var(ncid,13,t42PS,start,count)
!C  Retrieve data for Variable 'OH', 'HO2', 'O3', 'NO3', 'H2O2'
         start(4) = 1
         count(1) = 128
         count(2) = 26
         count(3) = 64
         count(4) = 12

         status=nf90_get_var(ncid,14,t42OH,start,count)
         status=nf90_get_var(ncid,15,t42HO2,start,count)
         status=nf90_get_var(ncid,16,t42O3,start,count)
         status=nf90_get_var(ncid,17,t42NO3,start,count)
         status=nf90_get_var(ncid,18,t42H2O2,start,count)
         
!      nyear=idate1/1000000
!      month=(idate1-nyear*1000000)/10000
!      nday =(idate1-nyear*1000000-month*10000)/100
!      nhour=MOD(idate1,100)
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100

      write(*,*)MONTH,nday,idate


         do j=1,64
         do i=1,128
           POXID_2(i,j) = T42PS(i,j,MONTH)*0.01
         enddo
         enddo


         IF((MONTH.eq.1.or.MONTH.eq.3..or.MONTH.eq.5.or.MONTH.eq.7 &
         .or.MONTH.eq.8.or.MONTH.eq.10).and.NDAY.GE.16) THEN
            IF(MONTH.eq.1) THEN
               do k=1,26
               do j=1,64
               do i=1,128
                  oh2(i,j,k) =                                    &
                      T42OH(i,k,j,MONTH)*(1.-float(NDAY-16)/30.)  &
                     +T42OH(i,k,j,MONTH+1)*(float(NDAY-16)/30.)  
                  ho22(i,j,k)=                                     &  
                      T42HO2(i,k,j,MONTH)*(1.-float(NDAY-16)/30.)  &
                     +T42HO2(i,k,j,MONTH+1)*(float(NDAY-16)/30.)  
                  o32(i,j,k) =                                     & 
                      T42O3(i,k,j,MONTH)*(1.-float(NDAY-16)/30.)   &
                     +T42O3(i,k,j,MONTH+1)*(float(NDAY-16)/30.)   
                  no32(i,j,k)=                                     &  
                       T42NO3(i,k,j,MONTH)*(1.-float(NDAY-16)/30.) &
                      +T42NO3(i,k,j,MONTH+1)*(float(NDAY-16)/30.)  
                  h2o22(i,j,k)=                                    & 
                      T42H2O2(i,k,j,MONTH)*(1.-float(NDAY-16)/30.) &
                     +T42H2O2(i,k,j,MONTH+1)*(float(NDAY-16)/30.)  
               enddo
               enddo
               enddo
            ELSE IF(MONTH.eq.3.or.MONTH.eq.5.or.                   & 
                    MONTH.eq.8.or.MONTH.eq.10) THEN                 
               do k=1,26
               do j=1,64
               do i=1,128
                  oh2(i,j,k) =                                     &
                       T42OH(i,k,j,MONTH)*(1.-float(NDAY-16)/31.)  &
                      +T42OH(i,k,j,MONTH+1)*(float(NDAY-16)/31.)   
                  ho22(i,j,k)=                                     &
                       T42HO2(i,k,j,MONTH)*(1.-float(NDAY-16)/31.) & 
                      +T42HO2(i,k,j,MONTH+1)*(float(NDAY-16)/31.)  
                  o32(i,j,k) =                                     &
                       T42O3(i,k,j,MONTH)*(1.-float(NDAY-16)/31.)  &
                      +T42O3(i,k,j,MONTH+1)*(float(NDAY-16)/31.)   
                  no32(i,j,k)=                                     &
                      T42NO3(i,k,j,MONTH)*(1.-float(NDAY-16)/31.)  &
                     +T42NO3(i,k,j,MONTH+1)*(float(NDAY-16)/31.)  
                  h2o22(i,j,k)=                                    &
                      T42H2O2(i,k,j,MONTH)*(1.-float(NDAY-16)/31.) &
                     +T42H2O2(i,k,j,MONTH+1)*(float(NDAY-16)/31.)  
               enddo
               enddo
               enddo

            ELSE IF(MONTH.eq.7) THEN
               do k=1,26
               do j=1,64
               do i=1,128
                  oh2(i,j,k) =                                     &
                      T42OH(i,k,j,MONTH)*(1.-float(NDAY-16)/32.)  &
                     +T42OH(i,k,j,MONTH+1)*(float(NDAY-16)/32.)   
                  ho22(i,j,k)=                                     &
                      T42HO2(i,k,j,MONTH)*(1.-float(NDAY-16)/32.) & 
                     +T42HO2(i,k,j,MONTH+1)*(float(NDAY-16)/32.)  
                  o32(i,j,k) =                                     &
                      T42O3(i,k,j,MONTH)*(1.-float(NDAY-16)/32.)  &
                     +T42O3(i,k,j,MONTH+1)*(float(NDAY-16)/32.)   
                  no32(i,j,k)=                                     &
                      T42NO3(i,k,j,MONTH)*(1.-float(NDAY-16)/32.) & 
                     +T42NO3(i,k,j,MONTH+1)*(float(NDAY-16)/32.)  
                  h2o22(i,j,k)=                                    &
                     T42H2O2(i,k,j,MONTH)*(1.-float(NDAY-16)/32.) & 
                    +T42H2O2(i,k,j,MONTH+1)*(float(NDAY-16)/32.)  
               enddo
               enddo
               enddo
            ENDIF

         ELSE IF(MONTH.eq.2.and.NDAY.GE.15) THEN
            do k=1,26
            do j=1,64
            do i=1,128
               oh2(i,j,k) =                                       &
                   T42OH(i,k,j,MONTH)*(1.-float(NDAY-16)/30.)    &
                  +T42OH(i,k,j,MONTH+1)*(float(NDAY-16)/30.)     
               ho22(i,j,k)=                                       & 
                   T42HO2(i,k,j,MONTH)*(1.-float(NDAY-16)/30.)   &
                  +T42HO2(i,k,j,MONTH+1)*(float(NDAY-16)/30.)    
               o32(i,j,k) =                                       & 
                   T42O3(i,k,j,MONTH)*(1.-float(NDAY-16)/30.)    &
                  +T42O3(i,k,j,MONTH+1)*(float(NDAY-16)/30.)     
               no32(i,j,k)=                                       &
                   T42NO3(i,k,j,MONTH)*(1.-float(NDAY-16)/30.)   &
                  +T42NO3(i,k,j,MONTH+1)*(float(NDAY-16)/30.)    
               h2o22(i,j,k)=                                      &
                   T42H2O2(i,k,j,MONTH)*(1.-float(NDAY-16)/30.)  &
                  +T42H2O2(i,k,j,MONTH+1)*(float(NDAY-16)/30.)   
            enddo
            enddo
            enddo
         ELSE IF((MONTH.eq.4..or.MONTH.eq.6.or.MONTH.eq.9.or.     &
                 MONTH.eq.11).and.NDAY.GE.16) THEN
            do k=1,26
            do j=1,64
            do i=1,128
               oh2(i,j,k) =                                       & 
                   T42OH(i,k,j,MONTH)*(1.-float(NDAY-16)/31.)    &
                  +T42OH(i,k,j,MONTH+1)*(float(NDAY-16)/31.)     
               ho22(i,j,k)=                                       & 
                   T42HO2(i,k,j,MONTH)*(1.-float(NDAY-16)/31.)   & 
                  +T42HO2(i,k,j,MONTH+1)*(float(NDAY-16)/31.)
               o32(i,j,k) =                                       &
                   T42O3(i,k,j,MONTH)*(1.-float(NDAY-16)/31.)    &
                  +T42O3(i,k,j,MONTH+1)*(float(NDAY-16)/31.)
               no32(i,j,k)=                                       &
                   T42NO3(i,k,j,MONTH)*(1.-float(NDAY-16)/31.)   & 
                  +T42NO3(i,k,j,MONTH+1)*(float(NDAY-16)/31.)
               h2o22(i,j,k)=                                      &
                   T42H2O2(i,k,j,MONTH)*(1.-float(NDAY-16)/31.)  &
                  +T42H2O2(i,k,j,MONTH+1)*(float(NDAY-16)/31.)
            enddo
            enddo
            enddo
         ELSE IF(MONTH.eq.12..and.NDAY.GE.16) THEN
            do k=1,26
            do j=1,64
            do i=1,128
               oh2(i,j,k) =                                       &
                   T42OH(i,k,j,12)*(1.-float(NDAY-16)/32.)       &
                  +T42OH(i,k,j,1)*(float(NDAY-16)/32.)
               ho22(i,j,k)=                                       & 
                   T42HO2(i,k,j,12)*(1.-float(NDAY-16)/32.)      &
                  +T42HO2(i,k,j,1)*(float(NDAY-16)/32.)
               o32(i,j,k) =                                       &
                   T42O3(i,k,j,12)*(1.-float(NDAY-16)/32.)       &
                  +T42O3(i,k,j,1)*(float(NDAY-16)/32.)
               no32(i,j,k)=                                       & 
                   T42NO3(i,k,j,12)*(1.-float(NDAY-16)/32.)      &
                  +T42NO3(i,k,j,1)*(float(NDAY-16)/32.)
               h2o22(i,j,k)=                                      &
                   T42H2O2(i,k,j,12)*(1.-float(NDAY-16)/32.)     & 
                  +T42H2O2(i,k,j,1)*(float(NDAY-16)/32.)
            enddo
            enddo
            enddo
         ELSE IF((MONTH.eq.3.or.MONTH.eq.5.or.MONTH.eq.7.or.      &
                 MONTH.eq.8.or.MONTH.eq.10.or.MONTH.eq.12)       &
                .and.NDAY.LT.16) THEN
            IF(MONTH.eq.3) THEN
               do k=1,26
               do j=1,64
               do i=1,128
                  oh2(i,j,k) =                                    &
                      T42OH(i,k,j,MONTH-1)*(float(16-NDAY)/30.)  &
                     +T42OH(i,k,j,MONTH)*(1.-float(16-NDAY)/30.)
                  ho22(i,j,k)=                                    & 
                     T42HO2(i,k,j,MONTH-1)*(float(16-NDAY)/30.)  & 
                    +T42HO2(i,k,j,MONTH)*(1.-float(16-NDAY)/30.)  
                  o32(i,j,k) =                                    &
                     T42O3(i,k,j,MONTH-1)*(float(16-NDAY)/30.)   &
                    +T42O3(i,k,j,MONTH)*(1.-float(16-NDAY)/30.)
                  no32(i,j,k)=                                    &
                     T42NO3(i,k,j,MONTH-1)*(float(16-NDAY)/30.)  &
                    +T42NO3(i,k,j,MONTH)*(1.-float(16-NDAY)/30.)
                  h2o22(i,j,k)=                                   &
                     T42H2O2(i,k,j,MONTH-1)*(float(16-NDAY)/30.) &
                    +T42H2O2(i,k,j,MONTH)*(1.-float(16-NDAY)/30.)
               enddo
               enddo
               enddo
            ELSE IF(MONTH.eq.5.or.MONTH.eq.7.or.MONTH.eq.10.or.   &
                   MONTH.eq.12) THEN
               do k=1,26
               do j=1,64
               do i=1,128
                  oh2(i,j,k) =                                    & 
                      T42OH(i,k,j,MONTH-1)*(float(16-NDAY)/31.)  &
                     +T42OH(i,k,j,MONTH)*(1.-float(16-NDAY)/31.)
                  ho22(i,j,k)=                                    &
                     T42HO2(i,k,j,MONTH-1)*(float(16-NDAY)/31.)  &
                    +T42HO2(i,k,j,MONTH)*(1.-float(16-NDAY)/31.)
                  o32(i,j,k) =                                    & 
                     T42O3(i,k,j,MONTH-1)*(float(16-NDAY)/31.)   &
                    +T42O3(i,k,j,MONTH)*(1.-float(16-NDAY)/31.)
                  no32(i,j,k)=                                    & 
                     T42NO3(i,k,j,MONTH-1)*(float(16-NDAY)/31.)  &
                    +T42NO3(i,k,j,MONTH)*(1.-float(16-NDAY)/31.)
                  h2o22(i,j,k)=                                   & 
                     T42H2O2(i,k,j,MONTH-1)*(float(16-NDAY)/31.) &
                    +T42H2O2(i,k,j,MONTH)*(1.-float(16-NDAY)/31.)
               enddo
               enddo
               enddo
            ELSE IF(MONTH.eq.8) THEN
               do k=1,26
               do j=1,64
               do i=1,128
                  oh2(i,j,k) =                                      & 
                      T42OH(i,k,j,MONTH-1)*(float(16-NDAY)/32.)    &
                     +T42OH(i,k,j,MONTH)*(1.-float(16-NDAY)/32.)
                  ho22(i,j,k)=                                      &
                     T42HO2(i,k,j,MONTH-1)*(float(16-NDAY)/32.)    &
                    +T42HO2(i,k,j,MONTH)*(1.-float(16-NDAY)/32.)
                  o32(i,j,k) =                                      &
                     T42O3(i,k,j,MONTH-1)*(float(16-NDAY)/32.)     & 
                    +T42O3(i,k,j,MONTH)*(1.-float(16-NDAY)/32.)
                  no32(i,j,k)=                                      &
                     T42NO3(i,k,j,MONTH-1)*(float(16-NDAY)/32.)    &
                    +T42NO3(i,k,j,MONTH)*(1.-float(16-NDAY)/32.)
                  h2o22(i,j,k)=                                     &
                     T42H2O2(i,k,j,MONTH-1)*(float(16-NDAY)/32.)   &
                    +T42H2O2(i,k,j,MONTH)*(1.-float(16-NDAY)/32.)
               enddo
               enddo
               enddo
            ENDIF
         ELSE IF(MONTH.eq.2.and.NDAY.LT.15) THEN
            do k=1,26
            do j=1,64
            do i=1,128
               oh2(i,j,k) =                                        &
                   T42OH(i,k,j,MONTH-1)*(float(15-NDAY)/30.)      &
                  +T42OH(i,k,j,MONTH)*(1.-float(15-NDAY)/30.)    
               ho22(i,j,k)=                                        & 
                  T42HO2(i,k,j,MONTH-1)*(float(15-NDAY)/30.)      &
                 +T42HO2(i,k,j,MONTH)*(1.-float(15-NDAY)/30.)
               o32(i,j,k) =                                        &
                  T42O3(i,k,j,MONTH-1)*(float(15-NDAY)/30.)       &
                 +T42O3(i,k,j,MONTH)*(1.-float(15-NDAY)/30.)
               no32(i,j,k)=                                        &
                  T42NO3(i,k,j,MONTH-1)*(float(15-NDAY)/30.)      &
                 +T42NO3(i,k,j,MONTH)*(1.-float(15-NDAY)/30.)
               h2o22(i,j,k)=                                       & 
                  T42H2O2(i,k,j,MONTH-1)*(float(15-NDAY)/30.)     &
                 +T42H2O2(i,k,j,MONTH)*(1.-float(15-NDAY)/30.)
            enddo
            enddo
            enddo
         ELSE IF((MONTH.eq.4..or.MONTH.eq.6.or.MONTH.eq.9.or.      &
                 MONTH.eq.11).and.NDAY.LT.16) THEN
                 write(*,*)'hi'
            do k=1,26
            do j=1,64
            do i=1,128
               oh2(i,j,k) =                                        & 
                   T42OH(i,k,j,MONTH-1)*(float(16-NDAY)/31.)      &
                  +T42OH(i,k,j,MONTH)*(1.-float(16-NDAY)/31.)
               ho22(i,j,k)=                                        &
                  T42HO2(i,k,j,MONTH-1)*(float(16-NDAY)/31.)      &
                 +T42HO2(i,k,j,MONTH)*(1.-float(16-NDAY)/31.)
               o32(i,j,k) =                                        &
                  T42O3(i,k,j,MONTH-1)*(float(16-NDAY)/31.)       &
                 +T42O3(i,k,j,MONTH)*(1.-float(16-NDAY)/31.)
               no32(i,j,k)=                                        &
                  T42NO3(i,k,j,MONTH-1)*(float(16-NDAY)/31.)      &
                 +T42NO3(i,k,j,MONTH)*(1.-float(16-NDAY)/31.)
               h2o22(i,j,k)=                                       &
                  T42H2O2(i,k,j,MONTH-1)*(float(16-NDAY)/31.)     &
                 +T42H2O2(i,k,j,MONTH)*(1.-float(16-NDAY)/31.)
            enddo
            enddo
            enddo
         ELSE IF(MONTH.eq.1.and.NDAY.LT.16) THEN
            do k=1,26
            do j=1,64
            do i=1,128
               oh2(i,j,k) = T42OH(i,k,j,12)*(float(16-NDAY)/32.)   &
                         + T42OH(i,k,j,1)*(1.-float(16-NDAY)/32.)
               ho22(i,j,k)=T42HO2(i,k,j,12)*(float(16-NDAY)/32.)   &
                         +T42HO2(i,k,j,1)*(1.-float(16-NDAY)/32.)
               o32(i,j,k) = T42O3(i,k,j,12)*(float(16-NDAY)/32.)   &
                         + T42O3(i,k,j,1)*(1.-float(16-NDAY)/32.)
               no32(i,j,k)=T42NO3(i,k,j,12)*(float(16-NDAY)/32.)   &
                         +T42NO3(i,k,j,1)*(1.-float(16-NDAY)/32.)
               h2o22(i,j,k)=T42H2O2(i,k,j,12)*(float(16-NDAY)/32.) &
                        +T42H2O2(i,k,j,1)*(1.-float(16-NDAY)/32.)
            enddo
            enddo
            enddo
         ENDIF



         CALL bilinx2(oh3,oh2,XLON,XLAT,T42LON,T42LAT,               &
                     128,64,iy,jx,26)
         CALL bilinx2(ho23,ho22,XLON,XLAT,T42LON,T42LAT,             &
                     128,64,iy,jx,26)
         CALL bilinx2(o33,o32,XLON,XLAT,T42LON,T42LAT,               &
                     128,64,iy,jx,26)
         CALL bilinx2(no33,no32,XLON,XLAT,T42LON,T42LAT,             &
                     128,64,iy,jx,26)
         CALL bilinx2(h2o23,h2o22,XLON,XLAT,T42LON,T42LAT,           & 
                     128,64,iy,jx,26)
         CALL bilinx2(POXID_3,POXID_2,XLON,XLAT,T42LON,T42LAT,       &  
                     128,64,iy,jx,1)
         do i=1,iy
         do j=1,jx
            do l=1,KZ
!               prcm=((PS4(j,i)-PTOP)*SIGMA2(l)+PTOP)*10.
                prcm=((POXID_3(j,i)*0.1-PTOP)*SIGMA2(l)+PTOP)*10.
!                  if(i .eq. 20 .and. j .eq. 20)then
!                  write(*,*)POXID_3(j,i),SIGMA2(l),prcm,l
!                  end if

               k0 = -1
               do k=26,1,-1
                  pmpi= POXID_3(j,i)*t42hybm(k)+t42hyam(k)
                  k0=k
                  if(prcm.gt.pmpi) goto 1002
               enddo
 1002          continue
               if(k0.eq.26) then
                  pmpj=POXID_3(j,i)*t42hybm(26-1)+t42hyam(26-1)
                  pmpi=POXID_3(j,i)*t42hybm(26)+t42hyam(26)
!                if(prcm .lt. pmpi .or. pmpi .lt. pmpj)then
!                  write(*,*)'NVE',k0,POXID_3(j,i),t42hybm(26)
!                end if
                  oh4(j,i,l)=oh3(j,i,26)                           &
                          +(oh3(j,i,26)-oh3(j,i,26-1))            &
                               *(prcm-pmpi)/(pmpi-pmpj)
                  ho24(j,i,l)=ho23(j,i,26)                         &
                           +(ho23(j,i,26)-ho23(j,i,26-1))         &
                               *(prcm-pmpi)/(pmpi-pmpj)
                  o34(j,i,l)=o33(j,i,26)                           & 
                          +(o33(j,i,26)-o33(j,i,26-1))            &
                               *(prcm-pmpi)/(pmpi-pmpj)
                  no34(j,i,l)=no33(j,i,26)                         &
                           +(no33(j,i,26)-no33(j,i,26-1))         &
                               *(prcm-pmpi)/(pmpi-pmpj)
                  h2o24(j,i,l)=h2o23(j,i,26)                       &
                            +(h2o23(j,i,26)-h2o23(j,i,26-1))      & 
                               *(prcm-pmpi)/(pmpi-pmpj)
               else if(k0.ge.1) then
                  pmpj=POXID_3(j,i)*t42hybm(k0)+t42hyam(k0)
                  pmpi=POXID_3(j,i)*t42hybm(k0+1)+t42hyam(k0+1)
!                if(prcm .lt. pmpi .or. pmpi .lt. pmpj)then
!                  write(*,*)'NVE',k0
!                end if
                  oh4(j,i,l)=(oh3(j,i,k0+1)*(prcm-pmpj)            & 
                            +oh3(j,i,k0)*(pmpi-prcm))             &
                           /(pmpi-pmpj)
                  ho24(j,i,l)=(ho23(j,i,k0+1)*(prcm-pmpj)          &
                             +ho23(j,i,k0)*(pmpi-prcm))           &
                            /(pmpi-pmpj)
                  o34(j,i,l)=(o33(j,i,k0+1)*(prcm-pmpj)            &
                            +o33(j,i,k0)*(pmpi-prcm))             &
                           /(pmpi-pmpj)
                  no34(j,i,l)=(no33(j,i,k0+1)*(prcm-pmpj)          &
                             +no33(j,i,k0)*(pmpi-prcm))           &
                            /(pmpi-pmpj)
                  h2o24(j,i,l)=(h2o23(j,i,k0+1)*(prcm-pmpj)        &
                              +h2o23(j,i,k0)*(pmpi-prcm))         &
                             /(pmpi-pmpj)
!                if(i .eq. 20 .and. j .eq.  20)then
!                write(*,*)o33(i,j,k0+1),l,k0+1
!                end if
               endif
            enddo
         enddo
         enddo
      do k=18,1,-1
      write(*,'(4(f10.4,1x))')t42o3(5,k,50,month)*1e9,o32(5,50,k)*1e9,          &
        o33(32,17,k)*1e9, o34(32,17,k)*1e9
!      write(*,'(4(f11.2,1x))')t42oh(5,k,50,month),oh2(5,50,k),          &
!        oh3(32,17,k), oh4(32,17,k)
      end do

         call writeox(ptop,idate)

         end subroutine getmozart

         end module mod_oxidant

