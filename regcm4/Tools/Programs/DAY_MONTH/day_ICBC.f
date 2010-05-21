      implicit none
      integer im,jn,nvar
      parameter (im=192,jn=208,nvar=18*5+2)
      real*4  b(im,jn),c(im,jn,nvar)
      character*10 in_fn
      character*2 monCHA(12)
      character*4 yr_CHA(1950:1977)
      data monCHA/'01','02','03','04','05','06'
     &            ,'07','08','09','10','11','12'/
      data yr_CHA/'1950','1951','1952','1953','1954'
     &           ,'1955','1956','1957','1958','1959'
     &           ,'1960','1961','1962','1963','1964'
     &           ,'1965','1966','1967','1968','1969'
     &           ,'1970','1971','1972','1973','1974'
     &           ,'1975','1976','1977'/
      character*8 outfn
      integer i,j,n,mrec,nrec,nday,mday,month,nyear,nrecord
      integer IDATE
      do nyear=1950,1977
      do month=1,12
         outfn=yr_CHA(nyear)//monCHA(month)//'01'
         open(20,file='ICBC'//outfn,form='unformatted'
     &          ,recl=(im-2)*(jn-2)*4,access='direct')
         mrec=0
         in_fn=yr_CHA(nyear)//monCHA(month)//'0100'
         open(10,file='../EH5OM_Run/ICBC/ICBC'//in_fn
     &          ,form='unformatted'
     &          ,recl=im*jn*4,access='direct')
         nrec=nvar+1
         if(month.eq.1.or.month.eq.3.or.month.eq.5 .or.
     &      month.eq.7.or.month.eq.8.or.month.eq.10.or.
     &      month.eq.12) then
            nrecord=31
         else if(month.eq.4.or.month.eq.6.or.month.eq.9.or.
     &           month.eq.11) then
            nrecord=30
         else
            nrecord=28
            if(mod(nyear,4).eq.0)   nrecord=nrecord+1
            if(mod(nyear,100).eq.0) nrecord=nrecord-1
            if(mod(nyear,400).eq.0) nrecord=nrecord+1
         endif
         do nday=1,nrecord
            do n=1,nvar
               do j=1,jn
               do i=1,im
                  c(i,j,n)=0.0
               enddo
               enddo
            enddo
            do mday=1,4
               nrec=nrec+1
               read(10,rec=nrec) IDATE
               write(*,*) 'nrec = ',nrec ,'  IDATE =',IDATE
               do n=1,nvar
                  nrec=nrec+1
                  read(10,rec=nrec) b
                  do j=1,jn
                  do i=1,im
                     c(i,j,n)=c(i,j,n)+b(i,j)
                  enddo
                  enddo
               enddo
            enddo
            do n=1,nvar
               do j=1,jn
               do i=1,im
                  c(i,j,n)=c(i,j,n)/4.
               enddo
               enddo
               mrec=mrec+1
               write(20,rec=mrec) ((c(i,j,n),i=2,im-1),j=2,jn-1)
            enddo
         enddo
         close(10)
         close(20)
      enddo
      enddo
      stop
      end
