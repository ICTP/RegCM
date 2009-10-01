      implicit none
      character*90 fname
      character*12 oname
      integer*2  a(192,96),b(100,100)
      real*8  offset,scale
      real*4  lat(96),lon(192)
      real*4  lon0,lon1,lat0,lat1
      character*1 ctrl
      integer month,nyear,number
      integer i,j,i0,i1,j0,j1,isize,jsize,mrec,nrec
      integer n6hour,nv
c
      open(100,file='0WINDOW',form='formatted')
      write(*,*) 'Please input lon0,lon1,lat0,lat1'
      write(*,*) 'Note: lon0 < lon1, and lat0 < lat1'
      read(100,*) lon0,lon1,lat0,lat1
      do j=1,96
         lat(j)=-89.0625+(j-1)*1.875
      enddo
      do i=1,192
         lon(i)=(i-1)*1.875
      enddo
      j0=96
      j1=96
      do j=1,96
         if(lat(j).le.lat0) j0=j
         if(lat(j).lt.lat1) j1=j
      enddo
      if(j1.ne.96) j1=j1+1
c
      i0=193
      i1=193
      do i=1,192
         if(lon0.ge.0.) then
           if(lon(i).le.lon0) i0=i
         else
           if(lon(i)-360.0.le.lon0) i0=i 
         endif
         if(lon1.ge.0.) then
           if(lon(i).lt.lon1) i1=i 
         else
           if(lon(i)-360.0.lt.lon1) i1=i
         endif
      enddo
      if(i1.ne.193) i1=i1+1
      if(i0.eq.193) i0=1
      if(i1.eq.193) i1=1
      write(*,*) 'present longitude: ',lon(i0),lon(i1)
      write(*,*) 'present latitude : ',lat(j0),lat(j1)
      isize=i1-i0+1
      if(i1.lt.i0) isize=192-i0+1+i1
      if(isize.gt.100) write(*,*) 'Please enlarge array b'
      jsize=j1-j0+1
      if(jsize.gt.100) write(*,*) 'Please enlarge array b'
      write(*,*) isize,jsize
      write(*,*) 'continue or not: y/n'
      read(*,'(A)') ctrl
      if(ctrl.eq.'n') goto 200
c
 100  continue
      write(*,*) 'Please output file name:'
      read(*,'(A)')oname
      open(20,file=oname,form='unformatted'
     &       ,recl=isize*jsize*2+16,access='direct')
      mrec=0
      write(*,*) 'Please input file name:'
      read(*,'(A)')fname
      write(*,*) 'Please input which month, year'
      read(*,*) month,nyear
      open(10,file=fname,form='unformatted'
     &       ,recl=192*96*2+16,access='direct')
      nrec=0
      number=30*4
      if(month.eq.1.or.month.eq.3.or.month.eq.5.or.month.eq.7.or.
     &   month.eq.8.or.month.eq.10.or.month.eq.12) number=31*4
      if(month.eq.2) number=28*4
      if(month.eq.2.and.mod(nyear,4).eq.0) number=29*4
      if(month.eq.2.and.mod(nyear,100).eq.0) number=28*4
      if(month.eq.2.and.mod(nyear,400).eq.0) number=29*4
      do n6hour=1,number
        do nv=1,17*5
          nrec=nrec+1
      read(10,rec=nrec) offset,scale,((a(i,j),i=1,192),j=96,1,-1)
          mrec=mrec+1
          if(i0.lt.i1) then
c         write(*,*) 'i0,i1 = ',i0,i1,'  j0,j1 = ',j0,j1
      write(20,rec=mrec) offset,scale,((a(i,j),i=i0,i1),j=j0,j1)
          else
            do j=j0,j1
              do i=i0,192
                b(i-i0+1,j-j0+1)=a(i,j)
              enddo
              do i=1,i1
                b(i+isize-i1,j-j0+1)=a(i,j)
              enddo
            enddo
c           write(*,*) 'isize = ',isize,'  jsize = ',jsize
      write(20,rec=mrec) offset,scale,((b(i,j),i=1,isize),j=1,jsize)
          endif
          write(*,*) nv,nrec,mrec
        enddo
      enddo
      close(10)
      close(20)
c
      write(*,*) 'continue or not: y/n'
      read(*,'(A)') ctrl
      if(ctrl.eq.'y') goto 100
 200  continue
      stop
      end
