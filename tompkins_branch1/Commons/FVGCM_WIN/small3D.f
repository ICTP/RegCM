      implicit none
      character*90 fname
      character*12 oname
      integer*2  a(288,181),b(100,100)
      real*8  offset,scale
      real*4  lat(181),lon(288)
      real*4  lon0,lon1,lat0,lat1
      character*1 ctrl
      integer month,nyear,number
      integer i,j,i0,i1,j0,j1,isize,jsize,mrec,nrec
      integer n6hour,nv
c
      write(*,*) 'Please input lon0,lon1,lat0,lat1'
      write(*,*) 'Note: lon0 < lon1, and lat0 < lat1'
      read(*,*) lon0,lon1,lat0,lat1
      do j=1,181
         lat(j)=-90.+(j-1)*1.
      enddo
      do i=1,288
         lon(i)=0.+(i-1)*1.25
      enddo
      j0=181
      j1=181
      do j=1,181
         if(lat(j).le.lat0) j0=j
         if(lat(j).lt.lat1) j1=j
      enddo
      if(j1.ne.181) j1=j1+1
c
      i0=289
      i1=289
      do i=1,288
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
      if(i1.ne.289) i1=i1+1
      if(i0.eq.289) i0=1
      if(i1.eq.289) i1=1
      write(*,*) 'present longitude: ',lon(i0),lon(i1)
      write(*,*) 'present latitude : ',lat(j0),lat(j1)
      isize=i1-i0+1
      if(i1.lt.i0) isize=288-i0+1+i1
      if(isize.gt.100) write(*,*) 'Please enlarge array b'
      jsize=j1-j0+1
      if(jsize.gt.100) write(*,*) 'Please enlarge array b'
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
     &       ,recl=288*181*2+16,access='direct')
      nrec=0
      number=30*4
      if(month.eq.1.or.month.eq.3.or.month.eq.5.or.month.eq.7.or.
     &   month.eq.8.or.month.eq.10.or.month.eq.12) number=31*4
      if(month.eq.2) number=28*4
      if(month.eq.2.and.mod(nyear,4).eq.0) number=29*4
      do n6hour=1,number
        do nv=1,74
          nrec=nrec+1
        if(nv.gt.1.and.nv.lt.74) then
          read(10,rec=nrec) offset,scale,a
          mrec=mrec+1
          if(i0.lt.i1) then
      write(20,rec=mrec) offset,scale,((a(i,j),i=i0,i1),j=j0,j1)
          else
            do j=j0,j1
              do i=i0,288
                b(i-i0+1,j-j0+1)=a(i,j)
              enddo
              do i=1,i1
                b(i+isize-i1,j-j0+1)=a(i,j)
              enddo
            enddo
      write(20,rec=mrec) offset,scale,((b(i,j),i=1,isize),j=1,jsize)
          endif
          endif
        enddo
      enddo
      close(10)
c
      write(*,*) 'continue or not: y/n'
      read(*,'(A)') ctrl
      if(ctrl.eq.'y') goto 100
 200  continue
      stop
      end
