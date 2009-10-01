      implicit none
      character*90 fname
      real*4  a(288,181),b(100,100)
      real*4  lat(181),lon(288)
      real*4  lon0,lon1,lat0,lat1
      integer i,j,i0,i1,j0,j1,isize,jsize
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
      open(20,file='HT_SRF',form='unformatted'
     &       ,recl=isize*jsize*4,access='direct')
c
      open(10,file='HTgSRF',form='unformatted'
     &       ,recl=288*181*4,access='direct')
      read(10,rec=1) a
      if(i0.lt.i1) then
        write(20,rec=1) ((a(i,j),i=i0,i1),j=j0,j1)
      else
        do j=j0,j1
          do i=i0,288
            b(i-i0+1,j-j0+1)=a(i,j)
          enddo
          do i=1,i1
            b(i+isize-i1,j-j0+1)=a(i,j)
          enddo
        enddo
        write(20,rec=1) ((b(i,j),i=1,isize),j=1,jsize)
      endif
      close(10)
      close(20)
c
      stop
      end
