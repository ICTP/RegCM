      program smallps
      implicit none
!
! Local variables
!
      real(4) , dimension(288,181) :: a
      real(4) , dimension(100,100) :: b
      character(1) :: ctrl
      character(90) :: fname
      integer :: i , i0 , i1 , isize , j , j0 , j1 , jsize , month ,    &
               & mrec , n6hour , nrec , number , nyear
      real(4) , dimension(181) :: lat
      real(4) :: lat0 , lat1 , lon0 , lon1
      real(4) , dimension(288) :: lon
!
      write (*,*) 'Please input lon0,lon1,lat0,lat1'
      write (*,*) 'Note: lon0 < lon1, and lat0 < lat1'
      read (*,*) lon0 , lon1 , lat0 , lat1
      do j = 1 , 181
        lat(j) = -90. + (j-1)*1.
      end do
      do i = 1 , 288
        lon(i) = 0. + (i-1)*1.25
      end do
      j0 = 181
      j1 = 181
      do j = 1 , 181
        if ( lat(j)<=lat0 ) j0 = j
        if ( lat(j)<lat1 ) j1 = j
      end do
      if ( j1/=181 ) j1 = j1 + 1
!
      i0 = 289
      i1 = 289
      do i = 1 , 288
        if ( lon0>=0. ) then
          if ( lon(i)<=lon0 ) i0 = i
        else
          if ( lon(i)-360.0<=lon0 ) i0 = i
        end if
        if ( lon1>=0. ) then
          if ( lon(i)<lon1 ) i1 = i
        else
          if ( lon(i)-360.0<lon1 ) i1 = i
        end if
      end do
      if ( i1/=289 ) i1 = i1 + 1
      if ( i0==289 ) i0 = 1
      if ( i1==289 ) i1 = 1
      write (*,*) 'present longitude: ' , lon(i0) , lon(i1)
      write (*,*) 'present latitude : ' , lat(j0) , lat(j1)
      isize = i1 - i0 + 1
      if ( i1<i0 ) isize = 288 - i0 + 1 + i1
      if ( isize>100 ) write (*,*) 'Please enlarge array b'
      jsize = j1 - j0 + 1
      if ( jsize>100 ) write (*,*) 'Please enlarge array b'
      write (*,*) 'continue or not: y/n'
      read (*,'(A)') ctrl
      if ( ctrl/='n' ) then
        open (20,file='CUTREG',form='unformatted',recl=isize*jsize*4,   &
            & access='direct')
        mrec = 0
!
 50     continue
        write (*,*) 'Please input file name:'
        read (*,'(A)') fname
        write (*,*) 'Please input which month, year'
        read (*,*) month , nyear
        open (10,file=fname,form='unformatted',recl=288*181*4,          &
             &access='direct')
        nrec = 0
        number = 30*4
        if ( month==1 .or. month==3 .or. month==5 .or. month==7 .or.    &
           & month==8 .or. month==10 .or. month==12 ) number = 31*4
        if ( month==2 ) number = 28*4
        if ( month==2 .and. mod(nyear,4)==0 ) number = 29*4
        do n6hour = 1 , number
          nrec = nrec + 1
          read (10,rec=nrec) a
          mrec = mrec + 1
          if ( i0<i1 ) then
            write (20,rec=mrec) ((a(i,j),i=i0,i1),j=j0,j1)
          else
            do j = j0 , j1
              do i = i0 , 288
                b(i-i0+1,j-j0+1) = a(i,j)
              end do
              do i = 1 , i1
                b(i+isize-i1,j-j0+1) = a(i,j)
              end do
            end do
            write (20,rec=mrec) ((b(i,j),i=1,isize),j=1,jsize)
          end if
        end do
        close (10)
!
        write (*,*) 'continue or not: y/n'
        read (*,'(A)') ctrl
        if ( ctrl=='y' ) go to 50
      end if
      end program smallps
