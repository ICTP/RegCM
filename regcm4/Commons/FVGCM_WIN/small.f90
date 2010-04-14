!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      program small

      implicit none
!
! Local variables
!
      real(4) , dimension(288,181) :: a
      real(4) , dimension(100,100) :: b
      character(1) :: ctrl
      character(90) :: fname
      integer :: i , i0 , i1 , isize , j , j0 , j1 , jsize , month ,    &
               & mrec , n6hour , nrec , nmbr , nv , nyear
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
        open (20,file='CUTREG_01',form='unformatted',recl=isize*jsize*4,&
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
        nmbr = 30*4
        if ( month==1 .or. month==3 .or. month==5 .or. month==7 .or.    &
           & month==8 .or. month==10 .or. month==12 ) nmbr = 31*4
        if ( month==2 ) nmbr = 28*4
        if ( month==2 .and. mod(nyear,4)==0 ) nmbr = 29*4
        do n6hour = 1 , nmbr
          do nv = 1 , 131
            nrec = nrec + 1
            if ( .not.((nv>=1 .and. nv<=18) .or. (nv>=76 .and. nv<=93)  &
               & .or. (nv>=114 .and. nv<=131)) ) then
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
            end if
          end do
        end do
        close (10)
!
        write (*,*) 'continue or not: y/n'
        read (*,'(A)') ctrl
        if ( ctrl=='y' ) go to 50
      end if
      end program small
