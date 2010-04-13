      program mon_ICBC

      implicit none

      integer , parameter :: im = 192-2
      integer , parameter :: jn = 208-2
      integer , parameter :: nvar = 18*5+2

      real(4), dimension(im,jn) :: b
      real(4), dimension(im,jn,nvar) :: c

      character(10) :: in_fn
      character(8) :: outfn
      character(2) , dimension(12) :: monCHA
      character(4) , dimension(1950:1977) :: yr_CHA
      integer :: i , j , n , mrec , nrec , nday , mday , month , nyear ,&
                 nrecord , idate
!
      data monCHA/'01','02','03','04','05','06',                        &
                & '07','08','09','10','11','12'/
      data yr_CHA/'1950','1951','1952','1953','1954',                   &
            &     '1955','1956','1957','1958','1959',                   &
            &     '1960','1961','1962','1963','1964',                   &
            &     '1965','1966','1967','1968','1969',                   &
            &     '1970','1971','1972','1973','1974',                   &
            &     '1975','1976','1977'/
!
      do nyear = 1950 , 1951
        outfn = yr_CHA(nyear)//'.dat'
        open(20,file='ICBC'//outfn,form='unformatted',                  &
         &   recl=im*jn*4,access='direct')
        mrec = 0
        do month = 1 , 12
          in_fn = yr_CHA(nyear)//monCHA(month)//'01'
          open(10,file='ICBC'//in_fn,form='unformatted',                &
            &  recl=im*jn*4,access='direct')
          nrec = 0
          if ( month.eq.1 .or. month.eq.3 .or. month.eq.5 .or.          &
          &    month.eq.7 .or. month.eq.8 .or. month.eq.10.or.          &
          &    month.eq.12 ) then
            nrecord = 31
          else if ( month.eq.4 .or. month.eq.6 .or. month.eq.9 .or.     &
           &        month.eq.11 ) then
            nrecord = 30
          else
            nrecord = 28
            if( mod(nyear,4).eq.0 )   nrecord = nrecord+1
            if( mod(nyear,100).eq.0 ) nrecord = nrecord-1
            if( mod(nyear,400).eq.0 ) nrecord = nrecord+1
          end if

          c(:,:,:) = 0.0
          do nday = 1 , nrecord
            write (*,*) 'mrec = ' , mrec , ' nrec = ' , nrec
            do n = 1 , nvar
              nrec = nrec+1
              read(10,rec=nrec) b
              do j = 1 , jn
                do i = 1 , im
                  c(i,j,n) = c(i,j,n)+b(i,j)
                end do
              end do
            end do
          end do
          do n = 1 , nvar
            do j = 1 , jn
              do i = 1 , im
                c(i,j,n) = c(i,j,n)/float(nrecord)
              end do
            end do
            mrec = mrec+1
            write (20,rec=mrec) ((c(i,j,n),i=1,im),j=1,jn)
          end do
          close(10)
        end do
        close(20)
      end do
      end program mon_ICBC
