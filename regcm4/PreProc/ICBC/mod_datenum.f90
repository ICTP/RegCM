      module mod_datenum

      implicit none

      integer , dimension(289280) :: mdate

      contains

      subroutine finddate(npos,idate)
        implicit none
!
!       Dummy arguments
!
        integer :: idate , npos
        intent (in) idate
        intent (out) npos
!
!       Local variables
!
        integer :: i
!
        i = 0
 100    continue
        i = i + 1
        if ( mdate(i)==idate ) then
          npos = i
          go to 99999
        end if
        if ( i<=29824 ) go to 100
        write (*,*) 'ERROR IN FINDDATE'
        stop
99999   continue
      end subroutine finddate

      subroutine initdate
        implicit none
!
! Local variables
!
      integer :: i , m , mbase , mon , nbase , nday , nrec , nyear
!
        nrec = 0
        do nyear = 1989 , 2009
          mbase = nyear*1000000
          do mon = 1 , 12
            mbase = mbase + 10000
            if ( mon==1 .or. mon==3 .or. mon==5 .or. mon==7 .or.        &
               & mon==8 .or. mon==10 .or. mon==12 ) then
              nday = 31
            else if ( mon==4 .or. mon==6 .or. mon==9 .or. mon==11 ) then
              nday = 30
            else if ( mod(nyear,4)==0 ) then
              nday = 29
              if ( mod(nyear,100)==0 ) then
                nday = 28
                if ( mod(nyear,400)==0 ) nday = 29
              end if
            else
              nday = 28
            end if
            nbase = mbase
            do i = 1 , nday
              nbase = nbase + 100
              do m = 1 , 4
                nrec = nrec + 1
                if ( nrec>29824 ) go to 99999
                if ( m==1 ) then
                  mdate(nrec) = nbase
                else if ( m==2 ) then
                  mdate(nrec) = nbase + 6
                else if ( m==3 ) then
                  mdate(nrec) = nbase + 12
                else
                  mdate(nrec) = nbase + 18
                end if
              end do
            end do
          end do
        end do
99999   continue
      end subroutine initdate

      subroutine initdate_eh50(ssttyp)
        implicit none
!
!       Dummy arguments
!
        character(5) :: ssttyp
        intent (in) ssttyp
!
!       Local variables
!
        integer :: i , m , mbase , mon , nbase , nday , nrec , nyear
!
        nrec = 0
        if ( ssttyp=='EH5RF' ) then
          do nyear = 1941 , 2000
            mbase = nyear*1000000
            do mon = 1 , 12
              mbase = mbase + 10000
              if ( mon==1 .or. mon==3 .or. mon==5 .or. mon==7 .or.      &
                 & mon==8 .or. mon==10 .or. mon==12 ) then
                nday = 31
              else if ( mon==4 .or. mon==6 .or. mon==9 .or. mon==11 ) then
                nday = 30
              else if ( mod(nyear,4)==0 ) then
                nday = 29
                if ( mod(nyear,100)==0 ) then
                  nday = 28
                  if ( mod(nyear,400)==0 ) nday = 29
                end if
              else
                nday = 28
              end if
              nbase = mbase
              do i = 1 , nday
                nbase = nbase + 100
                do m = 1 , 4
                  nrec = nrec + 1
                  if ( m==1 ) then
                    mdate(nrec) = nbase
                  else if ( m==2 ) then
                    mdate(nrec) = nbase + 6
                  else if ( m==3 ) then
                    mdate(nrec) = nbase + 12
                  else
                    mdate(nrec) = nbase + 18
                  end if
                end do
              end do
            end do
          end do
          mdate(87661) = 2001010100
        else if ( ssttyp=='EH5A2' .or. ssttyp=='EH5B1' .or.             &
                 &ssttyp=='EHA1B' ) then
          do nyear = 2001 , 2100
            mbase = nyear*1000000
            do mon = 1 , 12
              mbase = mbase + 10000
              if ( mon==1 .or. mon==3 .or. mon==5 .or. mon==7 .or.      &
                 & mon==8 .or. mon==10 .or. mon==12 ) then
                nday = 31
              else if ( mon==4 .or. mon==6 .or. mon==9 .or. mon==11 ) then
                nday = 30
              else if ( mod(nyear,4)==0 ) then
                nday = 29
                if ( mod(nyear,100)==0 ) then
                  nday = 28
                  if ( mod(nyear,400)==0 ) nday = 29
                end if
              else
                nday = 28
              end if
              nbase = mbase
              do i = 1 , nday
                nbase = nbase + 100
                do m = 1 , 4
                  nrec = nrec + 1
                  if ( m==1 ) then
                    mdate(nrec) = nbase
                  else if ( m==2 ) then
                    mdate(nrec) = nbase + 6
                  else if ( m==3 ) then
                    mdate(nrec) = nbase + 12
                  else
                    mdate(nrec) = nbase + 18
                  end if
                end do
              end do
            end do
          end do
        else
        end if
      end subroutine initdate_eh50

      subroutine finddate_eh50(npos,idate)
        implicit none
!
!       Dummy arguments
!
        integer :: idate , npos
        intent (in) idate
        intent (out) npos
!
!       Local variables
!
        integer :: i
!
        i = 0
 100    continue
        i = i + 1
        if ( mdate(i)==idate ) then
          npos = i
          go to 99999
        end if
        if ( i<=146096 ) go to 100
        write (*,*) 'ERROR IN FINDDATE'
        stop
99999   continue
      end subroutine finddate_eh50

      end module mod_datenum
