!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_datenum

      implicit none

      integer , dimension(300000) :: mdate

      contains

      subroutine initdate_era
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
      end subroutine initdate_era

      subroutine finddate_era(npos,idate)
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
      end subroutine finddate_era

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
!
!-----------------------------------------------------------------------
!
      subroutine initdate_ccsm
        implicit none
!
!       Local variables
!
        integer :: i , m , mbase , mon , nbase , nday , nrec , nyear
!
        nrec = 0
        do nyear = 1948 , 2045
          mbase = nyear*1000000
          do mon = 1 , 12
            mbase = mbase + 10000
            if ( mon==1 .or. mon==3 .or. mon==5 .or. mon==7 .or.        &
               & mon==8 .or. mon==10 .or. mon==12 ) then
              nday = 31
            else if ( mon==4 .or. mon==6 .or. mon==9 .or. mon==11 ) then
              nday = 30
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
        write (*,*) 'nrec = ' , nrec
      end subroutine initdate_ccsm

      subroutine initdate_icbc
        implicit none
!
! Local variables
!
        integer :: i , m , mbase , mon , nbase , nday , nrec , nyear
!
        nrec = 0
        do nyear = 1941 , 2145
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
              if ( mod(nyear,100)==0 ) nday = nday - 1
              if ( mod(nyear,400)==0 ) nday = nday + 1
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
        write (*,*) 'NREC = ' , nrec
      end subroutine initdate_icbc

      subroutine finddate_icbc(npos,idate)
        implicit none
!
! Dummy arguments
!
        integer :: idate , npos
        intent (in) idate
        intent (out) npos
!
! Local variables
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
        if ( i<=299500 ) go to 100
        write (*,*) 'ERROR IN FINDDATE'
        stop
99999   continue
      end subroutine finddate_icbc

      end module mod_datenum
