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

      module mod_date
      implicit none
 
      integer :: ntimax
      integer :: nnnnnn , nnnchk , nnnend , nstart , nstrt0 , nnbase
      integer :: ldatez , lyear , lmonth , lday , lhour
      integer :: jyear , jyear0 , jyearr , ntime

      integer :: mdate , mdate0 , moutdate
      integer :: nmonth , nyear

      integer :: mmrec , ndate0 , ndate1
      integer :: idate0 , idate1 , idate2 , idatee , idates , idatex

      integer :: julday , julian
      real(8) :: declin , dectim , deltmx , gmt
      real(8) :: xtime
      integer :: ktau , ktaur

      real(8) :: calday , dtime , twodt
      logical :: doabsems , dolw , dosw
      integer :: mbdate , mbsec , mcdate , mcsec , mdbase , mdcur ,     &
               & msbase , mscur , nelapse , nestep , nnbdat , nnbsec ,  &
               & nndbas , nnsbas , nrstrt , nstep , nstepr , nstop

      contains

      subroutine initdate
      use mod_dynparam
      use mod_param1
      use mod_message
      implicit none
!
! Local variables
!
      integer :: i , m , mbase , mon , nbase , nday , nrec , nyear
!
      if ( ibdyfrq.eq.6 ) then
        nrec = 0
        do nyear = 1948 , 2145
          mbase = nyear*1000000
          do mon = 1 , 12
            mbase = mbase + 10000
            if ( mon.eq.1 .or. mon.eq.3 .or. mon.eq.5 .or. mon.eq.7 .or.&
               & mon.eq.8 .or. mon.eq.10 .or. mon.eq.12 ) then
              nday = 31
            else if ( mon.eq.4 .or. mon.eq.6 .or. mon.eq.9 .or.         &
                    & mon.eq.11 ) then
              nday = 30
            else
              nday = 28
              if ( mod(nyear,400).eq.0 .or.                             &
               & ( mod(nyear,4).eq.0 .and. mod(nyear,100).ne.0 ) )      &
               &  nday = nday + 1
            end if
            nbase = mbase
            do i = 1 , nday
              nbase = nbase + 100
              do m = 1 , 4
                nrec = nrec + 1
                if ( m.eq.1 ) then
                  mdatez(nrec) = nbase
                else if ( m.eq.2 ) then
                  mdatez(nrec) = nbase + 6
                else if ( m.eq.3 ) then
                  mdatez(nrec) = nbase + 12
                else
                  mdatez(nrec) = nbase + 18
                end if
              end do
            end do
          end do
        end do
      else if ( ibdyfrq.eq.3 ) then
        nrec = 0
        do nyear = 1948 , 2046
          mbase = nyear*1000000
          do mon = 1 , 12
            mbase = mbase + 10000
            if ( mon.eq.1 .or. mon.eq.3 .or. mon.eq.5 .or. mon.eq.7 .or.&
               & mon.eq.8 .or. mon.eq.10 .or. mon.eq.12 ) then
              nday = 31
            else if ( mon.eq.4 .or. mon.eq.6 .or. mon.eq.9 .or.         &
                    & mon.eq.11 ) then
              nday = 30
            else
              nday = 28
              if ( mod(nyear,400).eq.0 .or.                             &
               & ( mod(nyear,4).eq.0 .and. mod(nyear,100).ne.0 ) )      &
               &  nday = nday + 1
            end if
            nbase = mbase
            do i = 1 , nday
              nbase = nbase + 100
              do m = 1 , 8
                nrec = nrec + 1
                if ( m.eq.1 ) then
                  mdatez(nrec) = nbase
                else if ( m.eq.2 ) then
                  mdatez(nrec) = nbase + 3
                else if ( m.eq.3 ) then
                  mdatez(nrec) = nbase + 6
                else if ( m.eq.4 ) then
                  mdatez(nrec) = nbase + 9
                else if ( m.eq.5 ) then
                  mdatez(nrec) = nbase + 12
                else if ( m.eq.6 ) then
                  mdatez(nrec) = nbase + 15
                else if ( m.eq.7 ) then
                  mdatez(nrec) = nbase + 18
                else
                  mdatez(nrec) = nbase + 21
                end if
              end do
            end do
          end do
        end do
      else if ( ibdyfrq.eq.12 ) then
        nrec = 0
        do nyear = 1948 , 2145
          mbase = nyear*1000000
          do mon = 1 , 12
            mbase = mbase + 10000
            if ( mon.eq.1 .or. mon.eq.3 .or. mon.eq.5 .or. mon.eq.7 .or.&
               & mon.eq.8 .or. mon.eq.10 .or. mon.eq.12 ) then
              nday = 31
            else if ( mon.eq.4 .or. mon.eq.6 .or. mon.eq.9 .or.         &
                    & mon.eq.11 ) then
              nday = 30
            else
              nday = 28
              if ( mod(nyear,400).eq.0 .or.                             &
               & ( mod(nyear,4).eq.0 .and. mod(nyear,100).ne.0 ) )      &
               &  nday = nday + 1
            end if
            nbase = mbase
            do i = 1 , nday
              nbase = nbase + 100
              do m = 1 , 2
                nrec = nrec + 1
                if ( m.eq.1 ) then
                  mdatez(nrec) = nbase
                else if ( m.eq.2 ) then
                  mdatez(nrec) = nbase + 12
                else
                end if
              end do
            end do
          end do
        end do
      else
        write (aline,*) 'please double check your ibdyfrq'
        call say
        write (aline,*) 'we only support ibdyfrq = 3, 6, or 12'
        call say
        write (aline,*) 'and we use ibdyfrq = 6 as the default'
        call say
        call fatal(__FILE__,__LINE__,'UNSUPPORTED BOUNDARY FREQUENCY')
      end if
      write (aline,*) 'nrec = ' , nrec
      call say
      end subroutine initdate
!
!
!
      subroutine finddate(npos,idate)
      use mod_dynparam
      use mod_message
      implicit none
!
! Dummy arguments
!
      integer :: idate , npos
      intent (in) idate
      intent (inout) npos
!
! Local variables
!
      integer :: i
!
      i = 0
      do
        i = i + 1
        if ( mdatez(i).eq.idate ) then
          npos = i
          exit
        end if
        if ( i.gt.289276 ) then
 
          write (aline,*) 'error in finddate : ' , npos , idate
          call say
          call fatal(__FILE__,__LINE__,'Date not found')
          exit
        end if
      end do
      end subroutine finddate

      end module mod_date
