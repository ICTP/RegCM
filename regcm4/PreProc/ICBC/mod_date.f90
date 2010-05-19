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

      integer , dimension(300000) :: mdate
      integer , dimension(427+1045) :: wkday

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
            else if ( mod(nyear,400).eq.0 .or.                          &
               & ( mod(nyear,4).eq.0 .and. mod(nyear,100).ne.0 ) ) then
              nday = 29
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
              else if ( mod(nyear,400).eq.0 .or.                        &
               &   ( mod(nyear,4).eq.0 .and. mod(nyear,100).ne.0 ) ) then
                nday = 29
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
              else if ( mod(nyear,400).eq.0 .or.                        &
               & ( mod(nyear,4).eq.0 .and. mod(nyear,100).ne.0 ) ) then
                nday = 29
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
            else if ( mod(nyear,400).eq.0 .or.                          &
               & ( mod(nyear,4).eq.0 .and. mod(nyear,100).ne.0 ) ) then
              nday = 29
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
!
!-----------------------------------------------------------------------
!
      subroutine headwk
      implicit none
!
! Local variables
!
      integer :: i , mday , month , myear
!
      wkday(1) = 19811029
      do i = 2 , 427
        wkday(i) = wkday(i-1) + 7
        myear = wkday(i)/10000
        month = wkday(i)/100 - myear*100
        mday = mod(wkday(i),10000) - month*100
        if ( month==1 .or. month==3 .or. month==5 .or. month==7 .or.    &
           & month==8 .or. month==10 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = month + 1
          end if
        else if ( month==12 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = 1
            myear = myear + 1
          end if
        else if ( month==4 .or. month==6 .or. month==9 .or. month==11 ) &
                & then
          if ( mday>30 ) then
            mday = mday - 30
            month = month + 1
          end if
        else if ( mod(myear,400).eq.0 .or.                              &
           & ( mod(myear,4).eq.0 .and. mod(myear,100).ne.0 ) ) then
          if ( mday>29 ) then
            mday = mday - 29
            month = month + 1
          end if
        else
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        end if
        wkday(i) = myear*10000 + month*100 + mday
      end do
!
      wkday(428) = 19891231
      do i = 429 , 427 + 1045
        wkday(i) = wkday(i-1) + 7
        myear = wkday(i)/10000
        month = wkday(i)/100 - myear*100
        mday = mod(wkday(i),10000) - month*100
        if ( month==1 .or. month==3 .or. month==5 .or. month==7 .or.    &
           & month==8 .or. month==10 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = month + 1
          end if
        else if ( month==12 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = 1
            myear = myear + 1
          end if
        else if ( month==4 .or. month==6 .or. month==9 .or. month==11 ) &
                & then
          if ( mday>30 ) then
            mday = mday - 30
            month = month + 1
          end if
        else if ( mod(myear,400).eq.0 .or.                              &
           & ( mod(myear,4).eq.0 .and. mod(myear,100).ne.0 ) ) then
          if ( mday>29 ) then
            mday = mday - 29
            month = month + 1
          end if
        else
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        end if
        wkday(i) = myear*10000 + month*100 + mday
      end do
!
      end subroutine headwk
!
!-----------------------------------------------------------------------
!
      subroutine julian(mdate,nyrp,nmop,wt)
      implicit none
!
! Dummy arguments
!
      integer :: mdate , nmop , nyrp
      real(4) :: wt
      intent (in) mdate
      intent (out) nyrp , wt
      intent (inout) nmop
!
! Local variables
!
      real(4) :: fdenom , fnumer
      integer :: idate , iday , ileap , imo , iyr , j , julday , nmo ,  &
               & nyr
      integer , dimension(12) :: jprev , julmid , lenmon , midmon
! 
      data lenmon/31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , 31 , 30 ,&
         & 31/
      data midmon/16 , 14 , 16 , 15 , 16 , 15 , 16 , 16 , 15 , 16 , 15 ,&
         & 16/
 
!     ******           INITIALIZE NMOP, NYRP

      nmop = 1
      nyrp = 0
 
      idate = mdate/100
      iyr = idate/10000
      imo = (idate-iyr*10000)/100
      iday = mod(idate,100)
 
      ileap = mod(iyr,4)
      lenmon(2) = 28
      if ( ileap==0 ) lenmon(2) = 29
 
      jprev(1) = 0
      do j = 2 , 12
        jprev(j) = jprev(j-1) + lenmon(j-1)
      end do
      do j = 1 , 12
        julmid(j) = jprev(j) + midmon(j)
      end do
      julday = iday + jprev(imo)
 
!     PRINT *, 'MDATE, IYR, IMO, IDAY, JULDAY = '
!     A       ,  MDATE, IYR, IMO, IDAY, JULDAY
 
      do nyr = 1948 , 2145  !94
        do nmo = 1 , 12
 
          if ( (nyr==iyr) .and. (julmid(nmo)>julday) ) go to 100
          if ( nyr>iyr ) go to 100
 
          nmop = nmo
          nyrp = nyr
 
        end do
      end do
 
 100  continue
      fnumer = float(julday-julmid(nmop))
      if ( fnumer<0. ) fnumer = fnumer + 365.
      fdenom = float(julmid(nmo)-julmid(nmop))
      if ( fdenom<=0. ) fdenom = fdenom + 365.
      wt = fnumer/fdenom
 
!     PRINT *, 'JULMID(NMOP), JULDAY, JULMID(NMO), WT ='
!     A       ,  JULMID(NMOP), JULDAY, JULMID(NMO), WT
 
      end subroutine julian
!
      end module mod_date
