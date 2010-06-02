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
      integer , dimension(427+1097) :: wkday

      integer , dimension(12) :: mlen
      data mlen /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      contains

      function lleap(iyear)
        implicit none
        logical :: lleap
        integer , intent(in) :: iyear
        if ( mod(iyear,400).eq.0 .or.                                   &
          & ( mod(iyear,4).eq.0 .and. mod(iyear,100).ne.0 ) ) then
          lleap = .true.
        else
          lleap = .false.
        end if
      end function lleap

      function mdays(iyear, imon)
        implicit none
        integer :: mdays
        integer , intent(in) :: iyear , imon
        if (imon /= 2) then
          mdays = mlen(imon)
        else
          mdays = mlen(2)
          if (lleap(iyear)) then
            mdays = mdays + 1
          end if
        end if
      end function mdays

      function iidate(iy, im, id, ih)
        implicit none
        integer :: iidate
        integer , intent(in) :: iy , im , id , ih
        iidate = iy*1000000+im*10000+id*100+ih;
      end function iidate

      subroutine split_idate(idate, iy, im, id, ih)
        implicit none
        integer , intent(in) :: idate
        integer , intent(out) :: iy , im , id , ih
        integer :: base
        base = idate
        iy = base/1000000
        base = base-iy*1000000
        im = base/10000
        base = base-im*10000
        id = base/100
        base = base-id*100
        ih = base
      end subroutine split_idate

      subroutine addweek(idate)
        implicit none
        integer , intent(inout) :: idate
         integer :: nmd , basey , basem , based , baseh
        call split_idate(idate, basey, basem, based, baseh)
        based = based + 7
        nmd = mdays(basey, basem)
        if (based > nmd) then
          based = 1
          basem = basem + 1
          if (basem > 12) then
            basem = 1
            basey = basey + 1
          end if
        end if
        idate = iidate(basey, basem, based, baseh)
      end subroutine addweek

      subroutine addhours(idate, ihours)
        implicit none
        integer , intent(in) :: ihours
        integer , intent(inout) :: idate
        integer :: basey , basem , based , baseh , nmd
        integer :: nsum , isum , ilast

        call split_idate(idate, basey, basem, based, baseh)
        nsum  = ihours / 23
        ilast = ihours - nsum*23;
        do isum = 1 , nsum
          baseh = baseh + 23
          if (baseh > 23) then
            based = based + 1
            baseh = baseh - 24
            nmd = mdays(basey, basem)
            if (based > nmd) then
              based = 1
              basem = basem + 1
              if (basem > 12) then
                basem = 1
                basey = basey + 1
              end if
            end if
          end if
        end do
        baseh = baseh + ilast
        if (baseh > 23) then
          based = based + 1
          baseh = 0
          nmd = mdays(basey, basem)
          if (based > nmd) then
            based = 1
            basem = basem + 1
            if (basem > 12) then
              basem = 1
              basey = basey + 1
            end if
          end if
        end if
        idate = iidate(basey, basem, based, baseh)
      end subroutine addhours

      function lcaltype(iy, im, id)
        implicit none
        logical :: lcaltype
        integer :: icaltype
        integer , intent(in) :: iy , im , id
        ! Standard Julian/Gregorian switch
        ! Return true  if before 1582-10-04
        !        false if after  1582-10-15
        icaltype = 0
        if (iy < 1582) then
          icaltype = 1
        else if (iy == 1582) then
          if (im < 10) then
            icaltype = 1
          else if (im == 10) then
            if (id <= 4) then
              icaltype = 1
            end if
          end if
        end if
        if (iy > 1582) then
          icaltype = 2
        else if (iy == 1582) then
          if (im > 10) then
            icaltype = 2
          else if (im == 10) then
            if (id >= 15) then
              icaltype = 2
            end if
          end if
        end if
        if (icaltype == 0) then
          write (6, *) 'year  = ', iy
          write (6, *) 'month = ', im
          write (6, *) 'day   = ', id
          write (6, *) 'Day non existent, inside Julian/Gregorian jump'
          stop
        end if
        if (icaltype == 1) then
          lcaltype = .false.
        else
          lcaltype = .true.
        end if
      end function lcaltype

      function julianday(iy, im, id)
        implicit none
        real(8) :: julianday
        integer , intent(in) :: iy , im , id
        integer :: ia , ib , iiy , iim
        iiy = iy
        iim = im
        if (iim <= 2) then
          iiy = iiy - 1
          iim = iim + 12
        end if
        if (lcaltype(iy,im,id)) then
          ia = iiy/100
          ib = 2 - ia + ia / 4
        else
          ib = 0
        end if
        julianday = int(365.25D+00*(iiy+4716)) +                        &
          &         int(30.6001D+00*(iim+1))   +                        &
          &         id + ib - 1524.5D+00
      end function julianday

      function idatediff(idate2, idate1)
        !
        ! Returns number of hours between to integer dates in the format
        !
        !                      YYYYMMDDHH
        !
        ! If just YYYY is passed, difference is
        !                      from jan, 01 00:00:00 UTC
        ! If just YYYYMM is passed, difference is
        !                      from first day of month, 00:00:00 UTC
        ! If just YYYYMMDD is passed, difference is
        !                      from 00:00:00 UTC
        !
        implicit none
        integer :: idatediff
        integer , intent(in) :: idate2 , idate1
        integer :: iidate1 , iidate2
        integer :: iy1 , im1 , id1 , ih1
        integer :: iy2 , im2 , id2 , ih2
        integer :: jd1 , jd2
        iidate1 = idate1
        iidate2 = idate2
        if (iidate1 < 10000) iidate1 = iidate1*1000000+10100
        if (iidate2 < 10000) iidate2 = iidate2*1000000+10100
        if (iidate1 < 1000000) iidate1 = iidate1*10000+100
        if (iidate2 < 1000000) iidate2 = iidate2*10000+100
        if (iidate1 < 100000000) iidate1 = iidate1*100
        if (iidate2 < 100000000) iidate2 = iidate2*100
        call split_idate(iidate2, iy2, im2, id2, ih2)
        call split_idate(iidate1, iy1, im1, id1, ih1)
        jd2 = julianday(iy2, im2, id2)
        jd1 = julianday(iy1, im1, id1)
        idatediff = (jd2-jd1)*24+(ih2-ih1)
      end function idatediff

      function lsame_month(idate1, idate2)
        implicit none
        logical :: lsame_month
        integer , intent(in) :: idate1 , idate2
        lsame_month = .false.
        if (abs(idate1-idate2) < 10000) lsame_month = .true.
      end function

      function idayofweek(idate)
        implicit none
        integer :: idayofweek
        integer , intent(in) :: idate
        integer :: iy , im , id , ih , jd
        call split_idate(idate, iy, im, id, ih)
        jd = julianday(iy, im, id)
        idayofweek = ceiling(mod(jd+1.5D+00, 7.0))+1
      end function idayofweek

      function lsame_week(idate1, idate2)
        implicit none
        logical :: lsame_week
        integer , intent(in) :: idate1 , idate2
        integer :: iidate1 , iidate2
        integer :: idatewk1, ieofweek, idist
        lsame_week = .false.
        iidate1 = idate1/100*100
        iidate2 = idate2/100*100
        idist = idatediff(iidate2, iidate1)
        if (idist == 0) then
          lsame_week = .true.
        else
          idatewk1 = idayofweek(iidate1)
          ieofweek = 7-idatewk1
          if (idist > 0 .and. idist/24 <=  ieofweek) lsame_week = .true.
          if (idist < 0 .and. idist/24 >= -idatewk1) lsame_week = .true.
        end if
      end function lsame_week
      
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
      do i = 429 , 427 + 1097
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
