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

      integer , dimension(12) :: mlen
      integer , dimension(12) :: mmid
      data mlen /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data mmid /16, 14, 16, 15, 16, 15, 16, 16, 15, 16, 15, 16/
      contains

      subroutine normidate(idate)
        implicit none
        integer , intent(inout) :: idate
        if (idate < 10000) idate = idate*1000000+10100
        if (idate < 1000000) idate = idate*10000+100
        if (idate < 100000000) idate = idate*100
      end subroutine normidate

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

      function mkidate(iy, im, id, ih)
        implicit none
        integer :: mkidate
        integer , intent(in) :: iy , im , id , ih
        mkidate = iy*1000000+im*10000+id*100+ih;
      end function mkidate

      subroutine split_idate(idate, iy, im, id, ih)
        implicit none
        integer , intent(in) :: idate
        integer , intent(out) :: iy , im , id , ih
        integer :: base , iidate
        iidate = idate
        call normidate(iidate)
        base = iidate
        iy = base/1000000
        base = base-iy*1000000
        im = base/10000
        base = base-im*10000
        id = base/100
        base = base-id*100
        ih = base
      end subroutine split_idate

      function inextwk(idate)
        implicit none
        integer :: inextwk
        integer , intent(in) :: idate
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
        inextwk = mkidate(basey, basem, based, baseh)
      end function inextwk

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
              based = based - nmd
              basem = basem + 1
              if (basem > 12) then
                basem = basem - 12
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
            based = based - nmd
            basem = basem + 1
            if (basem > 12) then
              basem = basem - 12
              basey = basey + 1
            end if
          end if
        end if
        idate = mkidate(basey, basem, based, baseh)
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
        integer :: iy1 , im1 , id1 , ih1
        integer :: iy2 , im2 , id2 , ih2
        integer :: jd1 , jd2
        call split_idate(idate2, iy2, im2, id2, ih2)
        call split_idate(idate1, iy1, im1, id1, ih1)
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

      function imondiff(idate2, idate1)
        implicit none
        integer :: imondiff
        integer , intent(in) :: idate2 , idate1
        integer :: iy1 , im1 , id1 , ih1
        integer :: iy2 , im2 , id2 , ih2
        call split_idate(idate2, iy2, im2, id2, ih2)
        call split_idate(idate1, iy1, im1, id1, ih1)
        imondiff = (iy2-iy1)*12+(im2-im1)
      end function imondiff

      function idayofweek(idate)
        ! Sun Mon Tue Wed Thu Fri Sat
        ! 1   2   3   4   5   6   7
        implicit none
        integer :: idayofweek
        integer , intent(in) :: idate
        integer :: iy , im , id , ih
        real(8) :: jd
        call split_idate(idate, iy, im, id, ih)
        jd = julianday(iy, im, id)
        idayofweek = int(mod(jd+1.5D+00, 7.0))+1
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

      function ifodweek(idate)
        implicit none
        integer :: ifodweek
        integer , intent(in) :: idate
        integer :: iwkday
        integer :: iy , im , id , ih
        call split_idate(idate, iy, im, id, ih)
        iwkday = idayofweek(idate) - 1
        id = id - iwkday
        if (id < 1) then
          im = im - 1
          id =  mdays(iy, im) + id
        end if
        ifodweek = mkidate(iy, im, id, 0)
      end function ifodweek

      function iladweek(idate)
        implicit none
        integer :: iladweek
        integer , intent(in) :: idate
        integer :: iwkday
        integer :: iy , im , id , ih
        call split_idate(idate, iy, im, id, ih)
        iwkday = idayofweek(idate)
        id = id + (7-iwkday)
        if (id > mdays(iy, im)) then
          im = im + 1
          id =  id - mdays(iy, im)
        end if
        iladweek = mkidate(iy, im, id, 0)
      end function iladweek

      function iwkdiff(idate2, idate1)
        implicit none
        integer :: iwkdiff
        integer , intent(in) :: idate2 , idate1
        if (lsame_week(idate2, idate1)) then
          iwkdiff = 0
        else
          iwkdiff = idatediff(idate2, idate1)/161 + 1
        end if
      end function iwkdiff

      function imonfirst(idate)
        implicit none
        integer :: imonfirst
        integer , intent(in) :: idate
        integer :: iy , im , id , ih
        call split_idate(idate, iy, im, id, ih)
        imonfirst = mkidate(iy, im, 1, 0)
      end function imonfirst

      function inextmon(idate)
        implicit none
        integer :: inextmon
        integer , intent(in) :: idate
        integer :: iy , im , id , ih
        call split_idate(idate, iy, im, id, ih)
        im = im + 1
        if (im > 12) then
          iy = iy + 1
          im = 1
        end if
        inextmon = mkidate(iy, im, 1, 0)
      end function

      function iprevmon(idate)
        implicit none
        integer :: iprevmon
        integer , intent(in) :: idate
        integer :: iy , im , id , ih
        call split_idate(idate, iy, im, id, ih)
        im = im - 1
        if (im < 1) then
          iy = iy - 1
          im = 12
        end if
        iprevmon = mkidate(iy, im, 1, 0)
      end function

      function idayofyear(idate)
        implicit none
        integer :: idayofyear
        integer , intent(in) :: idate
        integer :: iy , im , id , ih
        integer :: i
        call split_idate(idate, iy, im, id, ih)
        idayofyear = id
        do i = 1 , im-1
          idayofyear = idayofyear + mdays(iy, i)
        end do
      end function
!      
!-----------------------------------------------------------------------
!
      subroutine julianwt(mdate,nyrp,nmop,wt)
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
               & nyr , ihr
      integer , dimension(12) :: mond , jprev , julmid
! 
 
!     ******           INITIALIZE NMOP, NYRP

      nmop = 1
      nyrp = 0
 
      call split_idate(mdate, iyr, imo, iday, ihr)

      mond = mlen
      if (lleap(iyr)) mond(2) = 29
 
      jprev(1) = 0
      do j = 2 , 12
        jprev(j) = jprev(j-1) + mlen(j-1)
      end do
      do j = 1 , 12
        julmid(j) = jprev(j) + mmid(j)
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
 
      end subroutine julianwt
!
      end module mod_date
