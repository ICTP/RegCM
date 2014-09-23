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

module mod_sst_eh5om

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_message
  use mod_dynparam
  use mod_sst_grid
  use mod_interp
  use netcdf

  private

  public :: sst_eh5om

  integer(ik4) , parameter :: a1  = 1941010106
  integer(ik4) , parameter :: a2  = 1961123118
  integer(ik4) , parameter :: a3  = 1962010100
  integer(ik4) , parameter :: a4  = 1993123118
  integer(ik4) , parameter :: a5  = 1994010100
  integer(ik4) , parameter :: a6  = 2001010100
  integer(ik4) , parameter :: a7  = 2001010100
  integer(ik4) , parameter :: a8  = 2029123118
  integer(ik4) , parameter :: a9  = 2030010100
  integer(ik4) , parameter :: a10 = 2061123118
  integer(ik4) , parameter :: a11 = 2062010100
  integer(ik4) , parameter :: a12 = 2093123118
  integer(ik4) , parameter :: a13 = 2094010100
  integer(ik4) , parameter :: a14 = 2100123118

  type(rcm_time_and_date) :: ieh5ostart

  contains
  !
  ! Comments on dataset sources and location:
  !
  ! EH5OM    EH5OM_SST is from the coupled model run by EC-MPI
  !          6 hourly frequncy, 1.875x1.875
  !          for 'RF'       run , from 1941 to 2000,
  !          for 'A2'/'B2'/'A1B', from 2001 to 2100.
  !
  !      ML = 1 is   0.0; ML = 2 is   1.875; => ML = 192 is 358.125E
  !      NL = 1 is  90.0; ML = 2 is  88.75 ; => ML = 96 is -89.0625
  !
  !
  subroutine sst_eh5om
    implicit none
    integer(ik4) , parameter :: ilon = 192 , jlat = 96
    integer(ik4) , parameter :: idtbc = 6
    integer(ik4) :: it_base
    integer(2) , dimension(ilon,jlat) :: ivar
    real(rk8) :: offset , xscale
    real(rk8) , dimension(ilon,jlat) :: sst
    type(rcm_time_and_date) :: idate
    integer(ik4) :: g1 , g2 , i1
    type(rcm_time_interval) :: tdiff , itbc
    integer(ik4) :: ieh5orec , nsteps
    integer(ik4) :: i , it , j
    integer(ik8) :: ilenrec
    real(rk8) , dimension(jlat) :: lati
    real(rk8) , dimension(ilon) :: loni
    character(len=256) :: fname
    character(3) :: code
    character(3) :: numcode
    logical :: there

    it_base = 0

    itbc = rcm_time_interval(idtbc,uhrs)

    g1  = toint10(globidate1)
    g2  = toint10(globidate2)

    there = .false.
    fname = 'unknown'
    if ( ssttyp(1:2) == 'EH' ) then
      if ( .not. date_in_scenario(globidate1,3) ) then
        if ( (g1 >= a1 .and. g1 <= a2) .or. &
             (g2 >= a1 .and. g2 <= a2) ) then
          fname = trim(inpglob)//'/SST/SST_20C_3_1941010106_1961123118'
        else if ( (g1 >= a3 .and. g1 <= a4) .or. &
                  (g2 >= a3 .and. g2 <= a4) ) then
          fname = trim(inpglob)//'/SST/SST_20C_3_1962010100_1993123118'
        else if ( (g1 >= a5 .and. g1 <= a6) .or. &
                  (g2 >= a5 .and. g2 <= a6) ) then
          fname = trim(inpglob)//'/SST/SST_20C_3_1994010100_2001010100'
        end if
        inquire (file=fname, exist=there)
        if ( .not.there ) then
          call die('sst_eh5om', trim(fname)//' is not available'// &
                 ' under '//trim(inpglob)//'/SST/',1)
        end if
      else
        if ( ssttyp(3:3) == '5' ) then
          code = ssttyp(4:5)
          numcode = '_1_'
        else
          code = ssttyp(3:5)
          numcode = '_3_'
        end if
        if ( (g1 >= a7 .and. g1 <= a8) .or. &
             (g2 >= a7 .and. g2 <= a8) ) then
          fname = trim(inpglob)//'/SST/SST_'//trim(code)// &
                       numcode//'2001010100_2029123118'
        else if ( (g1 >= a9 .and. g1 <= a10) .or. &
                  (g2 >= a9 .and. g2 <= a10) ) then
          fname = trim(inpglob)//'/SST/SST_'//trim(code)// &
                       numcode//'2030010100_2061123118'
        else if ( (g1 >= a11 .and. g1 <= a12) .or. &
                  (g2 >= a11 .and. g2 <= a12) ) then
          fname = trim(inpglob)//'/SST/SST_'//trim(code)// &
                       numcode//'2062010100_2093123118'
        else if ( (g1 >= a13 .and. g1 <= a14) .or. &
                  (g2 >= a13 .and. g2 <= a14) ) then
          fname = trim(inpglob)//'/SST/SST_'//trim(code)// &
                       numcode//'2094010100_2100123118'
        end if
        inquire (file=fname, exist=there)
        if ( .not.there ) then
          call die('sst_eh5om', trim(fname)//' is not available'// &
                 ' under '//trim(inpglob)//'/SST/',1)
        end if
      end if
    else
      write (stderr,*) 'PLEASE SET right SSTTYP in regcm.in'
      write (stderr,*) 'Supported types are EH5A2 EH5B1 EHA1B'
      call die('sst_eh5om')
    end if

    tdiff = globidate2 - globidate1
    nsteps = idnint(tohours(tdiff))/idtbc + 1
    write (stdout,*) 'GLOBIDATE1 : ' , tochar(globidate1)
    write (stdout,*) 'GLOBIDATE2 : ' , tochar(globidate2)
    write (stdout,*) 'NSTEPS     : ' , nsteps

    call open_sstfile(globidate1)

    do i = 1 , ilon
      loni(i) = float(i-1)*1.875
    end do
    do j = 1 , jlat
      lati(j) = -89.0625 + 1.875*float(j-1)
    end do

    idate = globidate1

    ilenrec = 0
    inquire(iolength=ilenrec) offset , xscale , ivar

    do it = 1 , nsteps
      i1 = toint10(idate)
      call getname(i1,fname,it_base)
      open (11,file=fname, form='unformatted',recl=ilenrec, &
        access='direct', action='read',status='old')
      tdiff = idate-ieh5ostart
      ieh5orec = idnint(tohours(tdiff))/ibdyfrq+1
      read (11,rec=ieh5orec-it_base) offset , xscale , ivar
      do j = 1 , jlat
        do i = 1 , ilon
          sst(i,j) = real(dble(ivar(i,jlat+1-j))*xscale + offset)
        end do
      end do
      call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,jx,iy,1)
      call writerec(idate)
      write (stdout,*) 'WRITING OUT SST DATA:' , tochar(idate)
      close(11)
      idate = idate + itbc
    end do
  end subroutine sst_eh5om

  subroutine getname(i1,fname,it_base)
    implicit none
    integer(ik4) , intent(in) :: i1
    character(len=256) , intent(out) :: fname
    integer(ik4) , intent(out) :: it_base
    character(3) :: code
    character(3) :: numcode
    if ( ssttyp(3:3) == '5' ) then
      code = ssttyp(4:5)
      numcode = '_1_'
    else
      code = ssttyp(3:5)
      numcode = '_3_'
    end if
    if ( (i1 >= a1 .and. i1 <= a2) ) then
      fname = trim(inpglob)//'/SST/SST_20C_3_1941010106_1961123118'
      it_base = 0
      ieh5ostart = 1941010106
    else if ( (i1 >= a3 .and. i1 <= a4) ) then
      fname = trim(inpglob)//'/SST/SST_20C_3_1962010100_1993123118'
      it_base = 30679
      ieh5ostart = 1941010106
    else if ( (i1 >= a5 .and. i1 <= a6) ) then
      fname = trim(inpglob)//'/SST/SST_20C_3_1994010100_2001010100'
      it_base = 30679 + 46752
      ieh5ostart = 1941010106
    else if ( (i1 >= a7 .and. i1 <= a8) ) then
      fname = trim(inpglob)//'/SST/SST_'//trim(code)// &
              numcode//'2001010100_2029123118'
      it_base = 0
      ieh5ostart = 2001010100
    else if ( (i1 >= a9 .and. i1 <= a10) )then
      fname = trim(inpglob)//'/SST/SST_'//trim(code)// &
              numcode//'2030010100_2061123118'
      it_base = 42368
      ieh5ostart = 2001010100
    else if ( (i1 >= a11 .and. i1 <= a12) ) then
      fname = trim(inpglob)//'/SST/SST_'//trim(code)// &
              numcode//'2062010100_2093123118'
      it_base = 42368 + 46752
      ieh5ostart = 2001010100
    else if ( (i1 >= a13 .and. i1 <= a14) ) then
      fname = trim(inpglob)//'/SST/SST_'//trim(code)// &
              numcode//'2094010100_2100123118'
      it_base = 42368 + 46752*2
      ieh5ostart = 2001010100
    end if
  end subroutine getname

end module mod_sst_eh5om

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
