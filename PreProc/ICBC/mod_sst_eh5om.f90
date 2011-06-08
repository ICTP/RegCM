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

  use m_realkinds
  use m_die
  use m_stdio
  use m_zeit
  use netcdf
  use mod_dynparam
  use mod_sst_grid
  use mod_date
  use mod_interp

  contains

  subroutine sst_eh5om

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comments on dataset sources and location:                          !
!                                                                    !
! EH5OM    EH5OM_SST is from the coupled model run by EC-MPI         !
!          6 hourly frequncy, 1.875x1.875                            !
!          for 'RF'       run , from 1941 to 2000,                   !
!          for 'A2'/'B2'/'A1B', from 2001 to 2100.                   !
!                                                                    !
!      ML = 1 is   0.0; ML = 2 is   1.875; => ML = 192 is 358.125E   !
!      NL = 1 is  90.0; ML = 2 is  88.75 ; => ML = 96 is -89.0625    !
!                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none
!
  integer , parameter :: ilon = 192 , jlat = 96
  integer , parameter :: idtbc = 6
!
  integer :: it_base
  integer(2) , dimension(ilon,jlat) :: ivar
  real(dp) :: offset , xscale
  real(sp) , dimension(ilon,jlat) :: sst
  type(rcm_time_and_date) :: idate , ieh5ostart
  type(rcm_time_and_date) :: a1 , a2 , a3 , a4 , a5 , a6 , a7 , a8 , a9 , &
                             a10  , a11 , a12 , a13 , a14
  type(rcm_time_interval) :: tdiff
  integer :: ieh5orec , nsteps
  integer :: i , it , j , nday , nhour , nmo , nyear
  real(sp) , dimension(jlat) :: lati
  real(sp) , dimension(ilon) :: loni
  logical :: there
!
  call zeit_ci('sst_eh5om')
!
  it_base = 0

  a1 = rcm_time_and_date(gregorian,1941,1,1,6,0,0)
  a2 = rcm_time_and_date(gregorian,1961,12,31,18,0,0)
  a3 = rcm_time_and_date(gregorian,1962,1,1,0,0,0)
  a4 = rcm_time_and_date(gregorian,1993,12,31,18,0,0)
  a5 = rcm_time_and_date(gregorian,1994,1,1,0,0,0)
  a6 = rcm_time_and_date(gregorian,2001,1,1,0,0,0)
  a7 = rcm_time_and_date(gregorian,2001,1,1,0,0,0)
  a8 = rcm_time_and_date(gregorian,2029,12,31,18,0,0)
  a9 = rcm_time_and_date(gregorian,2030,1,1,0,0,0)
  a10 = rcm_time_and_date(gregorian,2061,12,31,18,0,0)
  a11 = rcm_time_and_date(gregorian,2062,1,1,0,0,0)
  a12 = rcm_time_and_date(gregorian,2093,12,31,18,0,0)
  a13 = rcm_time_and_date(gregorian,2094,1,1,0,0,0)
  a14 = rcm_time_and_date(gregorian,2100,12,31,18,0,0)

  if ( ssttyp == 'EH5RF' ) then
    there = .false.
    if ( (globidate1 >= a1 .and. globidate1 <= a2) .or. &
         (globidate2 >= a1 .and. globidate2 <= a2) ) then
      inquire (file=trim(inpglob)//'/SST/SST_20C_3_1941010106_1961123118', &
               exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_20C_3_1941010106_1961123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (globidate1 >= a3 .and. globidate1 <= a4) .or. &
         (globidate2 >= a3 .and. globidate2 <= a4) ) then
      inquire (file=trim(inpglob)//'/SST/SST_20C_3_1962010100_1993123118', &
               exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_20C_3_1962010100_1993123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (globidate1 >= a5 .and. globidate1 <= a6) .or. &
         (globidate2 >= a5 .and. globidate2 <= a6) ) then
      inquire (file=trim(inpglob)//'/SST/SST_20C_3_1994010100_2001010100', &
               exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_20C_3_1994010100_2001010100 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( .not.there ) then
      call die('sst_eh5om', &
               'EH5RF SST is just available from 1941010106 to'// &
               ' 2001010100',1)
    end if
  else if ( ssttyp == 'EH5A2' ) then
    there = .false.
    if ( (globidate1 >= a7 .and. globidate1 <= a8) .or. &
         (globidate2 >= a7 .and. globidate2 <= a8) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A2_1_2001010100_2029123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A2_1_2001010100_2029123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (globidate1 >= a9 .and. globidate1 <= a10) .or. &
         (globidate2 >= a9 .and. globidate2 <= a10) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A2_1_2030010100_2061123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A2_1_2030010100_2061123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (globidate1 >= a11 .and. globidate1 <= a12) .or. &
         (globidate2 >= a11 .and. globidate2 <= a12) ) then
      inquire (file=trim(inpglob)//                                 &
          '/SST/SST_A2_1_2062010100_2093123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A2_1_2062010100_2093123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (globidate1 >= a13 .and. globidate1 <= a14) .or. &
         (globidate2 >= a13 .and. globidate2 <= a14) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A2_1_2094010100_2100123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A2_1_2094010100_2100123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( .not.there ) then
      call die('sst_eh5om', &
               'EH5A2 SST is just available from 2001010100 to'// &
               ' 2100123118',1)
    end if
  else if ( ssttyp == 'EH5B1' ) then
    there = .false.
    if ( (globidate1 >= a7 .and. globidate1 <= a8) .or. &
         (globidate2 >= a7 .and. globidate2 <= a8) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_B1_1_2001010100_2029123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_B1_1_2001010100_2029123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (globidate1 >= a9 .and. globidate1 <= a10) .or. &
         (globidate2 >= a9 .and. globidate2 <= a10) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_B1_1_2030010100_2061123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_B1_1_2030010100_2061123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (globidate1 >= a11 .and. globidate1 <= a12) .or. &
         (globidate2 >= a11 .and. globidate2 <= a12) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_B1_1_2062010100_2093123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_B1_1_2062010100_2093123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (globidate1 >= a13 .and. globidate1 <= a14) .or. &
         (globidate2 >= a13 .and. globidate2 <= a14) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_B1_1_2094010100_2100123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_B1_1_2094010100_2100123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( .not.there ) then
      call die('sst_eh5om', &
               'EH5B1 SST is just available from 2001010100 to'// &
               ' 2100123118',1)
    end if
  else if ( ssttyp == 'EHA1B' ) then
    there = .false.
    if ( (globidate1 >= a7 .and. globidate1 <= a8) .or. &
         (globidate2 >= a7 .and. globidate2 <= a8) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A1B_3_2001010100_2029123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A1B_3_2001010100_2029123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (globidate1 >= a9 .and. globidate1 <= a10) .or. &
         (globidate2 >= a9 .and. globidate2 <= a10) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A1B_3_2030010100_2061123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A1B_3_2030010100_2061123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (globidate1 >= a11 .and. globidate1 <= a12) .or. &
         (globidate2 >= a11 .and. globidate2 <= a12) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A1B_3_2062010100_2093123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A1B_3_2062010100_2093123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (globidate1 >= a13 .and. globidate1 <= a14) .or. &
         (globidate2 >= a13 .and. globidate2 <= a14) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A1B_3_2094010100_2100123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A1B_3_2094010100_2100123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( .not.there ) then
      call die('sst_eh5om', &
               'EHA1B SST is just available from 2001010100 to'// &
               ' 2100123118',1)
    end if
  else
    write (stderr,*) 'PLEASE SET right SSTTYP in regcm.in'
    write (stderr,*) 'Supported types are EH5RF EH5A2 EH5B1 EHA1B'
    call die('sst_eh5om')
  end if

  tdiff = globidate2 - globidate1
  nsteps = idnint(tdiff%hours())/idtbc + 1
  write (stdout,*) 'GLOBIDATE1 : ' , globidate1
  write (stdout,*) 'GLOBIDATE2 : ' , globidate2
  write (stdout,*) 'NSTEPS     : ' , nsteps

  call open_sstfile(globidate1)
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
  do i = 1 , ilon
    loni(i) = float(i-1)*1.875
  end do
  do j = 1 , jlat
    lati(j) = -89.0625 + 1.875*float(j-1)
  end do
 
!     **  REF  SST DATA, 1.875x1.1.25, AVAILABLE FROM 16/1/1959 TO
!     16/1/1991 ** A2&B2 SST DATA, 1.875x1.1.25, AVAILABLE FROM
!     16/1/2069 TO 16/1/2101
  idate = globidate1
  if ( ssttyp == 'EH5RF' ) then
    ieh5ostart = 1989010100
  else
    ieh5ostart = 2001010100
  end if
  do it = 1 , nsteps
    
    if ( ssttyp == 'EH5RF' ) then
      if ( idate >= a1 .and. idate <= a2 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_20C_3_1941010106_1961123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 1
      else if ( idate >= a3 .and. idate <= a4 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_20C_3_1962010100_1993123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 30680
      else if ( idate >= a5 .and. idate <= a6 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_20C_3_1994010100_2001010100',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 30680 + 46752
      else
      end if
    else if ( ssttyp == 'EH5A2' ) then
      if ( idate >= a7 .and. idate <= a8 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A2_1_2001010100_2029123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 0
      else if ( idate >= a9 .and. idate <= a10 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A2_1_2030010100_2061123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368
      else if ( idate >= a11 .and. idate <= a12 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A2_1_2062010100_2093123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752
      else if ( idate >= a13 .and. idate <= a14 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A2_1_2094010100_2100123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752*2
      else
      end if
    else if ( ssttyp == 'EH5B1' ) then
      if ( idate >= a7 .and. idate <= a8 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_B1_1_2001010100_2029123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 0
      else if ( idate >= a9 .and. idate <= a10 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_B1_1_2030010100_2061123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368
      else if ( idate >= a11 .and. idate <= a12 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_B1_1_2062010100_2093123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752
      else if ( idate >= a13 .and. idate <= a14 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_B1_1_2094010100_2100123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752*2
      else
      end if
    else if ( ssttyp == 'EHA1B' ) then
      if ( idate >= a7 .and. idate <= a8 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A1B_3_2001010100_2029123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 0
      else if ( idate >= a9 .and. idate <= a10 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A1B_3_2030010100_2061123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368
      else if ( idate >= a11 .and. idate <= a12 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A1B_3_2062010100_2093123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752
      else if ( idate >= a13 .and. idate <= a14 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A1B_3_2094010100_2100123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752*2
      end if
    end if

    tdiff = idate-ieh5ostart
    ieh5orec = idnint(tdiff%hours())/ibdyfrq+1
    read (11,rec=ieh5orec-it_base) offset , xscale , ivar
    do j = 1 , jlat
      do i = 1 , ilon
        sst(i,j) = real(dble(ivar(i,jlat+1-j))*xscale + offset)
        if ( sst(i,j) < 273.16 ) sst(i,j) = -9999.
      end do
    end do
 
    call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
    write (stdout,*) 'XLON,XLAT,SST = ' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
!       ******           WRITE OUT SST DATA ON MM4 GRID
    call writerec(idate,.false.)
    write (stdout,*) 'WRITING OUT SST DATA:' , nmo , nyear

    call addhours(idate, idtbc)

  end do

  call zeit_co('sst_eh5om')
 
  end subroutine sst_eh5om
!
end module mod_sst_eh5om
