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
  integer(ik4) , parameter :: ilon = 192 , jlat = 96
  integer(ik4) , parameter :: idtbc = 6
!
  integer(ik4) :: it_base
  integer(2) , dimension(ilon,jlat) :: ivar
  real(rk8) :: offset , xscale
  real(rk4) , dimension(ilon,jlat) :: sst
  type(rcm_time_and_date) :: idate , ieh5ostart
  integer(ik4) :: a1 , a2 , a3 , a4 , a5 , a6 , a7 , a8 , a9 , &
             a10  , a11 , a12 , a13 , a14 , g1 , g2 , i1
  type(rcm_time_interval) :: tdiff , itbc
  integer(ik4) :: ieh5orec , nsteps
  integer(ik4) :: i , it , j
  real(rk4) , dimension(jlat) :: lati
  real(rk4) , dimension(ilon) :: loni
  logical :: there
!
!
  it_base = 0

  itbc = rcm_time_interval(idtbc,uhrs)

  a1  = 1941010106
  a2  = 1961123118
  a3  = 1962010100
  a4  = 1993123118
  a5  = 1994010100
  a6  = 2001010100
  a7  = 2001010100
  a8  = 2029123118
  a9  = 2030010100
  a10 = 2061123118
  a11 = 2062010100
  a12 = 2093123118
  a13 = 2094010100
  a14 = 2100123118
  g1  = toint10(globidate1)
  g2  = toint10(globidate2)

  if ( ssttyp == 'EH5RF' ) then
    there = .false.
    if ( (g1 >= a1 .and. g1 <= a2) .or. (g2 >= a1 .and. g2 <= a2) ) then
      inquire (file=trim(inpglob)//'/SST/SST_20C_3_1941010106_1961123118', &
               exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_20C_3_1941010106_1961123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (g1 >= a3 .and. g1 <= a4) .or. (g2 >= a3 .and. g2 <= a4) ) then
      inquire (file=trim(inpglob)//'/SST/SST_20C_3_1962010100_1993123118', &
               exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_20C_3_1962010100_1993123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (g1 >= a5 .and. g1 <= a6) .or. (g2 >= a5 .and. g2 <= a6) ) then
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
    if ( (g1 >= a7 .and. g1 <= a8) .or. (g2 >= a7 .and. g2 <= a8) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A2_1_2001010100_2029123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A2_1_2001010100_2029123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (g1 >= a9 .and. g1 <= a10) .or. (g2 >= a9 .and. g2 <= a10) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A2_1_2030010100_2061123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A2_1_2030010100_2061123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (g1 >= a11 .and. g1 <= a12) .or. (g2 >= a11 .and. g2 <= a12) ) then
      inquire (file=trim(inpglob)//                                 &
          '/SST/SST_A2_1_2062010100_2093123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A2_1_2062010100_2093123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (g1 >= a13 .and. g1 <= a14) .or. (g2 >= a13 .and. g2 <= a14) ) then
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
    if ( (g1 >= a7 .and. g1 <= a8) .or. (g2 >= a7 .and. g2 <= a8) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_B1_1_2001010100_2029123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_B1_1_2001010100_2029123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (g1 >= a9 .and. g1 <= a10) .or. (g2 >= a9 .and. g2 <= a10) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_B1_1_2030010100_2061123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_B1_1_2030010100_2061123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (g1 >= a11 .and. g1 <= a12) .or. (g2 >= a11 .and. g2 <= a12) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_B1_1_2062010100_2093123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_B1_1_2062010100_2093123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (g1 >= a13 .and. g1 <= a14) .or. (g2 >= a13 .and. g2 <= a14) ) then
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
    if ( (g1 >= a7 .and. g1 <= a8) .or. (g2 >= a7 .and. g2 <= a8) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A1B_3_2001010100_2029123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A1B_3_2001010100_2029123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (g1 >= a9 .and. g1 <= a10) .or. (g2 >= a9 .and. g2 <= a10) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A1B_3_2030010100_2061123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A1B_3_2030010100_2061123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (g1 >= a11 .and. g1 <= a12) .or. (g2 >= a11 .and. g2 <= a12) ) then
      inquire (file=trim(inpglob)//                                 &
               '/SST/SST_A1B_3_2062010100_2093123118',exist=there)
      if ( .not.there ) then
        call die('sst_eh5om', &
               'SST_A1B_3_2062010100_2093123118 is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
    if ( (g1 >= a13 .and. g1 <= a14) .or. (g2 >= a13 .and. g2 <= a14) ) then
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
  nsteps = idnint(tohours(tdiff))/idtbc + 1
  write (stdout,*) 'GLOBIDATE1 : ' , tochar(globidate1)
  write (stdout,*) 'GLOBIDATE2 : ' , tochar(globidate2)
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
    ieh5ostart = 1941010106
  else
    ieh5ostart = 2001010100
  end if

  do it = 1 , nsteps
    
    i1 = toint10(idate)

    if ( ssttyp == 'EH5RF' ) then
      if ( i1 >= a1 .and. i1 <= a2 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_20C_3_1941010106_1961123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 0
      else if ( i1 >= a3 .and. i1 <= a4 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_20C_3_1962010100_1993123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 30679
      else if ( i1 >= a5 .and. i1 <= a6 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_20C_3_1994010100_2001010100',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 30679 + 46752
      else
      end if
    else if ( ssttyp == 'EH5A2' ) then
      if ( i1 >= a7 .and. i1 <= a8 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A2_1_2001010100_2029123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 0
      else if ( i1 >= a9 .and. i1 <= a10 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A2_1_2030010100_2061123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368
      else if ( i1 >= a11 .and. i1 <= a12 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A2_1_2062010100_2093123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752
      else if ( i1 >= a13 .and. i1 <= a14 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A2_1_2094010100_2100123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752*2
      else
      end if
    else if ( ssttyp == 'EH5B1' ) then
      if ( i1 >= a7 .and. i1 <= a8 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_B1_1_2001010100_2029123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 0
      else if ( i1 >= a9 .and. i1 <= a10 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_B1_1_2030010100_2061123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368
      else if ( i1 >= a11 .and. i1 <= a12 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_B1_1_2062010100_2093123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752
      else if ( i1 >= a13 .and. i1 <= a14 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_B1_1_2094010100_2100123118',            &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752*2
      else
      end if
    else if ( ssttyp == 'EHA1B' ) then
      if ( i1 >= a7 .and. i1 <= a8 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A1B_3_2001010100_2029123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 0
      else if ( i1 >= a9 .and. i1 <= a10 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A1B_3_2030010100_2061123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368
      else if ( i1 >= a11 .and. i1 <= a12 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A1B_3_2062010100_2093123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752
      else if ( i1 >= a13 .and. i1 <= a14 ) then
        open (11,file=trim(inpglob)//                           &
              '/SST/SST_A1B_3_2094010100_2100123118',           &
              form='unformatted',recl=(ilon*jlat/2+4)*ibyte,    &
              access='direct')
        it_base = 42368 + 46752*2
      end if
    end if

    tdiff = idate-ieh5ostart
    ieh5orec = idnint(tohours(tdiff))/ibdyfrq+1
    read (11,rec=ieh5orec-it_base) offset , xscale , ivar
    do j = 1 , jlat
      do i = 1 , ilon
        sst(i,j) = real(dble(ivar(i,jlat+1-j))*xscale + offset)
!       if ( sst(i,j) < 273.16 ) sst(i,j) = -9999.
      end do
    end do
 
    call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
    write (stdout,*) 'XLON,XLAT,SST = ' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
!       ******           WRITE OUT SST DATA ON MM4 GRID
    call writerec(idate,.false.)
    write (stdout,*) 'WRITING OUT SST DATA:' , tochar(idate)
    close(11)

    idate = idate + itbc

  end do

 
  end subroutine sst_eh5om
!
end module mod_sst_eh5om
