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

module mod_sst_fvgcm

  use m_realkinds
  use m_die
  use m_stdio
  use m_zeit
  use mod_dynparam
  use mod_sst_grid
  use mod_interp

  private

  public :: sst_fvgcm

  contains

  subroutine sst_fvgcm

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Comments on dataset sources and location:                          c
!                                                                    c
! FVGCM    HadAMH_SST in the original netCDF format.                 c
!          for 'RF'          run, from 1959 to 1991, 385 months      c
!          for 'A2' and 'B2' run, from 2069 to 2101, 385 months      c
!                                                                    c
!      ML = 1 is   0.0; ML = 2 is   1.875; => ML = 192 is 358.125E   c
!      NL = 1 is  90.0; ML = 2 is  88.75 ; => ML = 145 is -90.       c
!                                                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  implicit none
!
  integer , parameter :: ilon = 192 , jlat = 145
!
  integer :: i , it , j , k , ludom , lumax , nsteps , iv
  real(sp) , dimension(jlat) :: lati
  real(sp) , dimension(ilon) :: loni
  integer , dimension(20) :: lund
  real(sp) , dimension(ilon,jlat) :: temp
  real(sp) , dimension(ilon,jlat) :: sst
  type(rcm_time_and_date) :: idate , idateo , idatef
  integer :: year , month , day , hour
  logical :: there
!
  call zeit_ci('sst_fvgcm')

  call split_idate(idate,year,month,day,hour)

  if ( ssttyp == 'FV_RF' ) then
    inquire (file=trim(inpglob)//'/SST/Sst_1959_1991ref.dat',exist=there)
    if ( .not.there ) then
      write (stderr,*)                                       &
          'Sst_1959_1991ref.dat is not available under ',trim(inpglob),'/SST/'
      call die('sst_fvgcm')
    end if
    open (11,file=trim(inpglob)//'/SST/Sst_1959_1991ref.dat',       &
          form='unformatted',recl=ilon*jlat*ibyte,access='direct')
  else if ( ssttyp == 'FV_A2' ) then
    inquire (file=trim(inpglob)//'/SST/Sst_2069_2101_A2.dat',       &
             exist=there)
    if ( .not.there ) then
      write (stderr,*)                                              &
          'Sst_2069_2101_A2.dat is not available under ',trim(inpglob),'/SST/'
      call die('sst_fvgcm')
    end if
    open (11,file=trim(inpglob)//'/SST/Sst_2069_2101_A2.dat',       &
          form='unformatted',recl=ilon*jlat*ibyte,access='direct')
  else if ( ssttyp == 'FV_B2' ) then
    inquire (file=trim(inpglob)//'/SST/Sst_2069_2101_B2.dat',       &
             exist=there)
    if ( .not.there ) then
      write (stderr,*)                                              &
           'Sst_2069_2101_B2.dat is not available under ',trim(inpglob),'/SST/'
      call die('sst_fvgcm')
    end if
    open (11,file=trim(inpglob)//'/SST/Sst_2069_2101_B2.dat',       &
          form='unformatted',recl=ilon*jlat*ibyte,access='direct')
  else
    write (stderr,*) 'PLEASE SET right SSTTYP in regcm.in'
    write (stderr,*) 'Supported types are FV_RF FV_A2 FV_B2'
    call die('sst fvgcm')
  end if

  idateo = monfirst(globidate1)
  if (lfhomonth(globidate1)) then
    idateo = prevmon(globidate1)
  end if
  idatef = monfirst(globidate2)
  if (idatef < globidate2) then
    idatef = nextmon(idatef)
  end if
  nsteps = imondiff(idatef,idateo) + 1

  call open_sstfile(idateo)
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
  do i = 1 , ilon
    loni(i) = float(i-1)*1.875
  end do
  do j = 1 , jlat
    lati(j) = -90. + 1.25*float(j-1)
  end do
 
  idate = idateo
  do k = 1 , nsteps

    if ( ssttyp == 'FV_RF' ) then
      it = (year-1959)*12 + month
    else
      it = (year-2069)*12 + month
    end if

    read (11,rec=it) temp

    do j = 1 , jlat
      do i = 1 , ilon
        if ( temp(i,j) > -9000.0 .and. temp(i,j) < 10000.0 ) then
          sst(i,j) = temp(i,j)
        else
          sst(i,j) = -9999.
        end if
      end do
    end do
 
    call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
    write (stdout,*) 'XLON,XLAT,SST = ' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
    do j = 1 , jx
      do i = 1 , iy
        if ( sstmm(i,j) < -5000. .and. &
             (lu(i,j) > 13.5 .and. lu(i,j) < 15.5) ) then
          do iv = 1 , 20
            lund(iv) = 0
          end do
          lund(nint(lu(i-1,j-1))) = lund(nint(lu(i-1,j-1))) + 2
          lund(nint(lu(i-1,j))) = lund(nint(lu(i-1,j))) + 3
          lund(nint(lu(i-1,j+1))) = lund(nint(lu(i-1,j+1))) + 2
          lund(nint(lu(i,j-1))) = lund(nint(lu(i,j-1))) + 3
          lund(nint(lu(i,j+1))) = lund(nint(lu(i,j+1))) + 3
          lund(nint(lu(i+1,j-1))) = lund(nint(lu(i+1,j-1))) + 2
          lund(nint(lu(i+1,j))) = lund(nint(lu(i+1,j))) + 3
          lund(nint(lu(i+1,j+1))) = lund(nint(lu(i+1,j+1))) + 2
          ludom = 18
          lumax = 0
          do iv = 1 , 20
            if ( iv <= 13 .or. iv >= 16 ) then
              if ( lund(iv) > lumax ) then
                ludom = iv
                lumax = lund(k)
              end if
            end if
          end do
          lu(i,j) = float(ludom)
          write (stdout,*) ludom , sstmm(i,j)
        end if
        if ( sstmm(i,j) > -100. ) then
          sstmm(i,j) = sstmm(i,j)
        else
          sstmm(i,j) = -9999.
        end if
      end do
    end do
 
    call writerec(idate,.false.)

    write (stdout,*) 'WRITTEN OUT SST DATA : ' , tochar(idate)
    idate = nextmon(idate)

  end do
 
  call zeit_co('sst_fvgcm')

  end subroutine sst_fvgcm
!
end module mod_sst_fvgcm
