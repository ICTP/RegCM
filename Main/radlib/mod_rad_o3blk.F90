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
!    MDRCHANTABILITY or FITNDSS FOR A PARTICULAR PURPOSD.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_rad_o3blk

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_date
  use mod_interp
  use mod_vertint
  use mod_mppparam
  use mod_mpmessage
  use mod_dynparam
  use mod_memutil
  use mod_rad_common
  use mod_stdio
  use netcdf

  private

  public :: allocate_mod_rad_o3blk , o3data , read_o3data

  real(rk8) , dimension(31) :: o3ann , o3sum , o3win , o3wrk , ppann ,&
                              ppsum , ppwin , ppwrk
  real(rk8) , dimension(32) :: ppwrkh
  real(rk8) , pointer , dimension(:) :: prlevh

  real(rk8) , pointer , dimension(:,:) :: alon , alat , aps , laps
  real(rk8) , pointer , dimension(:,:,:) :: ozone1 , ozone2
  real(rk8) , pointer , dimension(:) :: asig
  real(rk8) , pointer , dimension(:,:,:) :: ozone
!
  data o3sum/5.297D-8 , 5.852D-8 , 6.579D-8 , 7.505D-8 , 8.577D-8 , &
             9.895D-8 , 1.175D-7 , 1.399D-7 , 1.677D-7 , 2.003D-7 , &
             2.571D-7 , 3.325D-7 , 4.438D-7 , 6.255D-7 , 8.168D-7 , &
             1.036D-6 , 1.366D-6 , 1.855D-6 , 2.514D-6 , 3.240D-6 , &
             4.033D-6 , 4.854D-6 , 5.517D-6 , 6.089D-6 , 6.689D-6 , &
             1.106D-5 , 1.462D-5 , 1.321D-5 , 9.856D-6 , 5.960D-6 , &
             5.960D-6/
  data ppsum      /955.890D0 , 850.532D0 , 754.599D0 , 667.742D0 , &
       589.841D0 , 519.421D0 , 455.480D0 , 398.085D0 , 347.171D0 , &
       301.735D0 , 261.310D0 , 225.360D0 , 193.419D0 , 165.490D0 , &
       141.032D0 , 120.125D0 , 102.689D0 ,  87.829D0 ,  75.123D0 , &
        64.306D0 ,  55.086D0 ,  47.209D0 ,  40.535D0 ,  34.795D0 , &
        29.865D0 ,  19.122D0 ,   9.277D0 ,   4.660D0 ,   2.421D0 , &
         1.294D0 ,   0.647D0/
!
  data o3win/4.629D-8 , 4.686D-8 , 5.017D-8 , 5.613D-8 , 6.871D-8 , &
             8.751D-8 , 1.138D-7 , 1.516D-7 , 2.161D-7 , 3.264D-7 , &
             4.968D-7 , 7.338D-7 , 1.017D-6 , 1.308D-6 , 1.625D-6 , &
             2.011D-6 , 2.516D-6 , 3.130D-6 , 3.840D-6 , 4.703D-6 , &
             5.486D-6 , 6.289D-6 , 6.993D-6 , 7.494D-6 , 8.197D-6 , &
             9.632D-6 , 1.113D-5 , 1.146D-5 , 9.389D-6 , 6.135D-6 , &
             6.135D-6/
  data ppwin      /955.747D0 , 841.783D0 , 740.199D0 , 649.538D0 , &
       568.404D0 , 495.815D0 , 431.069D0 , 373.464D0 , 322.354D0 , &
       277.190D0 , 237.635D0 , 203.433D0 , 174.070D0 , 148.949D0 , &
       127.408D0 , 108.915D0 ,  93.114D0 ,  79.551D0 ,  67.940D0 , &
        58.072D0 ,  49.593D0 ,  42.318D0 ,  36.138D0 ,  30.907D0 , &
        26.362D0 ,  16.423D0 ,   7.583D0 ,   3.620D0 ,   1.807D0 , &
         0.938D0 ,   0.469D0/

  contains

  subroutine allocate_mod_rad_o3blk
    implicit none
    call getmem1d(prlevh,1,kzp2,'mod_o3blk:prlevh')

    if ( iclimao3 == 1 ) then
      call getmem2d(laps,jci1,jci2,ici1,ici2,'mod_o3blk:laps')
      if ( myid == iocpu ) then
        call getmem1d(asig,1,kzp1,'mod_o3blk:asig')
        call getmem2d(alon,jcross1,jcross2,icross1,icross2,'mod_o3blk:alon')
        call getmem2d(alat,jcross1,jcross2,icross1,icross2,'mod_o3blk:alat')
        call getmem2d(aps,jcross1,jcross2,icross1,icross2,'mod_o3blk:aps')
        call getmem3d(ozone1,jcross1,jcross2, &
                             icross1,icross2,1,kzp1,'mod_o3blk:ozone1')
        call getmem3d(ozone2,jcross1,jcross2, &
                             icross1,icross2,1,kzp1,'mod_o3blk:ozone2')
        call getmem3d(ozone,jcross1,jcross2, &
                            icross1,icross2,1,kzp1,'mod_o3blk:ozone')
      end if
    end if
  end subroutine allocate_mod_rad_o3blk
!
!----------------------------------------------------------------------
!
  subroutine o3data
    implicit none
!
    integer(ik4) :: i , j , jj , k , kj
    real(rk8) :: pb1 , pb2 , pt1 , pt2
!
    do k = 1 , 31
      ppann(k) = ppsum(k)
    end do
    o3ann(1) = 0.5D0*(o3sum(1)+o3win(1))
!
    do k = 2 , 31
      o3ann(k) = o3win(k-1) + (o3win(k)-o3win(k-1)) / &
                 (ppwin(k)-ppwin(k-1))*(ppsum(k)-ppwin(k-1))
    end do
    do k = 2 , 31
      o3ann(k) = 0.5D0*(o3ann(k)+o3sum(k))
    end do
    do k = 1 , 31
      o3wrk(k) = o3ann(k)
      ppwrk(k) = ppann(k)
    end do
    !
    ! calculate half pressure levels for model and data levels
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        do k = kzp1 , 1 , -1
          kj = kzp1 - k + 1
          prlevh(kj) = (flev(k)*sfps(j,i)+ptp)*d_10
        end do
        ppwrkh(1) = 1100.0D0
        do k = 2 , 31
          ppwrkh(k) = (ppwrk(k)+ppwrk(k-1))*d_half
        end do
        ppwrkh(32) = d_zero
        do k = 1 , kz
          o3prof(j,i,k) = d_zero
          do jj = 1 , 31
            if ( (-(prlevh(k)-ppwrkh(jj))) >= d_zero ) then
              pb1 = d_zero
            else
              pb1 = prlevh(k) - ppwrkh(jj)
            end if
            if ( (-(prlevh(k)-ppwrkh(jj+1))) >= d_zero ) then
              pb2 = d_zero
            else
              pb2 = prlevh(k) - ppwrkh(jj+1)
            end if
            if ( (-(prlevh(k+1)-ppwrkh(jj))) >= d_zero ) then
              pt1 = d_zero
            else
              pt1 = prlevh(k+1) - ppwrkh(jj)
            end if
            if ( (-(prlevh(k+1)-ppwrkh(jj+1))) >= d_zero ) then
              pt2 = d_zero
            else
              pt2 = prlevh(k+1) - ppwrkh(jj+1)
            end if
            o3prof(j,i,k) = o3prof(j,i,k) + (pb2-pb1-pt2+pt1)*o3wrk(jj)
          end do
          o3prof(j,i,k) = o3prof(j,i,k)/(prlevh(k)-prlevh(k+1))
        end do
      end do
    end do
!
  end subroutine o3data
!
  subroutine read_o3data(idatex,scenario,xlat,xlon,ps,ptop,sigma)
    implicit none
    type (rcm_time_and_date) , intent(in) :: idatex
    real(rk8) , pointer , dimension(:,:) :: xlat , xlon , ps
    real(rk8) , pointer , dimension(:) :: sigma
    real(rk8) , intent(in) :: ptop
    character(len=8) , intent(in) :: scenario
!
    character(len=64) :: infile
    logical , save :: ifirst
    logical :: dointerp
    real(rk8) , dimension(72,37,24) :: xozone1 , xozone2
    real(rk8) , dimension(njcross,nicross,24) :: yozone
    real(rk8) , save , dimension(37) :: lat
    real(rk8) , save , dimension(72) :: lon
    real(rk8) , save , dimension(24) :: plev
    real(rk8) :: xfac1 , xfac2 , odist
    type (rcm_time_and_date) :: imonmidd
    integer(ik4) :: iyear , imon , iday , ihour
    integer(ik4) , save :: ncid = -1
    integer(ik4) :: im1 , iy1 , im2 , iy2
    integer(ik4) , save :: ism , isy
    type (rcm_time_and_date) :: iref1 , iref2
    type (rcm_time_interval) :: tdif
    data ifirst /.true./
!
    call split_idate(idatex,iyear,imon,iday,ihour)
    imonmidd = monmiddle(idatex)

    if ( (iyear < 1850 .and. iyear > 2099) .or. &
         (iyear == 2100 .and. imon /= 1 ) ) then
      write (stderr,*) 'NO CLIMATIC O3 DATA AVAILABLE FOR ',iyear*100+imon
      write (stderr,*) 'WILL USE TABULATED VALUES.'
      return
    end if

    if ( myid == iocpu ) then
      if ( ifirst ) then
        alon = real(xlon(jcross1:jcross2,icross1:icross2))
        alat = real(xlat(jcross1:jcross2,icross1:icross2))
        asig = real(sigma(1:kzp1))
        ifirst = .false.
      end if

      if ( iyear < 2010 ) then
        infile = 'Ozone_CMIP5_ACC_SPARC_RF.nc'
      else if ( scenario(4:6) == '2.6' ) then
        infile = 'Ozone_CMIP5_ACC_SPARC_RCP2.6.nc'
      else if ( scenario(4:6) == '4.5' ) then
        infile = 'Ozone_CMIP5_ACC_SPARC_RCP4.5.nc'
      else if ( scenario(4:6) == '8.5' ) then
        infile = 'Ozone_CMIP5_ACC_SPARC_RCP8.5.nc'
      else
        call fatal(__FILE__,__LINE__,'ONLY ACCEPTED RCP SCENARIOS')
      end if
    end if

    im1 = imon
    iy1 = iyear
    im2 = imon
    iy2 = iyear
    if ( idatex > imonmidd ) then
      call inextmon(iy2,im2)
      ism = im1
      isy = iy1
      iref1 = imonmidd
      iref2 = monmiddle(nextmon(idatex))
    else
      call iprevmon(iy1,im1)
      ism = im1
      isy = iy1
      iref1 = monmiddle(prevmon(idatex))
      iref2 = imonmidd
    end if
    dointerp = .false.
    if ( ncid < 0 ) then
      if ( myid == iocpu ) then
        call init_o3data(infile,ncid,lat,lon,plev)
      else
        ncid = 0
      end if
      ism = im1
      isy = iy1
      dointerp = .true.
    else
      if ( ism /= im1 .or. isy /= iy1 ) then
        ism = im1
        isy = iy1
        dointerp = .true.
      end if
    end if

    if ( dointerp ) then
      ! We need pressure
      laps = real((ps(jci1:jci2,ici1:ici2)))
      call grid_collect(laps,aps,jci1,jci2,ici1,ici2)
      if ( myid == iocpu ) then
        write (stderr,*) 'Reading Ozone Data...'
        call readvar3d_pack(ncid,iy1,im1,'ozone',xozone1)
        call readvar3d_pack(ncid,iy2,im2,'ozone',xozone2)
        call bilinx2(yozone,xozone1,alon,alat,lon,lat,72,37,njcross,nicross,24)
        call intv0(ozone1,yozone,aps,asig,plev,ptop,njcross,nicross,kzp1,24)
        call bilinx2(yozone,xozone2,alon,alat,lon,lat,72,37,njcross,nicross,24)
        call intv0(ozone2,yozone,aps,asig,plev,ptop,njcross,nicross,kzp1,24)
      end if
    end if

    if ( myid == iocpu ) then
      tdif = idatex-iref1
      xfac1 = tohours(tdif)
      tdif = idatex-iref2
      xfac2 = tohours(tdif)
      odist = xfac1 - xfac2
      xfac1 = xfac1/odist
      xfac2 = d_one-xfac1
      ozone = (ozone1*xfac2+ozone2*xfac1)*1.0D-06
    end if
    call grid_distribute(ozone,o3prof,jci1,jci2,ici1,ici2,1,kzp1)
  end subroutine read_o3data
!
  subroutine inextmon(iyear,imon)
    implicit none
    integer(ik4) , intent(inout) :: iyear , imon
    imon = imon + 1
    if ( imon > 12 ) then
      imon = 1
      iyear = iyear + 1
    end if
  end subroutine inextmon
!
  subroutine iprevmon(iyear,imon)
    implicit none
    integer(ik4) , intent(inout) :: iyear , imon
    imon = imon - 1
    if ( imon < 1 ) then
      imon = 12
      iyear = iyear - 1
    end if
  end subroutine iprevmon
!
  subroutine init_o3data(o3file,ncid,lat,lon,plev)
    implicit none
    character(len=*) , intent(in) :: o3file
    integer(ik4) , intent(out) :: ncid
    real(rk8) , intent(out) , dimension(:) :: lat , lon , plev
    integer(ik4) :: iret
    iret = nf90_open(o3file,nf90_nowrite,ncid)
    if ( iret /= nf90_noerr ) then
      write (stderr, *) nf90_strerror(iret) , o3file
      call fatal(__FILE__,__LINE__,'CANNOT OPEN OZONE FILE')
    end if
    call readvar1d(ncid,'latitude',lat)
    call readvar1d(ncid,'longitude',lon)
    call readvar1d(ncid,'level',plev)
    plev(:) = plev(:) * 0.001
  end subroutine init_o3data

  subroutine readvar3d_pack(ncid,iyear,imon,vname,val)
    implicit none
    integer(ik4) , intent(in) :: ncid , iyear , imon
    character(len=*) , intent(in) :: vname
    real(rk8) , intent(out) , dimension(:,:,:) :: val
    real(rk8) , save :: xscale , xfact
    integer(ik4) , save :: ilastncid , icvar
    integer(ik4) , save , dimension(4) :: istart , icount
    integer(ik4) :: iret , irec
    data ilastncid /-1/
    data icvar /-1/
    data xscale /1.0D0/
    data xfact /0.0D0/
    data istart  /  1 ,  1 ,  1 ,  1/
    data icount  / 72 , 37 , 24 ,  1/

    irec = ((iyear-1849)*12+imon-12)+1
    if ( ncid /= ilastncid ) then
      iret = nf90_inq_varid(ncid,vname,icvar)
      if ( iret /= nf90_noerr ) then
        write (stderr, *) nf90_strerror(iret)
        call fatal(__FILE__,__LINE__,'CANNOT READ FROM OZONE FILE')
      end if
      iret = nf90_get_att(ncid,icvar,'scale_factor',xscale)
      if ( iret /= nf90_noerr ) then
        write (stderr, *) nf90_strerror(iret)
        call fatal(__FILE__,__LINE__,'CANNOT READ FROM OZONE FILE')
      end if
      iret = nf90_get_att(ncid,icvar,'add_offset',xfact)
      if ( iret /= nf90_noerr ) then
        write (stderr, *) nf90_strerror(iret)
        call fatal(__FILE__,__LINE__,'CANNOT READ FROM OZONE FILE')
      end if
    end if
    istart(4) = irec
    iret = nf90_get_var(ncid,icvar,val,istart,icount)
    if ( iret /= nf90_noerr ) then
      write (stderr, *) nf90_strerror(iret)
      call fatal(__FILE__,__LINE__,'CANNOT READ FROM OZONE FILE')
    end if
    val = dble(val*xscale+xfact)
  end subroutine readvar3d_pack

  subroutine readvar1d(ncid,vname,val)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , intent(out) , dimension(:) :: val
    integer(ik4) :: icvar , iret
    iret = nf90_inq_varid(ncid,vname,icvar)
    if ( iret /= nf90_noerr ) then
      write (stderr, *) nf90_strerror(iret)
      call fatal(__FILE__,__LINE__,'CANNOT READ FROM OZONE FILE')
    end if
    iret = nf90_get_var(ncid,icvar,val)
    if ( iret /= nf90_noerr ) then
      write (stderr, *) nf90_strerror(iret)
      call fatal(__FILE__,__LINE__,'CANNOT READ FROM OZONE FILE')
    end if
  end subroutine readvar1d
!
end module mod_rad_o3blk
