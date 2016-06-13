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
  use mod_runparams
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
  use mod_regcm_types
  use netcdf

  implicit none

  private

  public :: allocate_mod_rad_o3blk , o3data , read_o3data

  real(rkx) , dimension(31) :: o3ann , o3sum , o3win , ppann ,&
                              ppsum , ppwin

  real(rkx) , pointer , dimension(:) :: plev
  real(rkx) , pointer , dimension(:,:) :: alon , alat , aps
  real(rkx) , pointer , dimension(:,:,:) :: ozone1 , ozone2
  real(rkx) , pointer , dimension(:,:,:) :: ozone , pp3d
  real(rkx) , pointer , dimension(:,:,:) :: yozone

  data o3sum/5.297e-8_rkx , 5.852e-8_rkx , 6.579e-8_rkx , 7.505e-8_rkx , &
             8.577e-8_rkx , 9.895e-8_rkx , 1.175e-7_rkx , 1.399e-7_rkx , &
             1.677e-7_rkx , 2.003e-7_rkx , 2.571e-7_rkx , 3.325e-7_rkx , &
             4.438e-7_rkx , 6.255e-7_rkx , 8.168e-7_rkx , 1.036e-6_rkx , &
             1.366e-6_rkx , 1.855e-6_rkx , 2.514e-6_rkx , 3.240e-6_rkx , &
             4.033e-6_rkx , 4.854e-6_rkx , 5.517e-6_rkx , 6.089e-6_rkx , &
             6.689e-6_rkx , 1.106e-5_rkx , 1.462e-5_rkx , 1.321e-5_rkx , &
             9.856e-6_rkx , 5.960e-6_rkx , 5.960e-6_rkx/
  data ppsum       / 955.890_rkx , 850.532_rkx , 754.599_rkx , 667.742_rkx , &
       589.841_rkx , 519.421_rkx , 455.480_rkx , 398.085_rkx , 347.171_rkx , &
       301.735_rkx , 261.310_rkx , 225.360_rkx , 193.419_rkx , 165.490_rkx , &
       141.032_rkx , 120.125_rkx , 102.689_rkx ,  87.829_rkx ,  75.123_rkx , &
        64.306_rkx ,  55.086_rkx ,  47.209_rkx ,  40.535_rkx ,  34.795_rkx , &
        29.865_rkx ,  19.122_rkx ,   9.277_rkx ,   4.660_rkx ,   2.421_rkx , &
         1.294_rkx ,   0.647_rkx/

  data o3win/4.629e-8_rkx , 4.686e-8_rkx , 5.017e-8_rkx , 5.613e-8_rkx , &
             6.871e-8_rkx , 8.751e-8_rkx , 1.138e-7_rkx , 1.516e-7_rkx , &
             2.161e-7_rkx , 3.264e-7_rkx , 4.968e-7_rkx , 7.338e-7_rkx , &
             1.017e-6_rkx , 1.308e-6_rkx , 1.625e-6_rkx , 2.011e-6_rkx , &
             2.516e-6_rkx , 3.130e-6_rkx , 3.840e-6_rkx , 4.703e-6_rkx , &
             5.486e-6_rkx , 6.289e-6_rkx , 6.993e-6_rkx , 7.494e-6_rkx , &
             8.197e-6_rkx , 9.632e-6_rkx , 1.113e-5_rkx , 1.146e-5_rkx , &
             9.389e-6_rkx , 6.135e-6_rkx , 6.135e-6_rkx/
  data ppwin       / 955.747_rkx , 841.783_rkx , 740.199_rkx , 649.538_rkx , &
       568.404_rkx , 495.815_rkx , 431.069_rkx , 373.464_rkx , 322.354_rkx , &
       277.190_rkx , 237.635_rkx , 203.433_rkx , 174.070_rkx , 148.949_rkx , &
       127.408_rkx , 108.915_rkx ,  93.114_rkx ,  79.551_rkx ,  67.940_rkx , &
        58.072_rkx ,  49.593_rkx ,  42.318_rkx ,  36.138_rkx ,  30.907_rkx , &
        26.362_rkx ,  16.423_rkx ,   7.583_rkx ,   3.620_rkx ,   1.807_rkx , &
         0.938_rkx ,   0.469_rkx/

  contains

  subroutine allocate_mod_rad_o3blk
    implicit none
    if ( iclimao3 == 1 ) then
      if ( myid == iocpu ) then
        call getmem2d(alon,jcross1,jcross2,icross1,icross2,'mod_o3blk:alon')
        call getmem2d(alat,jcross1,jcross2,icross1,icross2,'mod_o3blk:alat')
        call getmem2d(aps,jcross1,jcross2,icross1,icross2,'mod_o3blk:aps')
        call getmem3d(ozone1,jcross1,jcross2, &
                             icross1,icross2,1,kzp1,'mod_o3blk:ozone1')
        call getmem3d(ozone2,jcross1,jcross2, &
                             icross1,icross2,1,kzp1,'mod_o3blk:ozone2')
        call getmem3d(ozone,jcross1,jcross2, &
                            icross1,icross2,1,kzp1,'mod_o3blk:ozone')
        call getmem3d(pp3d,jcross1,jcross2, &
                           icross1,icross2,1,kzp1,'mod_o3blk:pp3d')
      end if
    end if
  end subroutine allocate_mod_rad_o3blk

  subroutine o3data(m2r)
    implicit none
    type(mod_2_rad) , intent(in) :: m2r
    integer(ik4) :: k
    real(rkx) , dimension(kzp1) :: ozprnt
    real(rkx) , pointer , dimension(:) :: o3wrk , ppwrk
    allocate(o3wrk(31),ppwrk(31))
    do k = 1 , 31
      ppann(k) = ppsum(k)
    end do
    o3ann(1) = 0.5_rkx*(o3sum(1)+o3win(1))
    do k = 2 , 31
      o3ann(k) = o3win(k-1) + (o3win(k)-o3win(k-1)) / &
                 (ppwin(k)-ppwin(k-1))*(ppsum(k)-ppwin(k-1))
    end do
    do k = 2 , 31
      o3ann(k) = 0.5_rkx*(o3ann(k)+o3sum(k))
    end do
    do k = 1 , 31
      o3wrk(k) = o3ann(k)
      ppwrk(k) = ppann(k)
    end do
    ppwrk(:) = ppwrk(:) * d_100
    call intlinprof(o3prof,o3wrk,m2r%psatms,m2r%pfatms,jci1,jci2, &
                    ici1,ici2,kzp1,ppwrk,31)
    if ( myid == italk ) then
      ozprnt = o3prof(3,3,:)
      call vprntv(ozprnt,kzp1,'Ozone profile at (3,3)')
    end if
    deallocate(o3wrk,ppwrk)
  end subroutine o3data

  subroutine read_o3data(idatex,scenario,m2r)
    implicit none
    type (rcm_time_and_date) , intent(in) :: idatex
    type(mod_2_rad) , intent(in) :: m2r
    character(len=8) , intent(in) :: scenario
    character(len=64) :: infile
    logical , save :: ifirst
    logical :: dointerp
    real(rkx) , dimension(kzp1) :: ozprnt
    real(rkx) , dimension(72,37,24) :: xozone1 , xozone2
    real(rkx) , save , dimension(37) :: lat
    real(rkx) , save , dimension(72) :: lon
    real(rkx) :: xfac1 , xfac2 , odist
    type (rcm_time_and_date) :: imonmidd
    integer(ik4) :: iyear , imon , iday , ihour
    integer(ik4) , save :: ncid = -1
    integer(ik4) :: im1 , iy1 , im2 , iy2
    integer(ik4) , save :: ism , isy
    type (rcm_time_and_date) :: iref1 , iref2
    type (rcm_time_interval) :: tdif
    data ifirst /.true./

    call split_idate(idatex,iyear,imon,iday,ihour)
    imonmidd = monmiddle(idatex)

    if ( (iyear < 1850 .and. iyear > 2099) .or. &
         (iyear == 2100 .and. imon /= 1 ) ) then
      write (stderr,*) 'NO CLIMATIC O3 DATA AVAILABLE FOR ',iyear*100+imon
      write (stderr,*) 'WILL USE TABULATED VALUES.'
      return
    end if

    if ( ifirst ) then
      call grid_collect(m2r%xlon,alon,jci1,jci2,ici1,ici2)
      call grid_collect(m2r%xlat,alat,jci1,jci2,ici1,ici2)
      ifirst = .false.
    end if

    if ( myid == iocpu ) then
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
        call getmem1d(plev,1,24,'ozone:plev')
        call getmem3d(yozone,1,njcross,1,nicross,1,24,'ozone:yozone')
        call init_o3data(infile,ncid,lat,lon)
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
      call grid_collect(m2r%psatms,aps,jci1,jci2,ici1,ici2)
      call grid_collect(m2r%pfatms,pp3d,jci1,jci2,ici1,ici2,1,kzp1)
      if ( myid == iocpu ) then
        write (stdout,*) 'Reading Ozone Data...'
        call readvar3d_pack(ncid,iy1,im1,'ozone',xozone1)
        call readvar3d_pack(ncid,iy2,im2,'ozone',xozone2)
        call bilinx2(yozone,xozone1,alon,alat,lon,lat,72,37,njcross,nicross,24)
        call intlinreg(ozone1,yozone,aps,pp3d,1,njcross,1,nicross,kzp1,plev,24)
        call bilinx2(yozone,xozone2,alon,alat,lon,lat,72,37,njcross,nicross,24)
        call intlinreg(ozone2,yozone,aps,pp3d,1,njcross,1,nicross,kzp1,plev,24)
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
      ozone = (ozone1*xfac2+ozone2*xfac1)*1.0e-6_rkx
    end if
    call grid_distribute(ozone,o3prof,jci1,jci2,ici1,ici2,1,kzp1)
    if ( myid == italk .and. dointerp ) then
      ozprnt = o3prof(3,3,:)
      call vprntv(ozprnt,kzp1,'Updated ozone profile at (3,3)')
    end if
  end subroutine read_o3data

  subroutine inextmon(iyear,imon)
    implicit none
    integer(ik4) , intent(inout) :: iyear , imon
    imon = imon + 1
    if ( imon > 12 ) then
      imon = 1
      iyear = iyear + 1
    end if
  end subroutine inextmon

  subroutine iprevmon(iyear,imon)
    implicit none
    integer(ik4) , intent(inout) :: iyear , imon
    imon = imon - 1
    if ( imon < 1 ) then
      imon = 12
      iyear = iyear - 1
    end if
  end subroutine iprevmon

  subroutine init_o3data(o3file,ncid,lat,lon)
    implicit none
    character(len=*) , intent(in) :: o3file
    integer(ik4) , intent(out) :: ncid
    real(rkx) , intent(out) , dimension(:) :: lat , lon
    integer(ik4) :: iret
    iret = nf90_open(o3file,nf90_nowrite,ncid)
    if ( iret /= nf90_noerr ) then
      write (stderr, *) nf90_strerror(iret) , o3file
      call fatal(__FILE__,__LINE__,'CANNOT OPEN OZONE FILE')
    end if
    call readvar1d(ncid,'latitude',lat)
    call readvar1d(ncid,'longitude',lon)
    call readvar1d(ncid,'level',plev)
    plev(:) = plev(:) * d_100 ! Put it in Pa
  end subroutine init_o3data

  subroutine readvar3d_pack(ncid,iyear,imon,vname,val)
    implicit none
    integer(ik4) , intent(in) :: ncid , iyear , imon
    character(len=*) , intent(in) :: vname
    real(rkx) , intent(out) , dimension(:,:,:) :: val
    real(rkx) , save :: xscale , xfact
    integer(ik4) , save :: ilastncid , icvar
    integer(ik4) , save , dimension(4) :: istart , icount
    integer(ik4) :: iret , irec
    data ilastncid /-1/
    data icvar /-1/
    data xscale /1.0_rkx/
    data xfact /0.0_rkx/
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
    val = real(val*xscale+xfact,rkx)
  end subroutine readvar3d_pack

  subroutine readvar1d(ncid,vname,val)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rkx) , intent(out) , dimension(:) :: val
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

end module mod_rad_o3blk

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
