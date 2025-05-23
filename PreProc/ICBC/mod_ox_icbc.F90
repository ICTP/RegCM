!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_ox_icbc

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_memutil
  use mod_stdio
  use mod_grid
  use mod_wrtoxd
  use mod_kdinterp
  use mod_date
  use mod_nchelper
  use netcdf

  private

  integer(ik4) :: l

  integer(ik4) :: oxilon, oxjlat, oxilev, oxitime
  real(rkx), pointer, contiguous, dimension(:) :: oxt42lon
  real(rkx), pointer, contiguous, dimension(:) :: oxt42lat
  real(rkx), pointer, contiguous, dimension(:) :: oxt42hyam, oxt42hybm
!
! Oxidant climatology variables
!
  real(rkx) :: p0, r4pt
  real(rkx), pointer, contiguous, dimension(:,:) :: xps
  real(rkx), pointer, contiguous, dimension(:,:,:,:,:) :: oxv2
  real(rkx), pointer, contiguous, dimension(:,:,:) :: xps2, xinp
  real(rkx), pointer, contiguous, dimension(:,:) :: poxid_3, xps3
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: oxv3

  real(rkx) :: prcm, pmpi, pmpj
  integer(ik4) :: ncid, istatus
  integer(ik4) :: ncicbc, ivarps, irec

  type(h_interpolator) :: hint
  type(rcm_time_and_date) :: iodate

  public :: init_ox_icbc, get_ox_icbc, close_ox_icbc

  data ncid /-1/
  data ncicbc /-1/

  contains

  subroutine init_ox_icbc(idate)
    implicit none
    type(rcm_time_and_date) :: idate
    integer(ik4) :: ivarid, istatus, dimid, is
    character(len=256) :: oxifile, icbcfilename

    r4pt = real(ptop)

    iodate = idate

    oxifile = trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
              'oxid_3d_64x128_L26_c030722.nc'
    write (icbcfilename,'(a,a,a,a,a,a)') trim(dirglob), pthsep, &
           trim(domname), '_ICBC.', trim(tochar10(idate)), '.nc'

    write(stdout,*) 'Opening ',trim(oxifile)
    istatus = nf90_open(oxifile, nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open oxid file')

    istatus = nf90_open(icbcfilename,nf90_nowrite, ncicbc)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open ICBC file '//trim(icbcfilename))
    istatus = nf90_inq_varid(ncicbc,'ps',ivarps)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var ps in icbc file '//trim(icbcfilename))
    irec = 1

    istatus = nf90_inq_dimid(ncid,'lon',dimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim lon')
    istatus = nf90_inquire_dimension(ncid,dimid,len=oxilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lon')
    istatus = nf90_inq_dimid(ncid,'lat',dimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim lat')
    istatus = nf90_inquire_dimension(ncid,dimid,len=oxjlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lat')
    istatus = nf90_inq_dimid(ncid,'lev',dimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim lev')
    istatus = nf90_inquire_dimension(ncid,dimid,len=oxilev)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lev')
    istatus = nf90_inq_dimid(ncid,'time',dimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim time')
    istatus = nf90_inquire_dimension(ncid,dimid,len=oxitime)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire time')

    call getmem1d(oxt42lon,1,oxilon,'mod_ox_icbc:oxilon')
    call getmem1d(oxt42lat,1,oxjlat,'mod_ox_icbc:oxjlat')
    call getmem1d(oxt42hyam,1,oxilev,'mod_ox_icbc:oxt42hyam')
    call getmem1d(oxt42hybm,1,oxilev,'mod_ox_icbc:oxt42hybm')

    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var lon')
    istatus = nf90_get_var(ncid,ivarid,oxt42lon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var lon')
    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var lat')
    istatus = nf90_get_var(ncid,ivarid,oxt42lat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var lat')
    istatus = nf90_inq_varid(ncid,'hyam',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var hyam')
    istatus = nf90_get_var(ncid,ivarid,oxt42hyam)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var hyam')
    istatus = nf90_inq_varid(ncid,'hybm',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var hybm')
    istatus = nf90_get_var(ncid,ivarid,oxt42hybm)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var hybm')
    istatus = nf90_inq_varid(ncid,'P0',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var P0')
    istatus = nf90_get_var(ncid,ivarid,p0)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var P0')

    call h_interpolator_create(hint,oxt42lat,oxt42lon,xlat,xlon)

    call getmem2d(poxid_3,1,jx,1,iy,'mod_ox_icbc:poxid_3')
    call getmem2d(xps3,1,jx,1,iy,'mod_ch_icbc:xps3')
    call getmem2d(xps,1,oxilon,1,oxjlat,'mod_ox_icbc:xps')
    call getmem3d(xps2,1,oxilon,1,oxjlat,1,oxitime,'mod_ox_icbc:xps2')
    call getmem3d(xinp,1,oxilon,1,oxjlat,1,oxilev,'mod_ox_icbc:xinp')
    call getmem4d(oxv3,1,jx,1,iy,1,oxilev,1,noxsp,'mod_ox_icbc:oxv3')
    call getmem5d(oxv2,1,oxilon,1,oxilev,1,oxjlat,1,oxitime, &
                  1,noxsp,'mod_ox_icbc:oxv2')

    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var PS')
    istatus = nf90_get_var(ncid,ivarid,xps2)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var PS')
    do is = 1, noxsp
      istatus = nf90_inq_varid(ncid,oxspec(is),ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//oxspec(is))
      istatus = nf90_get_var(ncid,ivarid,oxv2(:,:,:,:,is))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//oxspec(is))
    end do
  end subroutine init_ox_icbc

  subroutine get_ox_icbc(idate)
    implicit none

    integer(ik4) :: i, is, j, k, k0
    type(rcm_time_and_date), intent(in) :: idate
    real(rkx) :: wt1, wt2
    character(len=256) :: icbcfilename
    integer(ik4), dimension(3) :: istart, icount
    type(rcm_time_and_date) :: d1, d2
    type(rcm_time_interval) :: t1, tt
    integer(ik4) :: m1, m2, istatus

    d1 = monfirst(idate)
    d2 = nextmon(d1)
    m1 = getmonth(d1)
    m2 = getmonth(d2)
    t1 = idate-d1
    tt = d2-d1
    wt1 = real(tohours(t1)/tohours(tt),rkx)
    wt2 = 1.0_rkx - wt1

    if (.not. lsamemonth(idate, iodate) ) then
      istatus = nf90_close(ncicbc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close ICBC file')
      write (icbcfilename,'(a,a,a,a,a,a)') trim(dirglob), pthsep, &
             trim(domname), '_ICBC.', trim(tochar10(idate)), '.nc'
      istatus = nf90_open(icbcfilename,nf90_nowrite, ncicbc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open ICBC file '//trim(icbcfilename))
      istatus = nf90_inq_varid(ncicbc,'ps',ivarps)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var ps in icbc file '//trim(icbcfilename))
      iodate = idate
      irec = 1
    end if

    istart(1) = 1
    istart(2) = 1
    istart(3) = irec
    icount(1) = jx
    icount(2) = iy
    icount(3) = 1
    istatus = nf90_get_var(ncicbc,ivarps,poxid_3,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var ps')
    irec = irec + 1

    do is = 1, noxsp
      do i = 1, oxjlat
        do j = 1, oxilon
          do l = 1, oxilev
            xinp(j,i,l) = oxv2(j,l,i,m1,is)*wt2+oxv2(j,l,i,m2,is)*wt1
          end do
        end do
      end do
      call h_interpolate_cont(hint,xinp,oxv3(:,:,:,is))
    end do

    do i = 1, oxjlat
      do j = 1, oxilon
        xps(j,i) = xps2(j,i,m1)*wt1+xps2(j,i,m2)*wt2
      end do
    end do

    call h_interpolate_cont(hint,xps,xps3)

    do i = 1, iy
      do j = 1, jx
        do l = 1, kz
          prcm = ((poxid_3(j,i)*0.1_rkx-r4pt)*sigmah(l)+r4pt)*1000.0_rkx
          k0 = -1
          do k = oxilev, 1, -1
            pmpi = oxt42hyam(k)*p0+xps3(j,i)*oxt42hybm(k)
            k0 = k
            if (prcm > pmpi) exit
          end do
          if (k0 == oxilev) then
            pmpj = oxt42hyam(oxilev-1)*p0+xps3(j,i)*oxt42hybm(oxilev-1)
            pmpi = oxt42hyam(oxilev  )*p0+xps3(j,i)*oxt42hybm(oxilev  )
            do is = 1, noxsp
              oxv4(j,i,l,is) = max(oxv3(j,i,oxilev,is) + &
                 (oxv3(j,i,oxilev-1,is) - oxv3(j,i,oxilev,is)) * &
                 (prcm-pmpi)/(pmpi-pmpj),d_zero)
            end do
          else if (k0 >= 1) then
            pmpj = oxt42hyam(k0  )*p0+xps3(j,i)*oxt42hybm(k0  )
            pmpi = oxt42hyam(k0+1)*p0+xps3(j,i)*oxt42hybm(k0+1)
            wt1 = (prcm-pmpj)/(pmpi-pmpj)
            wt2 = 1.0 - wt1
            do is = 1, noxsp
              oxv4(j,i,l,is) = oxv3(j,i,k0+1,is)*wt1 + oxv3(j,i,k0,is)*wt2
            end do
          end if
        end do
      end do
    end do

    call write_ox_icbc(idate)

  end subroutine get_ox_icbc

  subroutine close_ox_icbc
    use netcdf
    implicit none
    call h_interpolator_destroy(hint)
    if ( ncid > 0 ) then
      istatus = nf90_close(ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close oxid file')
    end if
    if ( ncicbc > 0 ) then
      istatus = nf90_close(ncicbc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close icbc file')
    end if
  end subroutine close_ox_icbc

end module mod_ox_icbc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
