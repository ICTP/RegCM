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

module mod_ox_icbc

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_grid
  use mod_wrtoxd
  use mod_interp
  use mod_date
  use mod_nchelper
  use netcdf

  private
!
  integer(ik4) :: nyear , month , nday , nhour
  integer(ik4) :: k , l
  integer(ik4) :: k0

  integer(ik4) , parameter :: oxilon = 128 , oxjlat = 64 , oxilev = 26 , oxitime = 12
  real(rk4) , dimension(oxilon) :: oxt42lon
  real(rk4) , dimension(oxjlat) :: oxt42lat
  real(rk4) , dimension(oxilev) :: oxt42hyam , oxt42hybm
!
! Oxidant climatology variables
!
  real(rk4) :: p0 , r4pt
  real(rk4) , dimension(oxilon,oxjlat) :: xps
  real(rk4) , dimension(oxilon,oxilev,oxjlat,oxitime,noxsp) :: oxv2
  real(rk4) , dimension(oxilon,oxjlat,oxitime) :: xps2
  real(rk4) , pointer , dimension(:,:) :: poxid_3
  real(rk4) , pointer , dimension(:,:,:,:) :: oxv3

  real(rk4) :: prcm , pmpi , pmpj
  integer(ik4) :: ncid , istatus

  public :: header_ox_icbc , get_ox_icbc , close_ox_icbc

  data ncid /-1/

  contains

  subroutine header_ox_icbc
    implicit none
    integer(ik4) :: ivarid , istatus , is

    r4pt = real(ptop)

    call getmem2d(poxid_3,1,jx,1,iy,'mod_ox_icbc:poxid_3')
    call getmem4d(oxv3,1,jx,1,iy,1,oxilev,1,noxsp,'mod_ox_icbc:oxv3')

    istatus = nf90_open(trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
                      'oxid_3d_64x128_L26_c030722.nc', nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error open oxid file')

    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lon')
    istatus = nf90_get_var(ncid,ivarid,oxt42lon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lon')
    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lat')
    istatus = nf90_get_var(ncid,ivarid,oxt42lat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lat')
    istatus = nf90_inq_varid(ncid,'hyam',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var hyam')
    istatus = nf90_get_var(ncid,ivarid,oxt42hyam)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var hyam')
    istatus = nf90_inq_varid(ncid,'hybm',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var hybm')
    istatus = nf90_get_var(ncid,ivarid,oxt42hybm)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var hybm')
    istatus = nf90_inq_varid(ncid,'P0',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var P0')
    istatus = nf90_get_var(ncid,ivarid,p0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var P0')
    p0 = p0*0.01
    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var PS')
    istatus = nf90_get_var(ncid,ivarid,xps2)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var PS')
    xps2 = xps2 * 0.01
    do is = 1 , noxsp
      istatus = nf90_inq_varid(ncid,oxspec(is),ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//oxspec(is))
      istatus = nf90_get_var(ncid,ivarid,oxv2(:,:,:,:,is))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//oxspec(is))
    end do
  end subroutine header_ox_icbc

  subroutine get_ox_icbc(idate)
    implicit none
!
    integer(ik4) :: i , is , j , k , k0
    type(rcm_time_and_date) , intent(in) :: idate
    real(rk4) , dimension(oxilon,oxjlat,oxilev) :: xinp
    real(rk4) :: wt1 , wt2
    type(rcm_time_and_date) :: d1 , d2
    type(rcm_time_interval) :: t1 , tt
    integer(ik4) :: m1 , m2

    d1 = monfirst(idate)
    d2 = nextmon(d1)
    m1 = getmonth(d1)
    m2 = getmonth(d2)
    t1 = idate-d1
    tt = d2-d1
    wt1 = real(tohours(t1)/tohours(tt))
    wt2 = 1.0 - wt1

    do is = 1 , noxsp
      do i = 1 , oxjlat
        do j = 1 , oxilon
          do l = 1 , oxilev
            xinp(j,i,l) = oxv2(j,l,i,m1,is)*wt2+oxv2(j,l,i,m2,is)*wt1
          end do
        end do
      end do
      call bilinx2(oxv3(:,:,:,is),xinp,xlon,xlat,oxt42lon,oxt42lat, &
                   oxilon,oxjlat,iy,jx,oxilev) 
    end do

    do i = 1 , oxjlat
      do j = 1 , oxilon
        xps(j,i) = xps2(j,i,m1)*wt1+xps2(j,i,m2)*wt2
      end do
    end do

    call bilinx2(poxid_3,xps,xlon,xlat,oxt42lon,oxt42lat, &
                 oxilon,oxjlat,iy,jx,1)
    do i = 1 , iy 
      do j = 1 , jx
        do l = 1 , kz
          prcm=((poxid_3(j,i)*0.1-r4pt)*sigma2(l)+r4pt)*10.
          k0 = -1
          do k = oxilev , 1 , -1
            pmpi = oxt42hyam(k)*p0+poxid_3(j,i)*oxt42hybm(k)
            k0 = k
            if (prcm > pmpi) exit
          end do
          if (k0 == oxilev) then
            pmpj = oxt42hyam(oxilev-1)*p0+poxid_3(j,i)*oxt42hybm(oxilev-1)
            pmpi = oxt42hyam(oxilev  )*p0+poxid_3(j,i)*oxt42hybm(oxilev  )
            do is = 1 , noxsp
              oxv4(j,i,l,is) = oxv3(j,i,oxilev,is) + &
                 (oxv3(j,i,oxilev-1,is) - oxv3(j,i,oxilev,is)) * &
                 (prcm-pmpi)/(pmpi-pmpj)
            end do
          else if (k0 >= 1) then
            pmpj = oxt42hyam(k0  )*p0+poxid_3(j,i)*oxt42hybm(k0  )
            pmpi = oxt42hyam(k0+1)*p0+poxid_3(j,i)*oxt42hybm(k0+1)
            wt1 = (prcm-pmpj)/(pmpi-pmpj)
            wt2 = 1.0 - wt1
            do is = 1 , noxsp
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
    if ( ncid > 0 ) then
      istatus = nf90_close(ncid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close oxid file')
    end if
  end subroutine close_ox_icbc

end module mod_ox_icbc
