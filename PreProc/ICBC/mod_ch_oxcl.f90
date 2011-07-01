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

module mod_ch_oxcl

  use mod_dynparam
  use mod_memutil
  use mod_grid
  use mod_wrtoxd
  use mod_interp
  use mod_date
  use m_die
  use m_realkinds
  use netcdf

  private
!
  integer :: nyear , month , nday , nhour
  integer :: k , l
  integer :: k0

  integer , parameter :: oxilon = 128 , oxjlat = 64 , oxilev = 26 , oxitime = 12
  real(sp) , dimension(oxilon) :: oxt42lon
  real(sp) , dimension(oxjlat) :: oxt42lat
  real(sp) , dimension(oxilev) :: oxt42hyam , oxt42hybm
!
! Oxidant climatology variables
!
  real(sp) :: p0
  real(sp) , dimension(oxilon,oxjlat) :: xps
  real(sp) , dimension(oxilon,oxjlat) :: poxid_2
  real(sp) , dimension(oxilon,oxilev,oxjlat,oxitime,noxsp) :: oxv2
  real(sp) , dimension(oxilon,oxjlat,oxitime) :: xps2
  real(sp) , pointer , dimension(:,:) :: poxid_3
  real(sp) , pointer , dimension(:,:,:,:) :: oxv3

  real(sp) :: prcm , pmpi , pmpj
  integer :: ncid , istatus

  public :: header_ch_oxcl , get_ch_oxcl , close_ch_oxcl

  contains

  subroutine header_ch_oxcl
    implicit none
    integer :: ivarid , istatus , is

    call getmem2d(poxid_3,1,jx,1,iy,'mod_ch_oxcl:poxid_3')
    call getmem4d(oxv3,1,jx,1,iy,1,oxilev,1,noxsp,'mod_ch_oxcl:oxv3')

    istatus = nf90_open(trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
                      'oxid_3d_64x128_L26_c030722.nc', nf90_nowrite, ncid)
    if ( istatus /= nf90_noerr ) then
      write (stderr,*) 'Cannot open input file'
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if

    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,oxt42lon)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,oxt42lat)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'hyam',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,oxt42hyam)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'hybm',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,oxt42hybm)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'P0',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,p0)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,xps2)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
    do is = 1 , noxsp
      istatus = nf90_inq_varid(ncid,oxspec(is),ivarid)
      if ( istatus /= nf90_noerr ) then
        call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
      end if
      istatus = nf90_get_var(ncid,ivarid,oxv2(:,:,:,:,is))
      if ( istatus /= nf90_noerr ) then
        call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
      end if
    end do
  end subroutine header_ch_oxcl

  subroutine get_ch_oxcl(idate)
    implicit none
!
    integer :: i , is , j , k , k0
    type(rcm_time_and_date) , intent(in) :: idate
    real(sp) , dimension(oxilon,oxjlat,oxilev) :: xinp
    real(sp) :: wt1 , wt2 , r4pt
    type(rcm_time_and_date) :: d1 , d2
    type(rcm_time_interval) :: t1 , tt
    integer :: m1 , m2

    d1 = monfirst(idate)
    d2 = nextmon(d1)
    m1 = d1%month
    m2 = d2%month
    t1 = idate-d1
    tt = d2-d1
    wt1 = real(t1%hours()/tt%hours())
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

    poxid_2 = xps*0.01
    p0 = p0*0.01
    r4pt = real(ptop)

    call bilinx2(poxid_3,poxid_2,xlon,xlat,oxt42lon,oxt42lat, &
                 oxilon,oxjlat,iy,jx,1)

    do i = 1 , iy 
      do j = 1 , jx
        do l = 1 , kz
          prcm=((poxid_3(j,i)*0.1-r4pt)*sigma2(l)+r4pt)*10.
          k0 = -1
          do k = oxilev , 1 , -1
            pmpi = poxid_3(j,i)*oxt42hybm(k)+oxt42hyam(k)*p0
            k0 = k
            if (prcm > pmpi) exit
          end do
          if (k0 == oxilev) then        
            pmpj = poxid_3(j,i)*oxt42hybm(oxilev-1)+oxt42hyam(oxilev-1)*p0
            pmpi = poxid_3(j,i)*oxt42hybm(oxilev)+oxt42hyam(oxilev)*p0

            do is = 1 , noxsp
              oxv4(j,i,l,is) = oxv3(j,i,oxilev,is) + &
                 (oxv3(j,i,oxilev,is) - oxv3(j,i,oxilev-1,is))*(prcm-pmpi)/(pmpi-pmpj)
            end do
          else if (k0 >= 1) then
            pmpj = poxid_3(j,i)*oxt42hybm(k0)+oxt42hyam(k0)*p0
            pmpi = poxid_3(j,i)*oxt42hybm(k0+1)+oxt42hyam(k0+1)*p0
            do is = 1 , noxsp
              oxv4(j,i,l,is) = (oxv3(j,i,k0+1,is)*(prcm-pmpj) + &
                                oxv3(j,i,k0,is)*(prcm-pmpi))/(pmpi-pmpj)
            end do
          end if            
        end do
      end do
    end do            

    call write_ch_oxcl(idate)

  end subroutine get_ch_oxcl

  subroutine close_ch_oxcl
    use netcdf
    implicit none
    istatus=nf90_close(ncid)
    if ( istatus/=nf90_noerr ) then
      write (stderr,*) 'Cannot close input file'
      call die('header_ch_oxcl',nf90_strerror(istatus),istatus)
    end if
  end subroutine close_ch_oxcl

end module mod_ch_oxcl
