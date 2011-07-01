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

module mod_ch_icbc

  use mod_dynparam
  use mod_grid
  use mod_wrtoxd
  use mod_interp
  use mod_date
  use m_die
  use m_realkinds
  use mod_memutil
  use netcdf

  private
!
  integer :: chilon , chjlat , chilev

  real(sp) , pointer , dimension(:) :: cht42lon
  real(sp) , pointer , dimension(:) :: cht42lat
  real(sp) , pointer , dimension(:) :: cht42hyam , cht42hybm
!
! Oxidant climatology variables
!
  real(sp) :: p0
  real(sp) , pointer , dimension(:,:) :: poxid_3
  real(sp) , pointer , dimension(:,:,:,:) :: chv3
  real(sp) , pointer , dimension(:,:) :: xps
  real(sp) , pointer , dimension(:,:) :: poxid_2
  real(sp) , pointer , dimension(:,:,:,:) :: xinp

  real(sp) :: prcm , pmpi , pmpj
  integer :: ncid , istatus

  public :: header_ch_icbc , get_ch_icbc , close_ch_icbc

  contains

  subroutine header_ch_icbc
    implicit none
    integer :: ivarid , idimid , is , istatus

    call getmem2d(poxid_3,1,jx,1,iy,'mod_ch_icbc:poxid_3')
    call getmem4d(chv3,1,jx,1,iy,1,chilev,1,nchsp,'mod_ch_icbc:chv3')

    istatus = nf90_open(trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
                      'mz4_avg_2000-2007_aug.nc', nf90_nowrite, ncid)
    if ( istatus /= nf90_noerr ) then
      write (stderr,*) 'Cannot open input file'
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if

    istatus = nf90_inq_dimid(ncid,'lon',idimid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inquire_dimension(ncid,idimid,len=chilon)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_dimid(ncid,'lat',idimid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inquire_dimension(ncid,idimid,len=chjlat)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_dimid(ncid,'lev',idimid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inquire_dimension(ncid,idimid,len=chilev)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if

    call getmem1d(cht42lon,1,chilon,'mod_ch_icbc:cht42lon')
    call getmem1d(cht42lat,1,chjlat,'mod_ch_icbc:cht42lat')
    call getmem1d(cht42hyam,1,chilev,'mod_ch_icbc:cht42hyam')
    call getmem1d(cht42hybm,1,chilev,'mod_ch_icbc:cht42hybm')
    call getmem2d(xps,1,chilon,1,chjlat,'mod_ch_icbc:xps')
    call getmem2d(poxid_2,1,chilon,1,chjlat,'mod_ch_icbc:poxid_2')
    call getmem4d(xinp,1,chilon,1,chjlat,1,chilev,1,nchsp,'mod_ch_icbc:xinp')

    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,cht42lon)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,cht42lat)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'hyam',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,cht42hyam)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'hybm',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,cht42hybm)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'P0',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,p0)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if

    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,xps)
    if ( istatus /= nf90_noerr ) then
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if

    do is = 1 , nchsp
      istatus = nf90_inq_varid(ncid,trim(chspec(is))//'_VMR_avrg',ivarid)
      if ( istatus /= nf90_noerr ) then
        call die('header_ch_icbc',nf90_strerror(istatus),istatus)
      end if
      istatus = nf90_get_var(ncid,ivarid,xinp(:,:,:,is))
      if ( istatus /= nf90_noerr ) then
        call die('header_ch_icbc',nf90_strerror(istatus),istatus)
      end if
    end do

  end subroutine header_ch_icbc

  subroutine get_ch_icbc(idate)
    implicit none
!
    integer :: i , is , j , k , l , k0
    type(rcm_time_and_date) , intent(in) :: idate
    real(sp) :: r4pt

    do is = 1 , nchsp
      call bilinx2(chv3(:,:,:,is),xinp(:,:,:,is),xlon,xlat,cht42lon,cht42lat, &
                   chilon,chjlat,iy,jx,chilev) 
    end do

    poxid_2 = xps*0.01
    p0 = p0*0.01
    r4pt = real(ptop)

    call bilinx2(poxid_3,poxid_2,xlon,xlat,cht42lon,cht42lat,chilon,chjlat,iy,jx,1)

    do i = 1 , iy 
      do j = 1 , jx
        do l = 1 , kz
          prcm=((poxid_3(j,i)*0.1-r4pt)*sigma2(l)+r4pt)*10.0
          k0 = -1
          do k = chilev , 1 , -1
            pmpi = poxid_3(j,i)*cht42hybm(k)+cht42hyam(k)*p0
            k0 = k
            if (prcm > pmpi) exit
          end do
          if (k0 == chilev) then        
            pmpj = poxid_3(j,i)*cht42hybm(chilev-1)+cht42hyam(chilev-1)*p0
            pmpi = poxid_3(j,i)*cht42hybm(chilev)+cht42hyam(chilev)*p0

            do is = 1 , nchsp
              chv4(j,i,l,is) = chv3(j,i,chilev,is) + &
                 (chv3(j,i,chilev,is) - chv3(j,i,chilev-1,is))*(prcm-pmpi)/(pmpi-pmpj)
            end do
          else if (k0 >= 1) then
            pmpj = poxid_3(j,i)*cht42hybm(k0)+cht42hyam(k0)*p0
            pmpi = poxid_3(j,i)*cht42hybm(k0+1)+cht42hyam(k0+1)*p0
            do is = 1 , nchsp
              chv4(j,i,l,is) = (chv3(j,i,k0+1,is)*(prcm-pmpj) + &
                                chv3(j,i,k0,is)*(prcm-pmpi))/(pmpi-pmpj)
            end do
          end if            
        end do
      end do
    end do            

    call write_ch_icbc(idate)

  end subroutine get_ch_icbc

  subroutine close_ch_icbc
    use netcdf
    implicit none
    istatus=nf90_close(ncid)
    if ( istatus/=nf90_noerr ) then
      write (stderr,*) 'Cannot close input file'
      call die('header_ch_icbc',nf90_strerror(istatus),istatus)
    end if
  end subroutine close_ch_icbc

end module mod_ch_icbc
