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
  use netcdf

  private
!
  integer :: nyear , month , nday , nhour
  integer :: k , l
  integer :: k0

  integer , parameter :: chilon = 128 , chjlat = 64 , chilev = 28
  real(sp) , dimension(chilon) :: cht42lon
  real(sp) , dimension(chjlat) :: cht42lat
  real(sp) , dimension(chilev) :: cht42hyam , cht42hybm
  real(sp) , dimension(chilon,chjlat) :: xps
!
! Oxidant climatology variables
!
  real(sp) :: p0
  real(sp) , dimension(chilon,chjlat) :: poxid_2
  real(sp) , allocatable, dimension(:,:) :: poxid_3
  real(sp) , allocatable, dimension(:,:,:,:) :: chv3
  real(sp) , dimension(chilon,chjlat,chilev,nchsp) :: xinp

  real(sp) :: prcm , pmpi , pmpj
  integer :: ncid , istatus

  public :: headermozart_ch_icbc , getmozart_ch_icbc , freemozart_ch_icbc

  contains

  subroutine headermozart_ch_icbc
    implicit none
    integer :: ivarid , is , istatus

    allocate(poxid_3(jx,iy))
    allocate(chv3(jx,iy,chilev,nchsp))

    istatus = nf90_open(trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
                      'mz4_avg_2000-2007_aug.nc', nf90_nowrite, ncid)
    if ( istatus /= nf90_noerr ) then
      write (stderr,*) 'Cannot open input file'
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if

    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,cht42lon)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,cht42lat)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'hyam',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,cht42hyam)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'hybm',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,cht42hybm)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(ncid,'P0',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,p0)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if

    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_var(ncid,ivarid,xps)
    if ( istatus /= nf90_noerr ) then
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if

    do is = 1 , nchsp
      print *, chspec(is)//'_VMR_avrg'
      istatus = nf90_inq_varid(ncid,chspec(is)//'_VMR_avrg',ivarid)
      if ( istatus /= nf90_noerr ) then
        call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
      end if
      istatus = nf90_get_var(ncid,ivarid,xinp(:,:,:,is))
      if ( istatus /= nf90_noerr ) then
        call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
      end if
    end do

  end subroutine headermozart_ch_icbc

  subroutine getmozart_ch_icbc(idate)
    implicit none
!
    integer :: i , is , j , k , k0
    type(rcm_time_and_date) , intent(in) :: idate

    do is = 1 , nchsp
      call bilinx2(chv3(:,:,:,is),xinp(:,:,:,is),xlon,xlat,cht42lon,cht42lat, &
                   chilon,chjlat,iy,jx,chilev) 
    end do

    poxid_2 = xps*0.01
    p0 = p0*0.01

    call bilinx2(poxid_3,poxid_2,xlon,xlat,cht42lon,cht42lat,chilon,chjlat,iy,jx,1)

    do i = 1 , iy 
      do j = 1 , jx
        do l = 1 , kz
          prcm=((poxid_3(j,i)*0.1-ptop)*sigma2(l)+ptop)*10.0
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

  end subroutine getmozart_ch_icbc

  subroutine freemozart_ch_icbc
    use netcdf
    implicit none
    istatus=nf90_close(ncid)
    if ( istatus/=nf90_noerr ) then
      write (stderr,*) 'Cannot close input file'
      call die('headermozart_ch_icbc',nf90_strerror(istatus),istatus)
    end if
    deallocate(poxid_3)
    deallocate(chv3)
  end subroutine freemozart_ch_icbc

end module mod_ch_icbc
