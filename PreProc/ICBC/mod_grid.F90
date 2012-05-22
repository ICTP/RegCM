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

module mod_grid

  use mod_memutil
  use mod_realkinds
  use mod_stdio
  use mod_message
  use mod_nchelper

  real(sp) , pointer , dimension(:,:) :: coriol , dlat , dlon , &
         msfx , topogm , xlandu , xlat , xlon
  real(sp) , pointer , dimension(:,:) :: pa , tlayer , za
  real(sp) , pointer , dimension(:,:) :: b3pd
  real(sp) , pointer , dimension(:) :: dsigma , sigma2
  real(sp) , pointer , dimension(:) :: sigmaf
  real(dp) :: delx , grdfac
  integer :: i0 , i1 , j0 , j1
  real(dp) :: lat0 , lat1 , lon0 , lon1

  contains

  subroutine init_grid(iy,jx,kz)
    implicit none
    integer , intent(in) :: iy , jx , kz
    call getmem2d(coriol,1,jx,1,iy,'mod_grid:coriol')
    call getmem2d(dlat,1,jx,1,iy,'mod_grid:dlat')
    call getmem2d(dlon,1,jx,1,iy,'mod_grid:dlon')
    call getmem2d(msfx,1,jx,1,iy,'mod_grid:msfx')
    call getmem2d(topogm,1,jx,1,iy,'mod_grid:topogm')
    call getmem2d(xlandu,1,jx,1,iy,'mod_grid:xlandu')
    call getmem2d(xlat,1,jx,1,iy,'mod_grid:xlat')
    call getmem2d(xlon,1,jx,1,iy,'mod_grid:xlon')
    call getmem2d(pa,1,jx,1,iy,'mod_grid:pa')
    call getmem2d(tlayer,1,jx,1,iy,'mod_grid:tlayer')
    call getmem2d(za,1,jx,1,iy,'mod_grid:za')
    call getmem2d(b3pd,1,jx,1,iy,'mod_grid:b3pd')
    call getmem1d(dsigma,1,kz,'mod_grid:dsigma')
    call getmem1d(sigma2,1,kz,'mod_grid:sigma2')
    call getmem1d(sigmaf,1,kz+1,'mod_grid:sigmaf')

    call read_domain

  end subroutine init_grid

  subroutine read_domain
    use netcdf
    use mod_dynparam
    implicit none
    integer :: istatus
    integer :: incin
    integer :: idimid
    integer :: ivarid
    character(256) :: fname
    integer :: jx_in , iy_in , kz_in , k

    fname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'

    istatus = nf90_open(fname, nf90_nowrite, incin)
    call checkncerr(istatus,__FILE__,__LINE__,'Error open domain file '//trim(fname))

    istatus = nf90_inq_dimid(incin, "iy", idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error searching iy dimension')
    istatus = nf90_inquire_dimension(incin, idimid, len=iy_in)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading iy dimension')
    istatus = nf90_inq_dimid(incin, "jx", idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error searching jx dimension')
    istatus = nf90_inquire_dimension(incin, idimid, len=jx_in)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading jx dimension')
    istatus = nf90_inq_dimid(incin, "kz", idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error searching kz dimension')
    istatus = nf90_inquire_dimension(incin, idimid, len=kz_in)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading kz dimension')

    if ( iy_in /= iy .or. jx_in /= jx .or. kz_in /= kz+1 ) then
      write(stderr,*) 'IMPROPER DIMENSION SPECIFICATION'
      write(stderr,*) '  namelist  : ' , iy , jx , kz
      write(stderr,*) '  DOMAIN    : ' , iy_in , jx_in , kz_in
      call die('read_domain','Dimension mismatch',1)
    end if

    istatus = nf90_get_att(incin, NF90_GLOBAL,"grid_factor",grdfac)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading grid_factor attribute')

    istatus = nf90_inq_varid(incin, "topo", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error searching topo variable')
    istatus = nf90_get_var(incin, ivarid, topogm)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading topo variable')
    istatus = nf90_inq_varid(incin, "landuse", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error searching landuse variable')
    istatus = nf90_get_var(incin, ivarid, xlandu)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading landuse variable')
    istatus = nf90_inq_varid(incin, "xlat", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error searching xlat variable')
    istatus = nf90_get_var(incin, ivarid, xlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading xlat variable')
    istatus = nf90_inq_varid(incin, "xlon", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error searching xlon variable')
    istatus = nf90_get_var(incin, ivarid, xlon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading xlon variable')
    istatus = nf90_inq_varid(incin, "dlat", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error searchin dlat variable')
    istatus = nf90_get_var(incin, ivarid, dlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading dlat variable')
    istatus = nf90_inq_varid(incin, "dlon", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error searching dlon variable')
    istatus = nf90_get_var(incin, ivarid, dlon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading dlon variable')
    istatus = nf90_inq_varid(incin, "xmap", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error searching xmap variable')
    istatus = nf90_get_var(incin, ivarid, msfx)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading xmap variable')
    istatus = nf90_inq_varid(incin, "sigma", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error searching sigma variable')
    istatus = nf90_get_var(incin, ivarid, sigmaf)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading sigma variable')

    istatus = nf90_close(incin)
    call checkncerr(istatus,__FILE__,__LINE__,'Error closing Domain file '//trim(fname))
!
    do k = 1 , kz
      sigma2(k) = 0.5*(sigmaf(k+1)+sigmaf(k))
      dsigma(k) = sigmaf(k+1) - sigmaf(k)
    end do

!
  end subroutine read_domain
!
end module mod_grid
