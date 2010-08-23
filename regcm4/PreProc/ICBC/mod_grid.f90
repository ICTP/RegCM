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

      implicit none

      real(4) , allocatable , dimension(:,:) :: coriol , dlat , dlon ,  &
           & msfx , snowcv , topogm , toposdgm , xlandu , xlat , xlon
      real(4) , allocatable , dimension(:,:) :: pa , sst1 , sst2 ,      &
           & tlayer , za , ice1 , ice2
      real(4) , allocatable , dimension(:,:) :: b3pd
      real(4) , allocatable , dimension(:) :: dsigma , sigma2
      real(4) , allocatable , dimension(:) :: sigmaf
      real(4) :: delx , grdfac
      integer :: i0 , i1 , j0 , j1
      real(4) :: lat0 , lat1 , lon0 , lon1

      contains

      subroutine init_grid(iy,jx,kz)
        implicit none
        integer , intent(in) :: iy , jx , kz
        allocate(coriol(jx,iy))
        allocate(dlat(jx,iy))
        allocate(dlon(jx,iy))
        allocate(msfx(jx,iy))
        allocate(snowcv(jx,iy))
        allocate(topogm(jx,iy))
        allocate(toposdgm(jx,iy))
        allocate(xlandu(jx,iy))
        allocate(xlat(jx,iy))
        allocate(xlon(jx,iy))
        allocate(pa(jx,iy))
        allocate(sst1(jx,iy))
        allocate(sst2(jx,iy))
        allocate(tlayer(jx,iy))
        allocate(za(jx,iy))
        allocate(ice1(jx,iy))
        allocate(ice2(jx,iy))
        allocate(b3pd(jx,iy))
        allocate(dsigma(kz))
        allocate(sigma2(kz))
        allocate(sigmaf(kz+1))

        call read_domain

      end subroutine init_grid

      subroutine free_grid
        implicit none
        deallocate(coriol)
        deallocate(dlat)
        deallocate(dlon)
        deallocate(msfx)
        deallocate(snowcv)
        deallocate(topogm)
        deallocate(toposdgm)
        deallocate(xlandu)
        deallocate(xlat)
        deallocate(xlon)
        deallocate(pa)
        deallocate(sst1)
        deallocate(sst2)
        deallocate(tlayer)
        deallocate(za)
        deallocate(ice1)
        deallocate(ice2)
        deallocate(b3pd)
        deallocate(dsigma)
        deallocate(sigma2)
        deallocate(sigmaf)
      end subroutine free_grid

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
        if ( istatus /= nf90_noerr) then
          write (6,*) 'Error Opening Domain file ', trim(fname)
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        istatus = nf90_inq_dimid(incin, "iy", idimid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Dimension iy missing'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inquire_dimension(incin, idimid, len=iy_in)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error dimension iy'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_dimid(incin, "jx", idimid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Dimension jx missing'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inquire_dimension(incin, idimid, len=jx_in)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error dimension jx'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_dimid(incin, "kz", idimid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Dimension kz missing'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inquire_dimension(incin, idimid, len=kz_in)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error dimension kz'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        if ( iy_in/=iy .or. jx_in/=jx .or. kz_in/=kz+1 ) then
          print * , 'IMPROPER DIMENSION SPECIFICATION'
          print * , '  namelist  : ' , iy , jx , kz
          print * , '  DOMAIN    : ' , iy_in , jx_in , kz_in
          stop 'Dimension mismatch'
        end if

        istatus = nf90_get_att(incin, NF90_GLOBAL,"grid_factor", &
                       &       grdfac)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading grid_factor'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        istatus = nf90_inq_varid(incin, "topo", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error topo variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_get_var(incin, ivarid, topogm)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading topo variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_varid(incin, "htsd", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error htsd variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_get_var(incin, ivarid, toposdgm)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading htsd variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_varid(incin, "landuse", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error landuse variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_get_var(incin, ivarid, xlandu)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading landuse variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_varid(incin, "xlat", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error xlat variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_get_var(incin, ivarid, xlat)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading xlat variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_varid(incin, "xlon", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error xlon variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_get_var(incin, ivarid, xlon)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading xlon variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_varid(incin, "dlat", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error dlat variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_get_var(incin, ivarid, dlat)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading dlat variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_varid(incin, "dlon", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error dlon variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_get_var(incin, ivarid, dlon)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading dlon variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_varid(incin, "xmap", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error xmap variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_get_var(incin, ivarid, msfx)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading xmap variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_varid(incin, "sigma", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error sigma variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_get_var(incin, ivarid, sigmaf)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading sigma variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        istatus = nf90_close(incin)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error closing Domain file ', trim(fname)
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        do k = 1 , kz
          sigma2(k) = 0.5*(sigmaf(k+1)+sigmaf(k))
          dsigma(k) = sigmaf(k+1) - sigmaf(k)
        end do

      end subroutine read_domain

      end module mod_grid
