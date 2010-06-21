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

      integer :: jx_in , iy_in , kz_in
      real(4) :: ptop_in
      real(4) :: clat_in , clon_in , plat_in , plon_in
      real(4) :: truelath_in , truelatl_in
      character(6) :: cgtype_in
      integer :: igrads_in , ibigend_in

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

      subroutine gridml

      use mod_dynparam
      implicit none
!
! Local variables
!
      integer :: ierr , k
      character(256) :: terfile
!
!     THIS SUBROUTINE CALLS ROUTINES TO PRODUCE THE MAP FACTORS
!     IT ALSO READS A FILE OF TOPOGRAPHY AND LANDUSE APPROPRIATE FOR
!     GRID FOR EXPLANATION OF VARIABLES SEE SUBROUTINE MAPRON.
!
!     READ APPROPRIATE FILE OF TERRAIN AND LANDUSE FOR THIS GRID
!
      write (terfile,99001)                                             &
        & trim(dirter), pthsep, trim(domname), '.INFO'
      open (10,file=terfile,form='unformatted',                         &
          & recl=jx*iy*ibyte,access='direct')
      if ( dattyp=='FVGCM' .or. dattyp=='NRP2W' .or.                    &
         & dattyp=='GFS11' .or. dattyp=='EH5OM' ) then
        read (10,rec=1,iostat=ierr) iy_in , jx_in , kz_in , delx ,      &
                                  & clat_in , clon_in , plat_in ,       &
                                  & plon_in , grdfac , cgtype_in ,      &
                                  & (sigmaf(k),k=1,kz+1) , ptop_in ,    &
                                  & igrads_in , ibigend_in ,            &
                                  & truelatl_in , truelath_in ,         &
                                  & lon0 , lon1 , lat0 , lat1
      else
        read (10,rec=1,iostat=ierr) iy_in , jx_in , kz_in , delx ,      &
                                  & clat_in , clon_in , plat_in ,       &
                                  & plon_in , grdfac , cgtype_in ,      &
                                  & (sigmaf(k),k=1,kz+1) , ptop_in ,    &
                                  & igrads_in , ibigend_in ,            &
                                  & truelatl_in , truelath_in
      end if
      if ( iy_in/=iy .or. jx_in/=jx .or. kz_in/=kz ) then
        print * , 'IMPROPER DIMENSION SPECIFICATION (ICBC.f)'
        print * , '  icbc.param: ' , iy , jx , kz
        print * , '  DOMAIN.INFO: ' , iy_in , jx_in , kz_in
        print * , '  Also check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'Dimensions (subroutine gridml)'
      end if
      read (10,rec=2,iostat=ierr) topogm
      read (10,rec=3,iostat=ierr) toposdgm
      read (10,rec=4,iostat=ierr) xlandu
      read (10,rec=5,iostat=ierr) xlat
      read (10,rec=6,iostat=ierr) xlon
      read (10,rec=7,iostat=ierr) dlat
      read (10,rec=8,iostat=ierr) dlon
      read (10,rec=9,iostat=ierr) msfx
      close (10)
      if ( ierr/=0 ) then
        print * , 'END OF FILE REACHED (ICBC.f)'
        print * , '  Check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'EOF (subroutine gridml)'
      end if
!
99001 format (a,a,a,a)
      end subroutine gridml

      subroutine commhead
      use mod_dynparam
      implicit none
!
! Local variables
!
      integer :: k
!
      call gridml
!
      do k = 1 , kz
        sigma2(k) = 0.5*(sigmaf(k+1)+sigmaf(k))
        dsigma(k) = sigmaf(k+1) - sigmaf(k)
      end do
      end subroutine commhead

      end module mod_grid
